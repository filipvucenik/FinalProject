import pandas as pd
import numpy as np
import os
import argparse
import pysam


def calcualte_if_overlaping(b1, b2, e1, e2):
    return max(b1, b2) < min(e1, e2)


def get_nanosim_read_info(name):
    chromosome = name.split('_')[0]
    chromosome = chromosome.replace("-", "_")

    start = int(name.split('_')[1])
    strand = name.split('_')[4]

    return chromosome, start, strand


def parse_bed(bed_path):
    bed_pd = pd.read_csv(bed_path, sep='\t', header=None,
                         names=["chromosome1", "start1", "end1", "strand", "chromosome2", "start2", "end2"])
    return bed_pd


def get_raft_read_info(name):
    tmp = name.split(',')
    chromosome, start_og, strand = get_nanosim_read_info(tmp[1])
    start_relative = int(tmp[2].split("=")[1].split("-")[0])

    return chromosome, start_og + start_relative, strand


def get_seqread_read_info(name):
    tmp = name.split(',')
    convert = {"forward": "+", "reverse": "-"}
    strand = convert[tmp[1]]
    chromosome = tmp[4]
    start = int(tmp[2].split("=")[1].split("-")[0])
    return chromosome, start, strand


def parse_reads(reads_path, read_source):
    list_of_reads = []
    with pysam.FastxFile(reads_path) as reads:
        for entry in reads:
            if read_source == "nanosim":
                chromosome, start, strand = get_nanosim_read_info(entry.name)
            elif read_source == "RAFT":
                chromosome, start, strand = get_raft_read_info(entry.name)
            elif read_source == "seqread":
                chromosome, start, strand = get_seqread_read_info(entry.name)

            list_of_reads.append([entry.name, chromosome, start, start + len(entry.sequence), strand])

    return pd.DataFrame(list_of_reads, columns=["read_name", "chromosome", "start", "end", "strand"])


def parse_paf(file_path, source="raven"):
    if source == "raven":
        names = ["qn", "ql", "qs", "qe", "strand", "tn", "tl", "ts", "te", "match", "aln_len", "mapq", "cigar",
                 "editdist"]

    return pd.read_csv(file_path, sep='\t', header=None, names=names)


def overlap_type(b1, b2, e1, e2):
    if b1 <= b2 and e1 >= e2:
        return 0  # 2 fully contained in 1
        # return e2-b2
    elif b1 >= b2 and e1 <= e2:
        return 1  # 1 fully contained in 2
        # return e1-b1
    elif b1 > b2 and e1 > e2:
        return 2  # 2 -> 1
        # return e2-b1
    elif b1 < b2 and e1 < e2:
        return 3  # 1 -> 2
        # return e1-b2


def calc_overlap_lenght(b1, b2, e1, e2):
    if b1 <= b2 and e1 >= e2:
        # return 0, e2-b2 # 2 fully contained in 1
        return e2 - b2
    elif b1 >= b2 and e1 <= e2:
        # return 1, e1-e2 # 1 fully contained in 2
        return e1 - b1
    elif b1 > b2 and e1 > e2:
        # return 2, e2-b1 # 2 -> 1
        return e2 - b1
    elif b1 < b2 and e1 < e2:
        # return 3, e1-b2 # 1 -> 2
        return e1 - b2


def get_gt_from_read_dict(read_dict, homozygous_areas, min_overlap):
    per_read_overlaps = []
    for index, row in read_dict.iterrows():
        print(f"Processing row {index} of {len(read_dict)}")
        # Filter out reads that are the same chromosome as the current one
        read_dict_tmp = read_dict.loc[read_dict['chromosome'] == row['chromosome']].reset_index(drop=True)
        read_dict_tmp = read_dict_tmp.loc[read_dict_tmp['read_name'] != row['read_name']].reset_index(drop=True)
        overlapping = read_dict_tmp.apply(
            lambda x: calcualte_if_overlaping(row['start'], x['start'], row['end'], x['end']), axis=1)
        overlapping_reads = read_dict_tmp[overlapping].reset_index(drop=True)
        ol_types = overlapping_reads.apply(lambda x: overlap_type(row['start'], x['start'], row['end'], x['end']),
                                           axis=1)
        overlap_lenghts = overlapping_reads.apply(
            lambda x: calc_overlap_lenght(row['start'], x['start'], row['end'], x['end']), axis=1)
        # long_overlaps = (overlap_lenghts > (row["end"] - row["start"])*min_overlap).astype(bool)
        long_overlaps = (overlap_lenghts > 500).astype(bool)
        overlap_lenghts = overlap_lenghts[long_overlaps].reset_index(drop=True).values
        ol_types = ol_types[long_overlaps].reset_index(drop=True).values
        same_ht_overlaps = overlapping_reads.loc[long_overlaps, "read_name"].reset_index(drop=True).values

        for overlapping_read, overlap_lenght, ol_type in zip(same_ht_overlaps, overlap_lenghts, ol_types):
            per_read_overlaps.append([row['read_name'], overlapping_read, overlap_lenght, ol_type, "same"])
            # Find overlaps with other haplotype
        other_haplotype = homozygous_areas.loc[homozygous_areas["chromosome1"] == row["chromosome"]].reset_index(
            drop=True)
        other_ht_overlapping = other_haplotype.apply(
            lambda x: calcualte_if_overlaping(row['start'], x['start1'], row['end'], x['end1']), axis=1)
        other_haplotype = other_haplotype[other_ht_overlapping].reset_index(drop=True)
        # greater_than_min_overlap = other_haplotype.apply(lambda x: calc_overlap_lenght(row['start'], x['start1'], row['end'], x['end1']), axis = 1) > (row["end"] - row["start"])*min_overlap
        greater_than_min_overlap = other_haplotype.apply(
            lambda x: calc_overlap_lenght(row['start'], x['start1'], row['end'], x['end1']), axis=1) > 500
        other_haplotype = other_haplotype[greater_than_min_overlap].reset_index(drop=True)
        if len(other_haplotype) != 0:
            print("Found other haplotype overlaps")
            for index2, row2 in other_haplotype.iterrows():
                read_dict_tmp2 = read_dict.loc[read_dict['chromosome'] == row2['chromosome2']].reset_index(drop=True)
                overlapping = read_dict_tmp2.apply(
                    lambda x: calcualte_if_overlaping(row2['start2'], x['start'], row2['end2'], x['end']), axis=1)
                overlapping_reads2 = read_dict_tmp2[overlapping].reset_index(drop=True)
                ol_types2 = overlapping_reads2.apply(
                    lambda x: overlap_type(row2['start2'], x['start'], row2['end2'], x['end']), axis=1)
                overlap_lenghts2 = overlapping_reads2.apply(
                    lambda x: calc_overlap_lenght(row2['start2'], x['start'], row2['end2'], x['end']), axis=1)
                # long_overlaps2 = (overlap_lenghts2 > (row["end"] - row["start"])*min_overlap).astype(bool)
                long_overlaps2 = (overlap_lenghts2 > 500).astype(bool)
                overlap_lenghts2 = overlap_lenghts2[long_overlaps2].reset_index(drop=True).values
                ol_types2 = ol_types2[long_overlaps2].reset_index(drop=True).values

                overlapping_reads2 = overlapping_reads2.loc[long_overlaps2, "read_name"].reset_index(drop=True).values

                for overlapping_read, overlap_lenght, ol_type in zip(overlapping_reads2, overlap_lenghts2, ol_types2):
                    per_read_overlaps.append([row['read_name'], overlapping_read, overlap_lenght, ol_type, "other"])

    per_read_overlaps = pd.DataFrame(per_read_overlaps,
                                     columns=["read_name", "overlapping_read", "overlap_lenght", "overlap_type",
                                              "strand"])
    return per_read_overlaps


def calc_metrics(gt_tmp, paf_tmp):
    overlapping_gt = gt_tmp["overlapping_read"].values
    n_overlapping_gt = len(overlapping_gt)
    overlapping_reads = paf_tmp["tn"].values
    n_overlapping_paf = len(overlapping_reads)
    true_positive = len(set(overlapping_reads).intersection(set(overlapping_gt)))
    false_positive = n_overlapping_paf - true_positive
    false_negative = n_overlapping_gt - true_positive

    return n_overlapping_gt, n_overlapping_paf, true_positive, false_positive, false_negative


def calc_stats(gt_pd, paf_pd):
    n_rows = len(gt_pd) - 1
    stats = []
    stats_pd = pd.DataFrame([], columns=["read", "true_positive", "false_positive", "false_negative"])
    count = 0
    for read in gt_pd['read_name'].unique():
        print(f"Processing read {count} out of {n_rows}")
        gt_tmp = gt_pd.loc[gt_pd["read_name"] == read].reset_index(drop=True)
        paf_tmp = paf_pd.loc[paf_pd["qn"] == read].reset_index(drop=True)

        n_overlapping_gt, n_overlapping_paf, true_positive, false_positive, false_negative = calc_metrics(gt_tmp,
                                                                                                          paf_tmp)

        n_overlapping_gt_same, n_overlapping_paf_same, true_positive_same, false_positive_same, false_negative_same = calc_metrics(
            gt_tmp[gt_tmp.strand == "same"], paf_tmp)
        n_overlapping_gt_other, n_overlapping_paf_other, true_positive_other, false_positive_other, false_negative_other = calc_metrics(
            gt_tmp[gt_tmp.strand == "other"], paf_tmp)

        stats.append({"read": read, "true_overlaps": n_overlapping_gt, "found_overlaps": n_overlapping_paf,
                      "true_positive": true_positive, "false_positive": false_positive,
                      "false_negative": false_negative,
                      "true_overlaps_same": n_overlapping_gt_same, "found_overlaps_same": n_overlapping_paf_same,
                      "true_positive_same": true_positive_same, "false_positive_same": false_positive_same,
                      "false_negative_same": false_negative_same,
                      "true_overlaps_other": n_overlapping_gt_other, "found_overlaps_other": n_overlapping_paf_other,
                      "true_positive_other": true_positive_other, "false_positive_other": false_positive_other,
                      "false_negative_other": false_negative_other})
        count += 1
    stats_pd = pd.DataFrame(stats)

    return stats_pd


def calc_global_stats(stats_pd):
    true_positive = stats_pd["true_positive"].sum()
    false_positive = stats_pd["false_positive"].sum()
    false_negative = stats_pd["false_negative"].sum()

    precision = true_positive / (true_positive + false_positive)
    recall = true_positive / (true_positive + false_negative)
    f1 = 2 * (precision * recall) / (precision + recall)

    true_positive_same = stats_pd["true_positive_same"].sum()
    false_positive_same = stats_pd["false_positive_same"].sum()
    false_negative_same = stats_pd["false_negative_same"].sum()

    precision_same = true_positive_same / (true_positive_same + false_positive_same)
    recall_same = true_positive_same / (true_positive_same + false_negative_same)
    f1_same = 2 * (precision_same * recall_same) / (precision_same + recall_same)

    true_positive_other = stats_pd["true_positive_other"].sum()
    false_positive_other = stats_pd["false_positive_other"].sum()
    false_negative_other = stats_pd["false_negative_other"].sum()

    precision_other = true_positive_other / (true_positive_other + false_positive_other)
    recall_other = true_positive_other / (true_positive_other + false_negative_other)
    f1_other = 2 * (precision_other * recall_other) / (precision_other + recall_other)

    return pd.DataFrame({"precision": precision, "recall": recall, "f1": f1,
                         "precision_same": precision_same, "recall_same": recall_same, "f1_same": f1_same,
                         "precision_other": precision_other, "recall_other": recall_other, "f1_other": f1_other},
                        index=[0])


def main(mode, reads_soruce, raw_reads_path, read_mapping, mapping_soruce, homozygous_areas, min_overlap, output,
         input=None):
    if (mode == 'gt' or mode == 'full'):
        print("Loading BED")
        homozygous_areas_pd = parse_bed(homozygous_areas)
        print("Parsing reads...")
        read_pd = parse_reads(raw_reads_path, reads_soruce)
        print("Creating ground truth...")
        gt_pd = get_gt_from_read_dict(read_pd, homozygous_areas_pd, min_overlap)
        print("Saving ground truth data...")
        gt_pd.to_csv(f"{output}.gt_reads.csv", index=False)
    if (mode == 'analyze' or mode == 'full'):
        if mode == 'analyze':
            print("Loading ground truth...")
            gt_pd = pd.read_csv(input, header=0)
        print("Parsing mapping...")
        paf_pd = parse_paf(read_mapping, mapping_soruce)
        print("Calculating statistics...")
        stats_pd = calc_stats(gt_pd, paf_pd)
        print("Saving per read statistics...")
        stats_pd.to_csv(f"{output}.stats.csv", index=False)
        print("Calculating global statistics...")
        global_stats_pd = calc_global_stats(stats_pd)
        print("Saving global statistics...")
        global_stats_pd.to_csv(f"{output}.global_stats.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Overlap GT construct')
    parser.add_argument('--mode', choices=["gt", "analyze", "full"], type=str, required=True)
    parser.add_argument('--read_source', choices=["nanosim", "RAFT", "seqread"], type=str)
    parser.add_argument('--read_mapping', type=str)
    parser.add_argument('--raw_reads', type=str)
    parser.add_argument('--mapping_source', choices=["raven", "hifiasm"], type=str)
    parser.add_argument('--output', type=str)
    parser.add_argument('--min_overlap', type=float, default=0.05)
    parser.add_argument('--homozygous_areas', type=str, default=None)
    parser.add_argument('--input', type=str, default=None)
    args = parser.parse_args()

    if args.mode == "full":
        main(args.mode, args.read_source, args.raw_reads, args.read_mapping, args.mapping_source, args.homozygous_areas,
             args.min_overlap, args.output, None)
    elif args.mode == "gt":
        main(args.mode, args.read_source, args.raw_reads, None, None, args.homozygous_areas, args.min_overlap,
             args.output, None)
    elif args.mode == "analyze":
        main(args.mode, None, None, args.read_mapping, args.mapping_source, None, args.min_overlap, args.output,
             args.input)
