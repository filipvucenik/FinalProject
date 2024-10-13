import sys
import matplotlib.pyplot as plt

n = 50

def plot_histogram(data_dict, GT_annos, SNP_annos):
    global n
    if n == 0:
        n = len(data_dict)
    keys = list(x.split("_")[0] for x in data_dict.keys())[:n]
    values = list(data_dict.values())[:n]

    plt.figure(figsize=(10, 6))
    
    colors = ['blue' if x in GT_annos and x in SNP_annos and
              len(GT_annos[x].intersection(SNP_annos[x])) >= 0.6*len(GT_annos[x].union(SNP_annos[x]))
              else 'red' 
              for x in list(data_dict.keys())[:n]]  
    
    plt.bar(keys, values, color=colors)
    plt.xlabel('reads')
    plt.ylabel('number of overlaps')
    plt.title('Histogram of SNP Data')
    plt.xticks(rotation=90)  # Rotate x-axis labels if needed
    plt.tight_layout()  # Adjust layout to fit labels
    plt.gcf().set_constrained_layout(True)
    plt.savefig('overlaps_histogram.png')
    plt.show()
    
def plot_scatterplot(data_dict, MQ):
    global n
    if n == 0:
        n = len(data_dict)
    values = list(data_dict.values())[:n]
    mqs = list(MQ[x] if x in MQ
               else 0 
               for x in list(data_dict.keys())[:n])

    plt.figure(figsize=(10, 6))
    plt.scatter(values, mqs)
    plt.xlabel('number of overlaps')
    plt.ylabel('mapping quality')
    plt.title('Scatterplot of Overlaps vs. Mapping Quality')
    plt.tight_layout()  # Adjust layout to fit labels
    plt.gcf().set_constrained_layout(True)
    plt.savefig('overlaps_vs_MQ.png')
    plt.show()

def main():
    if len(sys.argv) < 4:
        print("Usage: python snpHistogram.py <overlaps_file> <SNP_annotation_file> <GT_annotation_file>")
        sys.exit(1)
    
    overlaps_file = sys.argv[1]    
    SNP_annotation_file = sys.argv[2]
    GT_annotation_file = sys.argv[3]
    overlapsPerRead = {}
    
    with open(overlaps_file) as f:
        for line in f:
            line = line.strip()
            if line:
                read = line.split('\t')[0]
                if read in overlapsPerRead:
                    overlapsPerRead[read] += 1
                else:
                    overlapsPerRead[read] = 1
    
    SNP_annos = {}
    MQ = {}
    
    with open(SNP_annotation_file) as f:
        for line in f:
            line = line.strip()
            if line:
                split = line.split(" ")
                read = split[0]
                annos = split[1:-1]
                mq = float(split[-1].split(":")[1])
                SNP_annos[read] = set(x.split(":")[0] for x in annos)
                MQ [read] = mq 
                
    GT_annos = {}
    
    with open(GT_annotation_file) as f:
        for line in f:
            line = line.strip()
            if line:
                split = line.split(" ")
                read = split[0]
                annos = split[1:-1]
                GT_annos[read] = set(x.split(":")[0] for x in annos)
                    

                    
    plot_histogram(overlapsPerRead, GT_annos, SNP_annos)
    plot_scatterplot(overlapsPerRead, MQ)
    
    
if __name__ == "__main__":
    main()