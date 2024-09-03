#include "overlap_source.hpp"
#include "bioparser/fasta_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <filesystem>

#include <stdint.h>

#define KSW_NEG_INF -0x40000000

#define KSW_EZ_SCORE_ONLY 0x01  // don't record alignment path/cigar
#define KSW_EZ_RIGHT 0x02       // right-align gaps
#define KSW_EZ_GENERIC_SC 0x04  // without this flag: match/mismatch only; last symbol is a wildcard
#define KSW_EZ_APPROX_MAX 0x08  // approximate max; this is faster with sse
#define KSW_EZ_APPROX_DROP 0x10 // approximate Z-drop; faster with sse
#define KSW_EZ_EXTZ_ONLY 0x40   // only perform extension
#define KSW_EZ_REV_CIGAR 0x80   // reverse CIGAR in the output
#define KSW_EZ_SPLICE_FOR 0x100
#define KSW_EZ_SPLICE_REV 0x200
#define KSW_EZ_SPLICE_FLANK 0x400

#ifdef __cplusplus
extern "C"
{
#endif

    typedef struct
    {
        uint32_t max : 31, zdropped : 1;
        int max_q, max_t; // max extension coordinate
        int mqe, mqe_t;   // max score when reaching the end of query
        int mte, mte_q;   // max score when reaching the end of target
        int score;        // max score reaching both ends; may be KSW_NEG_INF
        int m_cigar, n_cigar;
        int reach_end;
        uint32_t *cigar;
    } ksw_extz_t;

    /**
     * NW-like extension
     *
     * @param km        memory pool, when used with kalloc
     * @param qlen      query length
     * @param query     query sequence with 0 <= query[i] < m
     * @param tlen      target length
     * @param target    target sequence with 0 <= target[i] < m
     * @param m         number of residue types
     * @param mat       m*m scoring mattrix in one-dimension array
     * @param gapo      gap open penalty; a gap of length l cost "-(gapo+l*gape)"
     * @param gape      gap extension penalty
     * @param w         band width (<0 to disable)
     * @param zdrop     off-diagonal drop-off to stop extension (positive; <0 to disable)
     * @param flag      flag (see KSW_EZ_* macros)
     * @param ez        (out) scores and cigar
     */
    void ksw_extz(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                  int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);

    void ksw_extz2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                       int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

    void ksw_extd(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                  int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2, int w, int zdrop, int flag, ksw_extz_t *ez);

    void ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                       int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

    void ksw_exts2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                       int8_t gapo, int8_t gape, int8_t gapo2, int8_t noncan, int zdrop, int flag, ksw_extz_t *ez);

    void ksw_extf2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t mch, int8_t mis, int8_t e, int w, int xdrop, ksw_extz_t *ez);

    /**
     * Global alignment
     *
     * (first 10 parameters identical to ksw_extz_sse())
     * @param m_cigar   (modified) max CIGAR length; feed 0 if cigar==0
     * @param n_cigar   (out) number of CIGAR elements
     * @param cigar     (out) BAM-encoded CIGAR; caller need to deallocate with kfree(km, )
     *
     * @return          score of the alignment
     */
    int ksw_gg(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);
    int ksw_gg2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);
    int ksw_gg2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

    void *ksw_ll_qinit(void *km, int size, int qlen, const uint8_t *query, int m, const int8_t *mat);
    int ksw_ll_i16(void *q, int tlen, const uint8_t *target, int gapo, int gape, int *qe, int *te);

#ifdef __cplusplus
}
#endif

/************************************
 *** Private macros and functions ***
 ************************************/

#ifdef HAVE_KALLOC
#include "kalloc.h"
#else
#include <stdlib.h>
#define kmalloc(km, size) malloc((size))
#define kcalloc(km, count, size) calloc((count), (size))
#define krealloc(km, ptr, size) realloc((ptr), (size))
#define kfree(km, ptr) free((ptr))
#endif

static inline uint32_t *ksw_push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len)
{
    if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1] & 0xf))
    {
        if (*n_cigar == *m_cigar)
        {
            *m_cigar = *m_cigar ? (*m_cigar) << 1 : 4;
            cigar = (uint32_t *)krealloc(km, cigar, (*m_cigar) << 2);
        }
        cigar[(*n_cigar)++] = len << 4 | op;
    }
    else
        cigar[(*n_cigar) - 1] += len << 4;
    return cigar;
}

// In the backtrack matrix, value p[] has the following structure:
//   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F, 3 for \tilde{E} and 4 for \tilde{F}
//   bit 3/0x08: 1 if a continuation on the E state (bit 5/0x20 for a continuation on \tilde{E})
//   bit 4/0x10: 1 if a continuation on the F state (bit 6/0x40 for a continuation on \tilde{F})
static inline void ksw_backtrack(void *km, int is_rot, int is_rev, int min_intron_len, const uint8_t *p, const int *off, const int *off_end, int n_col, int i0, int j0,
                                 int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{ // p[] - lower 3 bits: which type gets the max; bit
    int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
    uint32_t *cigar = *cigar_, tmp;
    while (i >= 0 && j >= 0)
    { // at the beginning of the loop, _state_ tells us which state to check
        int force_state = -1;
        if (is_rot)
        {
            r = i + j;
            if (i < off[r])
                force_state = 2;
            if (off_end && i > off_end[r])
                force_state = 1;
            tmp = force_state < 0 ? p[(size_t)r * n_col + i - off[r]] : 0;
        }
        else
        {
            if (j < off[i])
                force_state = 2;
            if (off_end && j > off_end[i])
                force_state = 1;
            tmp = force_state < 0 ? p[(size_t)i * n_col + j - off[i]] : 0;
        }
        if (state == 0)
            state = tmp & 7; // if requesting the H state, find state one maximizes it.
        else if (!(tmp >> (state + 2) & 1))
            state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
        if (state == 0)
            state = tmp & 7; // TODO: probably this line can be merged into the "else if" line right above; not 100% sure
        if (force_state >= 0)
            state = force_state;
        if (state == 0)
            cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
        else if (state == 1 || (state == 3 && min_intron_len <= 0))
            cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i; // deletion
        else if (state == 3 && min_intron_len > 0)
            cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 3, 1), --i; // intron
        else
            cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j; // insertion
    }
    if (i >= 0)
        cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, min_intron_len > 0 && i >= min_intron_len ? 3 : 2, i + 1); // first deletion
    if (j >= 0)
        cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
    if (!is_rev)
        for (i = 0; i < n_cigar >> 1; ++i) // reverse CIGAR
            tmp = cigar[i], cigar[i] = cigar[n_cigar - 1 - i], cigar[n_cigar - 1 - i] = tmp;
    *m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}

static inline void ksw_reset_extz(ksw_extz_t *ez)
{
    ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
    ez->max = 0, ez->score = ez->mqe = ez->mte = KSW_NEG_INF;
    ez->n_cigar = 0, ez->zdropped = 0, ez->reach_end = 0;
}

static inline int ksw_apply_zdrop(ksw_extz_t *ez, int is_rot, int32_t H, int a, int b, int zdrop, int8_t e)
{
    int r, t;
    if (is_rot)
        r = a, t = b;
    else
        r = a + b, t = a;
    if (H > (int32_t)ez->max)
    {
        ez->max = H, ez->max_t = t, ez->max_q = r - t;
    }
    else if (t >= ez->max_t && r - t >= ez->max_q)
    {
        int tl = t - ez->max_t, ql = (r - t) - ez->max_q, l;
        l = tl > ql ? tl - ql : ql - tl;
        if (zdrop >= 0 && ez->max - H > zdrop + l * e)
        {
            ez->zdropped = 1;
            return 1;
        }
    }
    return 0;
}

/*






*/
// ksw2_gg2_sse

#ifdef __SSE2__
#include <emmintrin.h>

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

int ksw_gg2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
    int r, t, n_col, n_col_, *off, tlen_, last_st, last_en, H0 = 0, last_H0_t = 0;
    uint8_t *qr, *mem, *mem2;
    __m128i *u, *v, *x, *y, *s, *p;
    __m128i q_, qe2_, zero_, flag1_, flag2_, flag8_, flag16_;

    zero_ = _mm_set1_epi8(0);
    q_ = _mm_set1_epi8(q);
    qe2_ = _mm_set1_epi8((q + e) * 2);
    flag1_ = _mm_set1_epi8(1);
    flag2_ = _mm_set1_epi8(2);
    flag8_ = _mm_set1_epi8(0x08);
    flag16_ = _mm_set1_epi8(0x10);

    if (w < 0)
        w = tlen > qlen ? tlen : qlen;
    n_col = w + 1 < tlen ? w + 1 : tlen; // number of columns in the backtrack matrix
    tlen_ = (tlen + 15) / 16;
    n_col_ = (n_col + 15) / 16 + 1;
    n_col = n_col_ * 16;

    mem = (uint8_t *)kcalloc(km, tlen_ * 5 + 1, 16);
    u = (__m128i *)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
    v = u + tlen_, x = v + tlen_, y = x + tlen_, s = y + tlen_;
    qr = (uint8_t *)kcalloc(km, qlen, 1);
    mem2 = (uint8_t *)kmalloc(km, ((qlen + tlen - 1) * n_col_ + 1) * 16);
    p = (__m128i *)(((size_t)mem2 + 15) >> 4 << 4);
    off = (int *)kmalloc(km, (qlen + tlen - 1) * sizeof(int));

    for (t = 0; t < qlen; ++t)
        qr[t] = query[qlen - 1 - t];

    for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r)
    {
        int st = 0, en = tlen - 1, st0, en0, st_, en_;
        int8_t x1, v1;
        __m128i x1_, v1_, *pr;
        // find the boundaries
        if (st < r - qlen + 1)
            st = r - qlen + 1;
        if (en > r)
            en = r;
        if (st < (r - w + 1) >> 1)
            st = (r - w + 1) >> 1; // take the ceil
        if (en > (r + w) >> 1)
            en = (r + w) >> 1; // take the floor
        st0 = st, en0 = en;
        st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
        off[r] = st;
        // set boundary conditions
        if (st > 0)
        {
            if (st - 1 >= last_st && st - 1 <= last_en)
                x1 = ((uint8_t *)x)[st - 1], v1 = ((uint8_t *)v)[st - 1]; // (r-1,s-1) calculated in the last round
            else
                x1 = v1 = 0; // not calculated; set to zeros
        }
        else
            x1 = 0, v1 = r ? q : 0;
        if (en >= r)
            ((uint8_t *)y)[r] = 0, ((uint8_t *)u)[r] = r ? q : 0;
        // loop fission: set scores first
        for (t = st0; t <= en0; ++t)
            ((uint8_t *)s)[t] = mat[target[t] * m + qr[t + qlen - 1 - r]];
        // core loop
        x1_ = _mm_cvtsi32_si128(x1);
        v1_ = _mm_cvtsi32_si128(v1);
        st_ = st >> 4, en_ = en >> 4;
        pr = p + r * n_col_ - st_;
        for (t = st_; t <= en_; ++t)
        {
            __m128i d, z, a, b, xt1, vt1, ut, tmp;

            z = _mm_add_epi8(_mm_load_si128(&s[t]), qe2_);

            xt1 = _mm_load_si128(&x[t]);                     // xt1 <- x[r-1][t..t+15]
            tmp = _mm_srli_si128(xt1, 15);                   // tmp <- x[r-1][t+15]
            xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); // xt1 <- x[r-1][t-1..t+14]
            x1_ = tmp;
            vt1 = _mm_load_si128(&v[t]);                     // vt1 <- v[r-1][t..t+15]
            tmp = _mm_srli_si128(vt1, 15);                   // tmp <- v[r-1][t+15]
            vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); // vt1 <- v[r-1][t-1..t+14]
            v1_ = tmp;
            a = _mm_add_epi8(xt1, vt1); // a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14]

            ut = _mm_load_si128(&u[t]);                  // ut <- u[t..t+15]
            b = _mm_add_epi8(_mm_load_si128(&y[t]), ut); // b <- y[r-1][t..t+15] + u[r-1][t..t+15]

            d = _mm_and_si128(_mm_cmpgt_epi8(a, z), flag1_); // d = a > z? 1 : 0
#ifdef __SSE4_1__
            z = _mm_max_epi8(z, a); // z = z > a? z : a (signed)
            tmp = _mm_cmpgt_epi8(b, z);
            d = _mm_blendv_epi8(d, flag2_, tmp); // d = b > z? 2 : d
#else                                            // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
            z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_)); // z = z > 0? z : 0;
            z = _mm_max_epu8(z, a);                         // z = max(z, a); this works because both are non-negative
            tmp = _mm_cmpgt_epi8(b, z);
            d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, flag2_)); // d = b > z? 2 : d; emulating blendv
#endif
            z = _mm_max_epu8(z, b);                       // z = max(z, b); this works because both are non-negative
            _mm_store_si128(&u[t], _mm_sub_epi8(z, vt1)); // u[r][t..t+15] <- z - v[r-1][t-1..t+14]
            _mm_store_si128(&v[t], _mm_sub_epi8(z, ut));  // v[r][t..t+15] <- z - u[r-1][t..t+15]

            z = _mm_sub_epi8(z, q_);
            a = _mm_sub_epi8(a, z);
            b = _mm_sub_epi8(b, z);
            tmp = _mm_cmpgt_epi8(a, zero_);
            d = _mm_or_si128(d, _mm_and_si128(flag8_, tmp));
            _mm_store_si128(&x[t], _mm_and_si128(a, tmp));
            tmp = _mm_cmpgt_epi8(b, zero_);
            d = _mm_or_si128(d, _mm_and_si128(flag16_, tmp));
            _mm_store_si128(&y[t], _mm_and_si128(b, tmp));
            _mm_store_si128(&pr[t], d);
        }
        if (r > 0)
        {
            if (last_H0_t >= st0 && last_H0_t <= en0)
                H0 += ((uint8_t *)v)[last_H0_t] - (q + e);
            else
                ++last_H0_t, H0 += ((uint8_t *)u)[last_H0_t] - (q + e);
        }
        else
            H0 = ((uint8_t *)v)[0] - 2 * (q + e), last_H0_t = 0;
        last_st = st, last_en = en;
        // for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%x\n", r, t, ((uint8_t*)u)[t], ((uint8_t*)v)[t], ((uint8_t*)x)[t], ((uint8_t*)y)[t], ((uint8_t*)(p + r * n_col_))[t-st]); // for debugging
    }
    kfree(km, mem);
    kfree(km, qr);
    ksw_backtrack(km, 1, 0, 0, (uint8_t *)p, off, 0, n_col, tlen - 1, qlen - 1, m_cigar_, n_cigar_, cigar_);
    kfree(km, mem2);
    kfree(km, off);
    return H0;
}
#endif // __SSE2__
std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

class PafOverlapSource : public OverlapSource
{
private:
    std::vector<std::vector<biosoup::Overlap>> overlaps;
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    std::string reads;
    std::string paf_file;
    std::unordered_map<std::string, std::size_t> sequence_id;
    std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> create_parser(std::string &path)
    {
        auto is_suffix = [](const std::string &s, const std::string &suff)
        {
            return s.size() >= suff.size() && s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
        };

        if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
            is_suffix(path, ".fna") || is_suffix(path, ".fna.gz") ||
            is_suffix(path, ".fa") || is_suffix(path, ".fa.gz"))
        {
            try
            {
                return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);
            }
            catch (const std::invalid_argument &exception)
            {
                std::cerr << exception.what() << std::endl;
            }
        }
        return nullptr;
    }

    void load_sequences()
    {
        auto parser = create_parser(reads);
        while (true)
        {
            std::vector<std::unique_ptr<biosoup::NucleicAcid>> buffer;
            try
            {
                buffer = parser->Parse(1U << 30);
            }
            catch (std::invalid_argument &exception)
            {
                std::cerr << exception.what() << std::endl;
                exit(1);
            }

            if (buffer.empty())
            {
                break;
            }
            sequences.reserve(sequences.size() + buffer.size());
            for (const auto &sequence : buffer)
            {
                sequences.push_back(std::make_unique<biosoup::NucleicAcid>(*sequence));
                sequence_id[sequence->name] = sequences.size() - 1;
            }
        }
    }

    void load_paf()
    {
        std::string line;
        std::ifstream fileStream(paf_file);
        std::cout << "Loading paf file" << std::endl;
        std::string tmp;

        if (!fileStream.is_open())
        {
            std::cerr << "Error opening file" << std::endl;
            std::filesystem::path currentPath = std::filesystem::current_path();
            std::cout << "Current directory: " << currentPath << std::endl;
            std::cout << "File path: " << paf_file << std::endl;
            exit(1);
        }

        while (std::getline(fileStream, line))
        {
            std::istringstream iss(line);
            std::string v;
            std::vector<std::string> variables;
            while (std::getline(iss, v, '\t'))
            {
                variables.push_back(v);
            }

            std::size_t id_l = sequence_id[variables[0]];
            std::size_t id_r = sequence_id[variables[5]];

            overlaps[id_l].emplace_back(id_l,
                                        std::stoi(variables[2]),
                                        std::stoi(variables[3]),
                                        id_r,
                                        std::stoi(variables[7]),
                                        std::stoi(variables[8]),
                                        255,
                                        tmp,
                                        variables[4] == "+");
        }

        std::cout << "Loaded paf file" << std::endl;

        auto ksw2_wrapper = [&](
                                std::uint32_t i,
                                const biosoup::Overlap &it,
                                const std::string &lhs,
                                const std::string &rhs) -> std::string
        {
            std::int8_t m = 3;
            std::int8_t n = -5;
            std::int8_t g = 4;
            std::int8_t e = 4;
            std::int8_t mn[25] = {
                m, n, n, n, 0,
                n, m, n, n, 0,
                n, n, m, n, 0,
                n, n, n, m, 0,
                0, 0, 0, 0, 0};

            std::unordered_map<char, std::uint8_t> transform = {
                {'A', 0}, {'a', 0}, {'C', 1}, {'c', 1}, {'G', 2}, {'g', 2}, {'T', 3}, {'t', 3}};

            auto lhs_ = new std::uint8_t[lhs.size()];
            for (std::size_t j = 0; j < lhs.size(); ++j)
            {
                lhs_[j] = transform[lhs[j]];
            }

            auto rhs_ = new std::uint8_t[rhs.size()];
            for (std::size_t j = 0; j < rhs.size(); ++j)
            {
                rhs_[j] = transform[rhs[j]];
            }

            int m_cigar = 0, n_cigar = 0;
            std::uint32_t *cigar = nullptr;

            auto score = ksw_gg2_sse(
                nullptr,    // void *km
                lhs.size(), // int qlen
                lhs_,       // const uint8_t *query
                rhs.size(), // int tlen
                rhs_,       // const uint8_t *target
                5,          // int8_t m
                mn,         // const int8_t *mat
                g,          // int8_t gapo
                e,          // int8_t gape
                500,        // int w
                &m_cigar,   // int *m_cigar_
                &n_cigar,   // int *n_cigar_
                &cigar);    // uint32_t **cigar_
            std::string cigar_string = "";
            if (n_cigar > 0)
            {

                std::size_t lhs_i = 0;
                std::size_t rhs_i = 0;
                for (std::size_t j = 0; j < static_cast<std::size_t>(n_cigar); ++j)
                {
                    std::size_t count = cigar[j] >> 4;
                    std::size_t op = cigar[j] & 15;
                    switch (op)
                    {
                    case 0:
                    { // M
                        bool match = true;
                        std::size_t inarw_count = 0;
                        for (std::size_t k = 0; k < count; k++)
                        {
                            if (lhs[lhs_i + k] != rhs[rhs_i + k] && match)
                            {
                                cigar_string += std::to_string(inarw_count) + "=";
                                inarw_count = 0;
                                match = false;
                            }
                            else if (lhs[lhs_i + k] == rhs[rhs_i + k] && !match)
                            {
                                cigar_string += std::to_string(inarw_count) + "X";
                                inarw_count = 0;
                                match = true;
                            }
                            inarw_count++;
                        }
                        if (match)
                            cigar_string += std::to_string(inarw_count) + "=";
                        else
                            cigar_string += std::to_string(inarw_count) + "X";
                        lhs_i += count;
                        rhs_i += count;
                        break;
                    }
                    case 1:
                    { // I
                        cigar_string += std::to_string(count) + "I";
                        lhs_i += count;
                        break;
                    }
                    case 2:
                    { // D
                        cigar_string += std::to_string(count) + "D";
                        rhs_i += count;
                        break;
                    }
                    default:
                        break;
                    }
                }
            }
            free(cigar);
            delete[] rhs_;
            delete[] lhs_;

            return cigar_string;
        };

        auto threads = std::make_shared<thread_pool::ThreadPool>(100);
        std::vector<std::future<void>> futures;
        for (std::size_t i = 0; i < overlaps.size(); i++)
        {
            futures.emplace_back(threads->Submit([&](std::size_t i) -> void
                                                 {
                for(std::size_t j = 0; j < overlaps[i].size(); j++){
                    auto lhs = sequences[i]->InflateData(overlaps[i][j].lhs_begin, overlaps[i][j].lhs_end - overlaps[i][j].lhs_begin);
                    biosoup::NucleicAcid rhs_ ("", sequences[overlaps[i][j].rhs_id]->InflateData(overlaps[i][j].rhs_begin, overlaps[i][j].rhs_end - overlaps[i][j].rhs_begin));
                    if(!overlaps[i][j].strand) rhs_.ReverseAndComplement();
                    auto rhs = rhs_.InflateData();

                    /*
                    if(sequences[i]->name == "read=1,reverse,position=2675766-2676093,length=327,NC_000913.3_mutated") {
                        std::cout << lhs << std::endl;
                        std::cout << rhs << std::endl;
                    }
                     */
                    overlaps[i][j].alignment = ksw2_wrapper(i, overlaps[i][j], lhs, rhs);
                    /*
                    if(sequences[i]->name == "read=1,reverse,position=2675766-2676093,length=327,NC_000913.3_mutated"){
                        std::cout<<overlaps[i][j].alignment<<std::endl;
                        std::cout<<"-----------------------------------------"<<std::endl;
                    }
                     */
                } }, i));
        }

        for (auto &future : futures)
        {
            future.wait();
        }
        std::cout << "pairwise alignment done" << std::endl;
    }

public:
    std::vector<std::vector<biosoup::Overlap>> *get_overlaps() override
    {
        if (overlaps.empty())
        {
            this->overlaps = std::vector<std::vector<biosoup::Overlap>>(sequences.size());
            load_paf();
        }

        return &overlaps;
    }
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> *get_sequences() override
    {
        return &sequences;
    }

    explicit PafOverlapSource(std::string &args)
    {
        std::cout << "args: " << args << std::endl;
        auto spacePos = std::find(args.begin(), args.end(), ' ');
        std::cout << "reads: " << args.substr(0, std::distance(args.begin(), spacePos)) << std::endl;
        std::cout << "paf_file: " << args.substr(std::distance(args.begin(), spacePos) + 1) << std::endl;
        this->reads = args.substr(0, std::distance(args.begin(), spacePos));
        this->paf_file = args.substr(std::distance(args.begin(), spacePos) + 1);
        load_sequences();
    }
};

extern "C" OverlapSource *__attribute__((visibility("default"))) create(std::string &args)
{
    return (OverlapSource *)new PafOverlapSource(args);
}
