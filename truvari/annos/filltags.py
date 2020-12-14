import sys
import numpy as np
import pysam

from collections import Counter

v = pysam.VariantFile(sys.argv[1])

def calc_hwe(nref, nalt, nhet):
    """
    Returns tuple of the HWE and the ExeHet
    """
    #void calc_hwe(args_t *args, int nref, int nalt, int nhet, float *p_hwe, float *p_exc_het)
    ngt = (nref + nalt) // 2 # also, just nsamples
    nrare = min(nref, nalt)
    # We got an hts_expand array for het probabilities
    probs = np.zeros(nrare + 1, dtype=np.dtype('d'))
    
    # start at midpoint
    mid = nrare * (nref + nalt - nrare) // (nref + nalt)
    
    # check to ensure that midpoint and rare alleles have same parity
    if ( (nrare & 1) ^ (mid & 1) ):
        mid += 1

    het = mid
    hom_r = (nrare - mid) // 2
    hom_c = ngt - het - hom_r
    probs[mid] = 1
    my_sum = 1

    while het > 1:
        probs[het - 2] = probs[het] * het * (het - 1.0) / (4.0 * (hom_r + 1.0) * (hom_c + 1.0))
        my_sum += probs[het - 2]

        # 2 fewer heterozygotes for next iteration -> add rare, one common homozygote
        hom_r += 1
        hom_c += 1
        # loop
        het -= 2
    
    het = mid
    hom_r = (nrare - mid) // 2
    hom_c = ngt - het - hom_r
    while het <= nrare - 2:
        probs[het + 2] = probs[het] * 4.0 * hom_r * hom_c / ((het + 2.0) * (het + 1.0));
        my_sum += probs[het + 2]
        
        # add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
        hom_r -= 1
        hom_c -= 1
        # loop
        het += 2

    probs = probs / my_sum
    
    p_exc_het = probs[nhet:].sum()

    p_hwe = probs[probs > probs[nhet]].sum()
    if p_hwe > 1:
        p_hwe = 1;
    #p_hwe = 1 - prob;
    return p_exc_het, 1 - p_hwe

def allele_annos(entry):
    data = np.array([int(_) for _ in list(entry.samples.values()[0]["FT"])], dtype=int)
    n_samps = len(data)
    alt_cnt = sum(data)
    ref_cnt = n_samps * 2 - alt_cnt
    n_het = len(data[data == 1])
    cnt = Counter()
    cnt[0] = len(data[data == 0]) * 2 + len(data[data == 1])
    cnt[1] = data.sum()
    af = cnt[1] / (n_samps * 2)
    srt = [v for k, v in sorted(cnt.items(), key=lambda item: item[1])]
    maf = 1 - (srt[-1] / (n_samps * 2))
    
    p_exc_het, p_hwe = calc_hwe(cnt[0], cnt[1], n_het)
    return af, maf, p_exc_het, p_hwe


def cnt_allele_annos(entry):
    n_samps = 0
    ref_cnt = 0
    alt_cnt = 0
    n_het = 0
    cnt = Counter() # 0 or 1 allele counts
    for i in entry.samples.values():
        n_samps += 1
        for j in i["GT"]:
            cnt[j] += 1
        n_het += 1 if i["GT"][0] != i["GT"][1] else 0

    af = cnt[1] / (n_samps * 2)
    srt = [v for k, v in sorted(cnt.items(), key=lambda item: item[1])]
    maf = 1 - (srt[-1] / (n_samps * 2))
    
    p_exc_het, p_hwe = calc_hwe(cnt[0], cnt[1], n_het)
    return af, maf, p_exc_het, p_hwe


tolerance = 1e-6
for entry in v:
    print(allele_annos(entry))
    continue
    #af, maf, p_exe_het, p_hwe = allele_annos(entry)
    error = False
    if abs(p_hwe - entry.info["HWE"][0]) > tolerance:
        error = True
    if abs(af - entry.info["AF"][0]) > tolerance:
        error = True
    if abs(maf - entry.info["MAF"]) > tolerance:
        error =True
    if abs(p_exe_het - entry.info["ExcHet"][0]) > tolerance:
        error = True
    if error:
        print('af', af, entry.info["AF"])
        print('maf', af, entry.info["MAF"])
        print('eh', p_exe_het, entry.info["ExcHet"])
        print('hwe', p_hwe, entry.info["HWE"])


