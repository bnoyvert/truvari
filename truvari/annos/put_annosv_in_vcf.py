"""
Given an AnnoSV *annotated.tsv and its input *vcf
Add all the new 'full' annotations to the vcf's INFO fields/headers
I'll output to stdout so I can pipe into vcf-sort and bgzip
"""
import sys
import glob
import gzip
import argparse

def parse_args(args):
    """
    parse args
    """
    parser = argparse.ArgumentParser(prog="put_annosv_in_vcf", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("anno", help="AnnoSV annotated.tsv (can be a glob)")
    parser.add_argument("vcf", help="VCF input to AnnoSV")
    parser.add_argument("-o", default=None, help="File to write (stdout)")
    args = parser.parse_args()
    return args

def make_anno_lookup(anno):
    """
    Create a dictionary of the variant and all the key/values of annotations put onto it
    """
    # The hard part is going to be getting the types correct so I can write the header well.. probably
    # Or I can just hardcode the header information, and then flat unstructure parse.
    # This is assuming that the INFO column is good and that the AnnoSV output doesnt change...
    # I think I might
    def make_key(dat):
        """
        make a key of chrom:start.ref.alt
        """
        return "%s:%s.%s.%s" % (data[1], data[2], data[7], data[8])
    
    def extract_anno(data, header):
        """
        Make key/values of each annotation... join with ';'. ensure they're 
        VCF compatiable formatting
        last 68 keys?
        """
        ret = []
        for v,k in zip(data[-len(header):], header):
            if v == "": continue
            if k == "AnnotSV_type" and v == "split": continue
            v = v.replace(";", ",").replace('/', ',').replace(' ','_')
            ret.append(k + "=" + v)
        return ";".join(ret)

    ret = {}
    for i in glob.glob(anno):
        with open(i) as fh:
            header = fh.readline().strip().split('\t')
            annos = [_.replace(' ', '_') for _ in header[header.index('AnnotSV type'):]]
            for line in fh:
                data = line.strip().split('\t')
                if data[-68] == "split": continue
                manno = extract_anno(data, annos)
                ret[make_key(data)] = extract_anno(data, annos)
    return ret

# no idea how to make the new header fields...
# Got 68 of these to translate

def main(args):
    """
    main runner
    """
    global n_header_info
    args = parse_args(args)
    anno_lookup = make_anno_lookup(args.anno)
    sys.stderr.write("%d annotated entries\n" % (len(anno_lookup)))
    v = gzip.GzipFile(args.vcf)
    cnt = 0
    nannod = 0
    for entry in v:
        cnt += 1
        entry = entry.decode()
        if entry.startswith("##"):
            sys.stdout.write(entry)
        elif entry.startswith("#CH"):
            sys.stdout.write(n_header_info)
            sys.stdout.write(entry)
        else:
            data = entry.strip().split('\t')
            # Ugh - AnnotSV has problems with grch38 nomenclature
            if data[0].startswith("chr"):
                chrom_key = data[0][3:]
            else:
                chrom_key = data[0]
            k = "%s:%s.%s.%s" % (chrom_key, data[1], data[3], data[4])
            if k in anno_lookup:
                data[7] += ';' + anno_lookup[k]
                nannod += 1
            sys.stdout.write("\t".join(data) + '\n')
    sys.stderr.write("%d entries written. %d entries annotated\n" % (cnt, nannod))

   
# Placing this here for safekeeping. 
n_header_info = """##INFO=<ID=Gene_name,Number=.,Type=String,Description="Gene symbol">
##INFO=<ID=NM,Number=.,Type=String,Description="Transcript symbol">
##INFO=<ID=CDS_length,Number=.,Type=Integer,Description="Length of the CoDing Sequence (CDS) (bp) overlapping with the SV">
##INFO=<ID=tx_length,Number=.,Type=Integer,Description="Length of transcript (bp) overlapping with the SV">
##INFO=<ID=location,Number=1,Type=String,Description="SV location in the gene (e.g.  txStart-exon1 ,  intron3-exon7 )">
##INFO=<ID=location2,Number=1,Type=String,Description="SV location in the gene's coding regions (e.g.  3'UTR-CDS )">
##INFO=<ID=intersectStart,Number=1,Type=Integer,Description="Start position of the intersection between the SV and the transcript">
##INFO=<ID=intersectEnd,Number=1,Type=Integer,Description="End position of the intersection between the SV and the transcript">
##INFO=<ID=DGV_GAIN_IDs,Number=.,Type=String,Description="DGV Gold Standard GAIN IDs overlapped with the annotated SV">
##INFO=<ID=DGV_GAIN_n_samples_with_SV,Number=1,Type=Integer,Description="Number of individuals with a shared DGV_GAIN_ID">
##INFO=<ID=DGV_GAIN_n_samples_tested,Number=1,Type=Integer,Description="Number of individuals tested">
##INFO=<ID=DGV_GAIN_Frequency,Number=1,Type=Float,Description="Relative GAIN Frequency: (DGV_GAIN_n_samples_with_SV / DGV_GAIN_n_samples_tested)">
##INFO=<ID=DGV_LOSS_IDs,Number=.,Type=String,Description="DGV Gold Standard LOSS overlapped with the annotated SV">
##INFO=<ID=DGV_LOSS_n_samples_with_SV,Number=1,Type=Integer,Description="Number of individuals with a shared DGV_LOSS_ID">
##INFO=<ID=DGV_LOSS_n_samples_tested,Number=1,Type=Integer,Description="Number of individuals tested">
##INFO=<ID=DGV_LOSS_Frequency,Number=1,Type=Float,Description="Relative LOSS Frequency: (DGV_LOSS_n_samples_with_SV / DGV_LOSS_n_samples_tested)">
##INFO=<ID=GD_ID,Number=.,Type=String,Description="gnomAD IDs overlapping the annotated SV with the same SV type">
##INFO=<ID=GD_AN,Number=.,Type=Integer,Description="gnomad total number of alleles genotyped (for biallelic sites) or individuals with copy-state estimates (for multiallelic sites)">
##INFO=<ID=GD_N_HET,Number=.,Type=Integer,Description="gnomAD number of individuals with heterozygous genotypes">
##INFO=<ID=GD_N_HOMALT,Number=.,Type=Integer,Description="gnomAD number of individuals with homozygous alternate genotypes">
##INFO=<ID=GD_AF,Number=1,Type=Float,Description="Maximum of the gnomAD allele frequency (for biallelic sites) and copy-state frequency (for multiallelic sites)">
##INFO=<ID=GD_POPMAX_AF,Number=1,Type=Float,Description="Maximum of the gnomAD maximum allele frequency across any population">
##INFO=<ID=GD_ID_others,Number=.,Type=String,Description="Other gnomAD IDs overlapping the annotated SV (with a different SV type)">
##INFO=<ID=DDD_SV,Number=.,Type=String,Description="Deciphering Developmental Disorders (DDD) SV coordinates from the DDD study (data control sets) overlapped with the annotated SV">
##INFO=<ID=DDD_DUP_n_samples_with_SV,Number=1,Type=Integer,Description="Number of individuals with a shared DDD_DUP">
##INFO=<ID=DDD_DUP_Frequency,Number=1,Type=Float,Description="DUP Frequency">
##INFO=<ID=DDD_DEL_n_samples_with_SV,Number=1,Type=Integer,Description="Number of individuals with a shared DDD_DEL">
##INFO=<ID=DDD_DEL_Frequency,Number=1,Type=Float,Description="DEL Frequency">
##INFO=<ID=1000g_event,Number=.,Type=String,Description="1000 genomes events (e.g. DEL, DUP, ALU...) overlapped with the annotated SV">
##INFO=<ID=1000g_AF,Number=1,Type=Float,Description="1000 genomes allele frequency">
##INFO=<ID=1000g_max_AF,Number=1,Type=Float,Description="Maximum observed allele frequency across the 1000 genomes populations">
##INFO=<ID=IMH_ID,Number=.,Type=String,Description="Ira M. Hall's lab IDs overlapping the annotated SV">
##INFO=<ID=IMH_AF,Number=1,Type=Float,Description="IMH Allele Frequency">
##INFO=<ID=IMH_ID_others,Number=.,Type=String,Description="Other IMH IDs overlapping the annotated SV (with a different SV type)">
##INFO=<ID=promoters,Number=.,Type=String,Description="List of the genes whose promoters are overlapped by the SV">
##INFO=<ID=dbVar_event,Number=.,Type=String,Description="dbVar NR SV event types (e.g. deletion, duplication)">
##INFO=<ID=dbVar_variant,Number=.,Type=String,Description="dbVar NR SV accession (e.g. nssv1415016)">
##INFO=<ID=dbVar_status,Number=1,Type=String,Description="dbVar NR SV clinical assertion (e.g. pathogenic, likely pathogenic)">
##INFO=<ID=TADcoordinates,Number=.,Type=String,Description="Coordinates of the TAD whose boundaries overlapped with the annotated SV (boundaries included in the coordinates)">
##INFO=<ID=ENCODEexperiments,Number=1,Type=String,Description="ENCODE experiments from where the TAD have been defined">
##INFO=<ID=GCcontent_left,Number=1,Type=String,Description="GC content around the left SV breakpoint (+/- 100bp)">
##INFO=<ID=GCcontent_right,Number=1,Type=String,Description="GC content around the right SV breakpoint (+/- 100bp)">
##INFO=<ID=Repeats_coord_left,Number=.,Type=String,Description="Repeats coordinates around the left SV breakpoint (+/- 100bp)">
##INFO=<ID=Repeats_type_left,Number=.,Type=String,Description="Repeats type around the left SV breakpoint (+/- 100bp). e.g. AluSp, L2b, L1PA2, LTR12C, SVA_D"> 
##INFO=<ID=Repeats_coord_right,Number=.,Type=String,Description="Repeats coordinates around the right SV breakpoint (+/- 100bp)">
##INFO=<ID=Repeats_type_right,Number=.,Type=String,Description="Repeats type around the right SV breakpoint (+/- 100bp). e.g. AluSp, L2b, L1PA2, LTR12C, SVA_D,">
##INFO=<ID=ACMG,Number=1,Type=String,Description="ACMG gene">
##INFO=<ID=synZ_ExAC,Number=1,Type=Float,Description="Positive synZ_ExAC (Z score) indicate gene intolerance to synonymous variation">
##INFO=<ID=misZ_ExAC,Number=1,Type=Float,Description="Positive misZ_ExAC (Z score) indicate gene intolerance to missense variation">
##INFO=<ID=pLI_ExAC,Number=1,Type=Float,Description="Score computed in the ExAc database indicating the probability that a gene is intolerant to a loss of function variation (Nonsense, splice acceptor and donor variants caused by SNV). ExAC consider pLI >= 0.9 as an extremely LoF intolerant set of genes">
##INFO=<ID=delZ_ExAC,Number=1,Type=Float,Description="Positive delZ_ExAC (Z score) indicate greater deletion intolerance">
##INFO=<ID=dupZ_ExAC,Number=1,Type=Float,Description="Positive dupZ_ExAC (Z score) indicate greater duplication intolerance">
##INFO=<ID=cnvZ_ExAC,Number=1,Type=Float,Description="Positive cnvZ_ExAC (Z score) indicate greater CNV intolerance">
##INFO=<ID=Phenotypes,Number=.,Type=String,Description="OMIM e.g. Charcot-Marie-Tooth disease">
##INFO=<ID=Inheritance,Number=.,Type=String,Description="e.g. AD (= 'Autosomal dominant')2">
##INFO=<ID=morbidGenesCandidates,Number=1,Type=String,Description="OMIM">
##INFO=<ID=Mim_Number,Number=.,Type=Integer,Description="OMIM unique six-digit identifier">
##INFO=<ID=morbidGenes,Number=1,Type=String,Description="OMIM">
##INFO=<ID=HI_DDDpercent,Number=1,Type=Float,Description="Haploinsufficiency ranks">
##INFO=<ID=DDD_status,Number=1,Type=String,Description="Deciphering Developmental Disorders (DDD) category. e.g. confirmed, probable, possible">
##INFO=<ID=DDD_mode,Number=1,Type=String,Description="Deciphering Developmental Disorders (DDD) allelic requirement. e.g. biallelic, hemizygous">
##INFO=<ID=DDD_consequence,Number=1,Type=String,Description="Deciphering Developmental Disorders (DDD) mutation consequence. e.g. loss of function, uncertain">
##INFO=<ID=DDD_disease,Number=.,Type=String,Description="Deciphering Developmental Disorders (DDD) disease name. e.g. 'OCULOAURICULAR SYNDROME'">
##INFO=<ID=DDD_pmids,Number=.,Type=String,Description="Deciphering Developmental Disorders (DDD) pmids">
##INFO=<ID=HI_CGscore,Number=1,Type=Integer,Description="HaploInsufficiency Score">
##INFO=<ID=TriS_CGscore,Number=1,Type=Integer,Description="TriploSensitivity Score">
##INFO=<ID=AnnotSV_ranking,Number=1,Type=Integer,Description="AnnotSV ACMG guideline classification">
"""

if __name__ == '__main__':
    main(sys.argv[1:])
 
