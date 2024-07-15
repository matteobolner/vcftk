
from cyvcf2 import VCF


a=VCF("/home/pelmo/work/workspace/PRSSV_candidate_genes/data/variants/gatk/vep_annotation/VIM/snvs_biallelic.vcf")


formats={}

for i in a:
    id=f"{i.CHROM}:{i.POS}_{i.POS}_{i.REF}_{",".join(i.ALT)}"

def build_var_ID(var):
    id=f"{var.CHROM}:{var.POS}_{var.REF}_{",".join(i.ALT)}"
    return id

def build_var_ID(df):
    ids = (
        df["CHROM"].astype(str)
        + ":"
        + df["POS"].astype(str)
        + "_"
        + df["REF"]
        + "_"
        + df["ALT"]
    )
    return ids


i.genotypes




help(a)

help(i)


for i in a:
    print(i.format("AD"))[0]



pd.DataFrame(i.format("AD")).transpose()

for i in a.ibd:
    print(i)


a.set_samples(['SS_69','SS_101'])






for var in a.variants():
    var.format("AD")
help(a)
