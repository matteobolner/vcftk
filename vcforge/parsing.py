import pandas as pd
import gzip
import numpy as np

from cyvcf2 import VCF


def get_cyvcf(vcf_path):
    return VCF(vcf_path)


def get_var_info_from_var(var):
    return {k: v for k, v in var.INFO}


def get_var_metadata_from_var(var):
    var_metadata = [
        var.CHROM,
        var.POS,
        var.ID,
        var.REF,
        ",".join(var.ALT),
        var.QUAL,
        var.FILTER,
        ":".join(var.FORMAT),
    ]
    return var_metadata


def get_vcf_metadata(cyvcf, add_info=False):

    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"]
    if add_info:
        vars_metadata = []
        vars_info = []
        for var in cyvcf:
            vars_metadata.append(get_var_metadata_from_var(var))
            vars_info.append(get_var_info_from_var(var))
        vars_info = pd.DataFrame(vars_info)
        vars_metadata_df = pd.DataFrame(vars_metadata, columns=columns)
        vars_metadata_df = pd.concat([vars_metadata_df, vars_info], axis=1)

    else:
        vars_metadata = [get_var_metadata_from_var(var) for var in cyvcf]
        vars_metadata_df = pd.DataFrame(vars_metadata, columns=columns)

    return vars_metadata_df


def add_variant_ids(vars_metadata):
    if vars_metadata["ID"].isna().all():
        if vars_metadata[["CHROM", "POS"]].duplicated().any():
            return build_var_ID(vars_metadata, alleles=True)
        else:
            return build_var_ID(vars_metadata, alleles=False)
    else:
        return vars_metadata["ID"]


def get_vcf_info(cyvcf, info_name):
    vars_info = [var.INFO[info_name] for var in cyvcf]
    vars_info = pd.Series(vars_info, name=info_name)
    return vars_info


def get_all_vcf_info(cyvcf):
    vars_info = [get_var_info_from_var(var) for var in cyvcf]
    vars_info = pd.DataFrame(vars_info)
    return vars_info


def get_var_format_from_vcf(cyvcf, format, allele):
    vars_format = []
    # ids = []
    for var in cyvcf:
        # ids.append(f"{var.CHROM}:{var.POS}")
        try:
            var_format = var.format(format).transpose()[allele]
        except:
            var_format = np.full(len(cyvcf.samples), np.nan)
        vars_format.append(var_format)
    var_format_df = pd.DataFrame(vars_format, columns=cyvcf.samples)
    var_format_df = var_format_df.replace(-2147483648, np.nan)
    return var_format_df

def build_var_ID(df, alleles=False):
    if alleles == True:
        ids = (
            df["CHROM"].astype(str)
            + ":"
            + df["POS"].astype(str)
            + "_"
            + df["REF"]
            + "_"
            + df["ALT"].apply(
                lambda x: ",".join(x) if isinstance(x, list) else x
            )
        )
    else:
        ids = df["CHROM"].astype(str) + ":" + df["POS"].astype(str)
    return ids
