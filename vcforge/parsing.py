import pandas as pd
import gzip
import numpy as np
from vcforge.utils import parse_table
from cyvcf2 import VCF


def get_cyvcf(vcf_path):
    return VCF(vcf_path)


def get_var_info_from_var(var):
    return {k: v for k, v in var.INFO}


def get_var_metadata_from_var(var):
    var_info = [
        var.CHROM,
        var.POS,
        var.ID,
        var.REF,
        ",".join(var.ALT),
        var.QUAL,
        var.FILTER,
        ":".join(var.FORMAT),
    ]
    return var_info


def get_vcf_metadata(cyvcf):
    vars_metadata = []
    for var in cyvcf:
        vars_metadata.append(get_var_metadata_from_var(var))
    vars_metadata = pd.DataFrame(
        vars_metadata,
        columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"],
    )
    return vars_metadata


def get_vcf_info_fields(cyvcf):
    vars_info = []
    for var in cyvcf:
        vars_info.append(get_var_info_from_var(var))
    return


def get_vcf_metadata_and_info(cyvcf):
    vars_info = []
    vars_metadata = []
    for var in cyvcf:
        vars_metadata.append(get_var_metadata_from_var(var))
        vars_info.append(get_var_info_from_var(var))
    vars_info = pd.DataFrame(vars_info)
    vars_metadata = pd.DataFrame(
        vars_metadata,
        columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"],
    )
    if len(vars_metadata["ID"].unique()) == 1:
        if vars_metadata[["CHROM", "POS"]].duplicated().any():
            vars_metadata["ID"] = build_var_ID(vars_metadata, alleles=True)
        else:
            vars_metadata["ID"] = build_var_ID(vars_metadata, alleles=False)
    df = pd.concat([vars_metadata, vars_info], axis=1)
    df = df.set_index("ID")
    return df


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
            + df["ALT"].str.join(",")
        )
    else:
        ids = df["CHROM"].astype(str) + ":" + df["POS"].astype(str)
    return ids
