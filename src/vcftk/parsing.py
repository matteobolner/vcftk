from numpy._core.fromnumeric import var
import pandas as pd
import gzip
import numpy as np

from cyvcf2 import VCF


def get_cyvcf(vcf_path):
    return VCF(vcf_path)


def get_var_info_from_var(var):
    return {k: v for k, v in var.INFO}


def setup_samples_and_vcf(input_vcf, input_sample_info=None, sample_id_column="sample"):
    """
    Setup the class with samples, and raw vcf dataframe.

    This function loads the data and sets up the class instance.
    Parameters:
        samples (DataFrame or str): DataFrame or path of file containing sample metadata.
        vcf (DataFrame or str): DataFrame or path of file with the vcf raw data.

    Raises:
        ValueError: If the sample ID column is not found in the data.
    """
    vcf = VCF(input_vcf)
    sample_info = pd.DataFrame(vcf.samples, columns=["test"])
    if input_sample_info is None:
        sample_info = pd.DataFrame(vcf.samples, columns=[sample_id_column])
    else:
        sample_info = parse_table(input_sample_info)
    if sample_info.index.name == sample_id_column:
        sample_info = sample_info.reset_index()
    if sample_info[sample_id_column].duplicated().any():
        raise ValueError(
            "Warning: there are duplicate values in the chosen sample column."
        )
    sample_info[sample_id_column] = sample_info[sample_id_column].astype(str)
    sample_info.set_index(sample_id_column, inplace=True)
    samples = [i for i in sample_info.index if i in vcf.samples]
    sample_info = sample_info.loc[samples]
    vcf.set_samples(samples)
    return sample_info, vcf


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


def get_variants_metadata(cyvcf):
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"]
    vars_metadata = [get_var_metadata_from_var(var) for var in cyvcf]
    vars_metadata_df = pd.DataFrame(vars_metadata, columns=columns)
    return vars_metadata_df


def get_variants_info(cyvcf):
    vars_info = [get_var_info_from_var(var) for var in cyvcf]
    vars_info = pd.DataFrame(vars_info)
    return vars_info


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
            + df["ALT"].apply(lambda x: ",".join(x) if isinstance(x, list) else x)
        )
    else:
        ids = df["CHROM"].astype(str) + ":" + df["POS"].astype(str)
    return ids


def build_var_ID_HGVS(var):
    # TODO: complete function
    chrom = var.CHROM
    pos = var.POS
    ref = var.REF
    alt = var.ALT
    if var.type == "snp":
        id = f"{chrom}:g.{pos}{ref}>{alt}"
    elif var.type == "mnp":
        id = f"{chrom}:g.{pos}{ref}>{','.join(alt)}"
    elif var.type == "indel":
        if var.subtype == "del":
            id = f"{chrom}:g.{pos}del"
        elif var.subtype == "ins":
            id = f"{chrom}:g.{pos}ins{alt}"

    return id
