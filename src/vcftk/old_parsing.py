import pandas as pd
import gzip
import numpy as np
from vcforge.utils import parse_table
from cyvcf2 import VCF


vcf_path="/home/pelmo/work/workspace/PRSSV_candidate_genes/data/variants/gatk/vep_annotation/VIM/snvs_multiallelic.vcf"
#vcf_path="/home/pelmo/work/workspace/PRSSV_candidate_genes/data/variants/gatk/filtered/snvs.vcf.gz"

def get_cyvcf(vcf_path):
    return VCF(vcf_path)

def get_var_info_from_var(var):
    return {k:v for k,v in var.INFO}

def get_var_metadata_from_var(var):
    var_info=[var.CHROM, var.POS,var.ID,var.REF,var.ALT,var.QUAL,var.FILTER,var.FORMAT]
    return var_info

def get_vcf_metadata(cyvcf):
    vars_metadata=[]
    for var in cyvcf:
        vars_metadata.append(get_var_metadata_from_var(var))
    vars_metadata=pd.DataFrame(vars_metadata, columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
    return vars_metadata

def get_vcf_info_fields(cyvcf):
    vars_info=[]
    for var in cyvcf:
        vars_info.append(get_var_info_from_var(var))
    return

def get_vcf_metadata_and_info(cyvcf):
    vars_info=[]
    vars_metadata=[]
    for var in cyvcf:
        vars_metadata.append(get_var_metadata_from_var(var))
        vars_info.append(get_var_info_from_var(var))
    vars_info=pd.DataFrame(vars_info)
    vars_metadata=pd.DataFrame(vars_metadata, columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
    return pd.concat([vars_metadata, vars_info], axis=1)


def get_var_format_from_var(cyvcf, format, allele):
    vars_format=[]
    for var in cyvcf:
        try:
            var_format=var.format("AD").transpose()[allele]
        except:
            var_format=np.full(len(cyvcf.samples), np.nan)
        vars_format.append(var_format)

    return pd.DataFrame(vars_format, columns=cyvcf.samples)












def read_vcf_to_dataframe(vcf_path):
    """
    Function to read a VCF file into a pandas DataFrame.

    Parameters:
    vcf_path (str): Path to the VCF file.

    Returns:
    pandas.DataFrame: DataFrame containing the VCF data.
    """
    # Read the VCF file into a DataFrame
    vcf=VCF(vcf_path)
    samples=vcf.samples
    allvars_metadata=[]
    allvars_info=[]
    allvars_data=[]
    for var in vcf:
        allvars_metadata.append(get_var_metadata_from_var(var))
        allvars_info=get_var_metadata_from_var(var)

    df=pd.DataFrame(allvars, columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])

    if gzipped:
        # Get the column names from the gzipped VCF file
        vcf_cols = get_vcf_colnames(vcf_path)

        # Read the gzipped VCF file into a DataFrame
        df = pd.read_csv(
            vcf_path, compression="gzip", header=None, sep="\t", comment="#"
        )
    else:
        # Get the column names from the VCF file
        vcf_cols = get_vcf_colnames(vcf_path, gzipped=False)

        # Read the VCF file into a DataFrame
        df = pd.read_csv(vcf_path, header=None, sep="\t", comment="#")

    # Strip the "#" from the first column name
    vcf_cols[0] = vcf_cols[0].lstrip("#")
    # Set the column names of the DataFrame
    df.columns = vcf_cols

    # Return the DataFrame
    return df

def get_vcf_colnames(vcf_path, gzipped=True):
    """
    Function to extract column names from a VCF file.

    Parameters:
    vcf_path (str): Path to the VCF file.
    gzipped (bool): Indicates whether the file is gzipped or not. Defaults to True.

    Returns:
    list: List of column names.
    """
    # Open the file in read mode
    if gzipped:
        with gzip.open(vcf_path, "rt") as ifile:
            # Iterate over the lines in the file
            for line in ifile:
                # Check if the line starts with #CHROM
                if line.startswith("#CHROM"):
                    # Extract column names from the line
                    vcf_names = [x.replace("\n", "") for x in line.split("\t")]
                    break
    else:
        with open(vcf_path, "r") as ifile:
            # Iterate over the lines in the file
            for line in ifile:
                # Check if the line starts with #CHROM
                if line.startswith("#CHROM"):
                    # Extract column names from the line
                    vcf_names = [x.replace("\n", "") for x in line.split("\t")]
                    break
    # Close the file
    ifile.close()
    # Return the list of column names
    return vcf_names


def read_vcf_to_dataframe(vcf_path):
    """
    Function to read a VCF file into a pandas DataFrame.

    Parameters:
    vcf_path (str): Path to the VCF file.

    Returns:
    pandas.DataFrame: DataFrame containing the VCF data.
    """
    # Check if the file is gzipped or not
    if vcf_path.endswith("gz"):
        gzipped = True
    else:
        gzipped = False

    # Read the VCF file into a DataFrame
    if gzipped:
        # Get the column names from the gzipped VCF file
        vcf_cols = get_vcf_colnames(vcf_path)

        # Read the gzipped VCF file into a DataFrame
        df = pd.read_csv(
            vcf_path, compression="gzip", header=None, sep="\t", comment="#"
        )
    else:
        # Get the column names from the VCF file
        vcf_cols = get_vcf_colnames(vcf_path, gzipped=False)

        # Read the VCF file into a DataFrame
        df = pd.read_csv(vcf_path, header=None, sep="\t", comment="#")

    # Strip the "#" from the first column name
    vcf_cols[0] = vcf_cols[0].lstrip("#")
    # Set the column names of the DataFrame
    df.columns = vcf_cols

    # Return the DataFrame
    return df


def parse_input(input_data):
    """
    Function to parse input data and return it as a pandas DataFrame.

    Parameters:
    input_data (pandas.DataFrame or str): Input data to be parsed.

    Returns:
    pandas.DataFrame: Input data as a pandas DataFrame.

    Raises:
    ValueError: If the input is not a pandas DataFrame or a file path to a VCF or gzipped VCF file.
    """
    # Check if the input is a pandas DataFrame
    if isinstance(input_data, pd.DataFrame):
        # Parse the input as a table
        return parse_table(input_data)
    # Check if the input is a file path
    elif isinstance(input_data, str):
        # Check if the file path ends with "vcf" or "vcf.gz"
        if input_data.endswith(("vcf", "vcf.gz")):
            # Read the VCF file into a pandas DataFrame and return it
            return read_vcf_to_dataframe(input_data)
        elif input_data.endswith(("tsv", "csv")):
            return parse_table(input_data)

    # If the input is neither a pandas DataFrame nor a file path, raise an error
    else:
        raise ValueError("Invalid format of input data or file path")


def get_info_fields(info_string):
    """
    Function to extract information fields from an INFO string.

    Parameters:
    info_string (str): The INFO string from a VCF entry.

    Returns:
    dict: A dictionary with the information field names as keys and their values as values.
    """
    # Initialize an empty dictionary to store the information fields
    info_dict = {}

    # Split the INFO string into individual fields
    split_info = info_string.split(";")

    # Iterate over each field
    for i in split_info:
        # Split the field into name and value
        info_name = i.split("=")[0]  # Name of the information field
        info_value = i.split("=")[1]  # Value of the information field

        # Add the name and value to the dictionary
        info_dict[info_name] = info_value

    # Return the dictionary of information fields
    return info_dict


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


def add_info_fields_to_row(row):
    info_dict = get_info_fields(row["INFO"])
    for k, v in info_dict.items():
        row[k] = v
    return row


def get_genotype_from_string(input_string):
    return input_string.split(":")[0]


def allele_frequency(genotypes_list):
    total_alleles = len(genotypes_list) * 2
    alleles = [i.split("/") for i in genotypes_list]
    alleles = sum(alleles, [])
    counts = np.unique(alleles, return_counts=True)
    freqs = {i: j / total_alleles for i, j in zip(counts[0], counts[1])}
    return freqs
