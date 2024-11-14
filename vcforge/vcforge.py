import pandas as pd
from vcforge.parsing import *
from vcforge.utils import *
from typing import Dict
from cyvcf2 import VCF, Writer


class VCFClass:
    def __init__(
        self,
        sample_id_column="sample",
        vcf_path=None,
        sample_info=None,
        add_info=False,
        create_ids_if_none=True,
        threads=1,
    ):
        self._vcf_path = vcf_path
        self._sample_id_column = sample_id_column
        self.sample_info, self.vcf = self._setup_data(sample_info, vcf_path)
        self.vcf.set_threads(threads)
        self.samples = self.vcf.samples
        self.variants = get_vcf_metadata(VCF(vcf_path), add_info=add_info)
        if create_ids_if_none:
            self.variants["ID"] = add_variant_ids(self.variants)
        self.variants = self.variants.set_index("ID")
        self.var_ids = list(self.variants.index)
        self.format_info = self._get_format_info()
        print(
            f"VCF contains {len(self.variants)} variants over {len(self.samples)} samples"
        )

    def _setup_data(self, input_sample_info, input_vcf) -> None:
        """
        Setup the class with samples, and raw vcf dataframe.

        This function loads the data and sets up the class instance.
        Parameters:
            samples (DataFrame or str): DataFrame or path of file containing sample metadata.
            vcf (DataFrame or str): DataFrame or path of file with the vcf raw data.

        Raises:
            ValueError: If the sample ID column is not found in the data.
        """
        parsed_samples = parse_table(input_sample_info)
        sample_info = parsed_samples
        if sample_info is None:
            raise ValueError("Sample metadata is not properly initialized.")
        if sample_info.index.name == self._sample_id_column:
            sample_info = sample_info.reset_index()
        if sample_info[self._sample_id_column].duplicated().any():
            raise ValueError(
                "Warning: there are duplicate values in the chosen sample column."
            )
        sample_info[self._sample_id_column] = sample_info[
            self._sample_id_column
        ].astype(str)
        sample_info.set_index(self._sample_id_column, inplace=True)
        vcf = VCF(input_vcf)
        samples = [i for i in sample_info.index if i in vcf.samples]
        vcf.set_samples(samples)
        sample_info = sample_info.loc[vcf.samples]
        return sample_info, vcf

    def _get_format_info(self):
        """
        Retrieve format information from the VCF file.

        This function extracts all unique format fields from the VCF file, retrieves
        their header information, and returns it as a DataFrame.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the header information for each unique format
            field in the VCF file.
        """
        all_formats = self.variants["FORMAT"].str.split(":").explode().unique()
        format_info = {}
        for i in all_formats:
            format_info[i] = self.vcf.get_header_type(i)
        format_info = pd.DataFrame(format_info).transpose()
        return format_info

    def split_by_sample_column(self, column: list) -> Dict[str, "VCF"]:
        """
        Split the dataset (data and sample metadata) in multiple independent VCF instances
        based on the values of one or more sample metadata columns.

        This function splits the dataset into multiple independent VCF instances, each
        containing a subset of the data based on the values of a sample metadata column. The
        function returns a dictionary containing the split data, where the dictionary keys are
        the unique values of the sample metadata column and the values are the VCF instances
        containing the split data.

        Args:
            column: The name of the column in the sample metadata DataFrame to use for splitting.

        Returns:
            A dictionary containing the split data, where the dictionary keys are the unique
            values of the sample metadata column and the values are the VCF instances
            containing the split data.
        """
        split_data: Dict[str, VCFClass] = {}
        for name, group in self.sample_info.groupby(by=column):
            print(name)
            tempclass = VCFClass(
                sample_id_column=self._sample_id_column,
                vcf_path=self._vcf_path,
                sample_info=group,
            )
            split_data[name] = tempclass
        return split_data

    # TODO: find a way to subset variants, cyvcf apparently cant do it
    """
    def subset_variants(self, variant_list):
        subsetted_variants=self.variants.loc[variant_list]
        subsetted = VCFClass(
            sample_id_column=self._sample_id_column,
            vcf_path=self._vcf_path,
            sample_info=self.sample_info,
        )
        subsetted.var_ids=variant_list
        subsetted.variants=subsetted.variants.loc[variant_ids]
        return subsetted
    """

    def reset_vcf_iterator(self):
        self.vcf = VCF(self._vcf_path)
        self.vcf.set_samples(self.samples)

    def format(self, format, allele):
        vars_format = get_var_format_from_vcf(self.vcf, format, allele)
        self.reset_vcf_iterator()
        return vars_format

    def subset_samples(self, samples):
        samples = self.sample_info.loc[samples]
        print(samples)
        return VCFClass(
            sample_id_column=self._sample_id_column,
            sample_info=samples,
            vcf_path=self._vcf_path,
        )

    def save_vcf(self, save_path, add_ids=False, var_ids=None, samples=None):
        w = Writer(save_path, self.vcf)
        vars_to_save = var_ids if var_ids is not None else self.var_ids
        for v, id in zip(self.vcf, self.var_ids):
            if id in vars_to_save:
                if add_ids is True:
                    v.ID = id
                w.write_record(v)
            else:
                pass
        w.close()
        self.reset_vcf_iterator()
        print(f"VCF saved to {save_path}")

    def get_var_stats(self, add_to_info=True):
        var_stats = {
            "NUM_CALLED": {},
            "CALL_RATE": {},
            "AA_FREQ": {},
            "NUCL_DIVERSITY": {},
            "VAR_TYPE": {},
            "VAR_SUBTYPE": {},
        }
        for var, id in zip(self.vcf, self.var_ids):
            var_stats["NUM_CALLED"][id] = var.num_called
            var_stats["CALL_RATE"][id] = var.call_rate
            var_stats["AA_FREQ"][id] = var.aaf
            var_stats["NUCL_DIVERSITY"][id] = var.nucl_diversity
            var_stats["VAR_TYPE"][id] = var.var_type
            var_stats["VAR_SUBTYPE"][id] = var.var_subtype
        var_stats = pd.DataFrame(var_stats)
        if add_to_info is True:
            self.variants = pd.concat([self.variants, var_stats], axis=1)
        self.reset_vcf_iterator()
        return var_stats

    def show_genotypes(self):
        """
        Return a DataFrame with the genotypes for each variant over the samples in the instance.

        The index of the DataFrame is the variant IDs, and the columns are the sample IDs.
        Each element of the DataFrame is a Genotype object, which can be used to access the genotype,
        phase, and read depths of the variant in the sample.

        Returns
        -------
        pandas.DataFrame
            DataFrame with the genotypes for each variant over the samples in the instance.
        """
        genotypes = []
        for var in self.vcf:
            genotypes.append([Genotype(i) for i in var.genotypes])
        genotypes = pd.DataFrame(genotypes, index=self.var_ids, columns=self.samples)
        self.reset_vcf_iterator()
        return genotypes

    def extract_vep_annotations(self, add_to_info=False):
        """
        Extract VEP annotations from the VCF file.
        This function explodes the "CSQ" column, which contains the VEP annotations, and
        creates a new DataFrame with one row per variant per transcript. The resulting
        DataFrame contains the VEP annotations for each variant, with the column names
        as described in the VCF header. If add_to_info is True, the variant info dataframe
        will be merged with the annotations. Keep in mind that there are likely multiple
        annotations per variant, therefore the resulting dataframe will have multiple rows
        per variant ID.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the VEP annotations for each variant. If add_to_info
            is True, the variant info dataframe will be merged with the annotations.

        Raises
        ------
        ValueError
            If the "CSQ" column is not found in the variants DataFrame.
        """
        if "CSQ" not in self.variants.columns:
            raise ValueError(
                "CSQ column not found in variants. This column is required for VEP annotations. Consider parsing VCF with add_info=True"
            )

        csq_info = (
            self.vcf.get_header_type("CSQ")["Description"]
            .split(" ")[6]
            .strip('"')
            .split("|")
        )
        csq_data = self.variants["CSQ"].str.split(",").explode()
        vep_annotations = csq_data.str.split("|", expand=True)
        vep_annotations.columns = csq_info
        vep_annotations = vep_annotations.replace("", np.nan)
        vep_annotations = (
            vep_annotations.reset_index().drop_duplicates().set_index("ID")
        )
        if add_to_info:
            vep_annotations = self.variants.drop(columns=["CSQ"]).merge(
                vep_annotations, left_index=True, right_index=True
            )

        return vep_annotations


class Genotype(object):
    __slots__ = ("alleles", "phased")

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)

    __repr__ = __str__
