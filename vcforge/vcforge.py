import pandas as pd
from vcforge.parsing import parse_input, build_var_ID
from typing import Dict


class VCFClass:
    def __init__(self, sample_id_column="sample", vcf=None, samples=None):
        self._sample_id_column = sample_id_column
        self.sample_info, self.vcf = self._setup_data(samples, vcf)
        self.samples = list(self.sample_info.index)
        self.var_info = self.vcf.drop(columns=self.samples)
        self.var_data = self.vcf[self.samples]

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
        try:
            # parse input data
            parsed_samples = parse_input(input_sample_info)
            parsed_vcf = parse_input(input_vcf)
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
            vcf = parsed_vcf
            if vcf.index.name == "ID":
                vcf.reset_index(inplace=True)
            vcf.columns = [str(i) for i in vcf.columns]
            if len(vcf["ID"].unique()) == 1:
                vcf["ID"] = build_var_ID(vcf)
            vcf.set_index("ID", inplace=True)
            # Drop samples not in data
            sample_info = sample_info[sample_info.index.isin(vcf.columns)]
        except ValueError as ve:
            raise ValueError(f"Error setting up data: {ve}")
        return sample_info, vcf

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
        split_data: Dict[str, VCF] = {}
        for name, group in self.sample_info.groupby(by=column):
            var_data = self.vcf[list(group.index)]
            vcf = pd.concat([self.var_info, var_data], axis=1)
            tempclass = VCFClass(
                sample_id_column=self._sample_id_column, vcf=vcf, samples=group
            )
            split_data[name] = tempclass
        return split_data


a = VCFClass(
    samples="/home/pelmo/work/workspace/PRSSV_candidate_genes/config/samples.tsv",
    vcf="/home/pelmo/work/workspace/PRSSV_candidate_genes/data/variants/gatk/vep_annotation/VIM/snvs_biallelic.vcf",
)




a.var_info.columns


t = a.split_by_sample_column("species")


a.var_info


var_info = df.drop(columns=samples["sample"].tolist())
var_info = var_info.set_index("ID")
var_info_first_cols = list(var_info.columns)
var_info_first_cols.remove("INFO")

var_info = var_info.apply(add_info_fields_to_row, axis=1)

var_info = var_info[
    var_info_first_cols + [i for i in var_info.columns if i not in var_info_first_cols]
]
var_info = var_info.drop(columns=["INFO"])

df = df.set_index("ID")

var_data = df[["FORMAT"] + samples["sample"].tolist()]


a.vcf
