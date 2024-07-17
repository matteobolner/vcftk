import pandas as pd
from vcforge.parsing import *
from typing import Dict
from cyvcf2 import VCF, Writer


class VCFClass:
    def __init__(
        self,
        sample_id_column="sample",
        vcf_path=None,
        sample_info=None,
        add_info=False,
        threads=1,
    ):
        self._vcf_path = vcf_path
        self._sample_id_column = sample_id_column
        self.sample_info, self.vcf = self._setup_data(sample_info, vcf_path)
        self.vcf.set_threads(threads)
        self.samples = self.vcf.samples
        if add_info is True:
            self.variants = get_vcf_metadata_and_info(VCF(vcf_path))
        else:
            self.variants = get_vcf_metadata(VCF(vcf_path))
        self.var_ids = list(self.variants.index)
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

    def split_by_sample_column(
        self, column: list, add_info: bool = False
    ) -> Dict[str, "VCF"]:
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
                add_info=add_info,
            )
            split_data[name] = tempclass
        return split_data

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
            add_info=False,
        )

    def save_vcf(self, save_path, add_ids=False):
        w = Writer(save_path, self.vcf)
        for v, id in zip(self.vcf, self.var_ids):
            if add_ids is True:
                v.ID = id
            w.write_record(v)
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
        genotypes=[]
        for var in self.vcf:
            genotypes.append([Genotype(i) for i in var.genotypes])
        genotypes=pd.DataFrame(genotypes, index=self.var_ids, columns=self.samples)
        self.reset_vcf_iterator()
        return genotypes


class Genotype(object):
    __slots__ = ("alleles", "phased")

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)

    __repr__ = __str__
