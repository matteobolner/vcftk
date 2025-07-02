import numpy as np
import pandas as pd
from vcftk.parsing import (
    build_var_ID,
    get_variants_metadata,
    get_variants_info,
    get_var_format_from_vcf,
    setup_samples_and_vcf,
    get_variants_stats,
)
from typing import Dict
from cyvcf2 import VCF, Writer


def set_pandas_display_options() -> None:
    """Set pandas display options."""
    # Ref: https://stackoverflow.com/a/52432757/
    display = pd.options.display
    display.max_columns = 1000
    display.max_rows = 10_000
    display.max_colwidth = 199
    display.width = 1000
    # display.precision = 2  # set as needed
    # display.float_format = lambda x: '{:,.2f}'.format(x)  # set as needed


set_pandas_display_options()


def setup(
    vcf_path,
    sample_info=None,
    sample_id_column="sample",
    threads=1,
    build_ids=True,
    # add_info=False,
):
    sample_info, vcf = setup_samples_and_vcf(
        input_vcf=vcf_path,
        input_sample_info=sample_info,
        sample_id_column=sample_id_column,
    )
    vcf.set_threads(threads)
    """
    variants = get_vcf_metadata(VCF(vcf_path), add_info)

    if create_ids_if_none:
        variants["ID"] = add_variant_ids(variants)
    variants = variants.set_index("ID")
    var_ids = list(variants.index)
    format_info = get_vcf_format_info(
        vcf=vcf, formats=variants["FORMAT"].str.split(":").explode().unique()
    )
    """
    return VCFClass(
        sample_id_column=sample_id_column,
        vcf=vcf,
        vcf_path=vcf_path,
        sample_info=sample_info,
        samples=vcf.samples,
        threads=threads,
        # variants=variants,
        # var_ids=var_ids,
        # format_info=format_info,
        # create_ids_if_none=create_ids_if_none,
        # add_info=add_info,
    )


class VCFClass:
    def __init__(
        self,
        vcf,
        vcf_path,
        sample_info,
        sample_id_column="sample",
        samples=[],
        threads=1,
        # ID BUILDING
        build_ids=False,
        verbose_ids=True,
        # variants=pd.DataFrame(),
        # var_ids=[],
        # format_info=pd.DataFrame(),
        # create_ids_if_none=False,
        # add_info=False,
    ):
        self.vcf = vcf
        self._vcf_path = vcf_path
        self._sample_id_column = sample_id_column
        self._threads = threads
        # self.variants = variants
        # self.var_ids = []
        # self.format_info = format_info
        self.sample_info = sample_info
        self.samples = samples
        self._build_ids = build_ids
        self._verbose_ids = verbose_ids
        # self._created_ids = create_ids_if_none
        # self._added_info = add_info
        # print(
        #    f"VCF contains {self.vcf.num_records} variants over {len(self.samples)} samples"
        # )

    def build_var_ids(self):
        ids = build_var_ID(self._variants_, alleles=self._verbose_ids)
        self._variants_["ID"] = ids
        self._variants_.set_index("ID", drop=True, inplace=True)
        return ids

    @property
    def variants(self):
        """Lazy initialization of variant metadata."""
        if not hasattr(self, "_variants_"):
            var_metadata = get_variants_metadata(self.vcf)
            if var_metadata["ID"].duplicated().any():
                Warning(
                    "There are duplicate / empty variant IDs - you must create unique IDs before proceeding, or problems will arise"
                )
            if self._build_ids:
                self.build_var_ids()
            else:
                self._variants_ = var_metadata.set_index("ID", drop=False)
        self.reset_vcf_iterator()
        return self._variants_

    @property
    def var_ids(self):
        if not hasattr(self, "_variants_"):
            _ = self.variants
        self._var_ids_ = list(self._variants_.index)
        self.reset_vcf_iterator()
        return self._var_ids_

    @property
    def var_info(self):
        """Lazy initialization of variant info fields."""
        if not hasattr(self, "_var_info_"):
            if not hasattr(self, "_var_ids_"):
                _ = self.var_ids
            var_info = get_variants_info(self.vcf)
            var_info["ID"] = self._var_ids_
            if var_info["ID"].duplicated().any():
                Warning(
                    "There are duplicate / empty variant IDs - you must create unique IDs before proceeding, or problems will arise"
                )
            self._var_info_ = var_info.set_index("ID", drop=False)
        self.reset_vcf_iterator()
        return self._var_info_

    def var_stats(self):
        if not hasattr(self, "_var_ids_"):
            _ = self.var_ids
        var_stats = get_variants_stats(self.vcf)
        var_stats["ID"] = self._var_ids_
        if var_stats["ID"].duplicated().any():
            Warning(
                "There are duplicate / empty variant IDs - you must create unique IDs before proceeding, or problems will arise"
            )
        self.reset_vcf_iterator()
        return var_stats

    def split(self, by="samples", columns=[]):
        if by == "samples":
            return self._split_by_samples(columns=columns)
        elif by == "variants":
            raise ValueError("Not yet implemented...")

    def _split_by_samples(self, columns) -> Dict[str, "VCF"]:
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
        for name, group in self.sample_info.groupby(by=columns):
            print(name)
            tempclass = setup(
                vcf_path=self._vcf_path,
                sample_info=group,
                sample_id_column=self._sample_id_column,
                threads=self._threads,
            )
            split_data[name] = tempclass
        return split_data

    def subset(self, what="samples", ids=[]):
        if what == "samples":
            return self._subset_samples(ids)
        elif what == "variants":
            raise ValueError("Not implemented yet...")

    def _subset_samples(self, ids):
        samples = self.sample_info.loc[ids]
        return setup(
            vcf_path=self._vcf_path,
            sample_info=samples,
            sample_id_column=self._sample_id_column,
            threads=self._threads,
        )

    # TODO: find a way to subset variants, cyvcf apparently cant do it
    """
    def _subset_variants(self, variant_list):
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

    def save_vcf(self, save_path, add_ids=False, var_ids=None):
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

    def display_genotypes(self):
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
        if not hasattr(self, "_var_ids_"):
            _ = self.var_ids
        allgts = [var.genotypes for var in self.vcf]
        gts = [[genotype_string(gt) for gt in i] for i in allgts]
        genotypes = pd.DataFrame(gts, index=self.var_ids, columns=self.samples)
        self.reset_vcf_iterator()
        return genotypes

    def get_genotypes(self, phasing_status=False):
        if not phasing_status:
            genotypes = pd.DataFrame(
                [[gt[:2] for gt in var.genotypes] for var in self.vcf]
            )
        else:
            genotypes = pd.DataFrame([[gt for gt in var.genotypes] for var in self.vcf])
        self.reset_vcf_iterator()
        genotypes.index = self.var_ids
        genotypes.columns = self.samples
        return genotypes

    def allele_frequency(self, var_ids="all", allele=1):
        genotypes = self.get_genotypes(phasing_status=False).transpose()
        if var_ids != "all":
            genotypes = genotypes.loc[var_ids]
        return compute_all_allele_frequencies(all_vars_gts=genotypes, allele=allele)

    def grouped_allele_frequency(self, var_ids="all", allele=1, groupcols=[]):
        genotypes = self.get_genotypes(phasing_status=False).transpose()
        if var_ids != "all":
            genotypes = genotypes[var_ids]
        genotypes = self.sample_info.merge(genotypes, left_index=True, right_index=True)
        afs = genotypes.groupby(by=groupcols).apply(
            compute_all_allele_frequencies, allele=allele
        )
        afs = afs.drop(columns=self.sample_info.columns).transpose()
        return afs


class Genotype(object):
    __slots__ = ("alleles", "phased")

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)

    __repr__ = __str__


def genotype_string(gt):
    sep = "/|"[int(gt[2])]
    gt_string = f"{gt[0]}{sep}{gt[1]}"
    return gt_string


def genotype_array(genotype_series):
    return np.array(genotype_series.values.tolist())


def compute_allele_frequency(var_gts, allele=1):
    ref_counts = np.sum(var_gts == 0)
    allele_counts = np.sum(var_gts == allele)
    if ref_counts == 0 and allele_counts == 0:
        return np.nan
    return allele_counts / (ref_counts + allele_counts)


def compute_all_allele_frequencies(all_vars_gts: pd.DataFrame, allele: int):
    afs = all_vars_gts.apply(
        lambda var_gts: compute_allele_frequency(
            genotype_array(var_gts), allele=allele
        ),
        axis=0,
    )
    return afs


def extract_vep_annotations(VCFClass, add_to_info=False, canonican_only=True):
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
    if "CSQ" not in VCFClass.variants.columns:
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
    vep_annotations = vep_annotations.reset_index().drop_duplicates().set_index("ID")
    if add_to_info:
        vep_annotations = self.variants.drop(columns=["CSQ"]).merge(
            vep_annotations, left_index=True, right_index=True
        )

    return vep_annotations
