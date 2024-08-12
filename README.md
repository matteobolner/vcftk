# vcforge
vcforge is a python library for working with VCF files.  

It provides a variety of tools and methods for importing, exploring, manipulating and analyzing Variant Call Format (VCF) data.

VCF parsing is handled by the cyvcf2 library (https://github.com/brentp/cyvcf2), which is based on the htslib library (https://github.com/samtools/htslib) written in C, resulting in very fast parsing.

vcforge utilizes pandas for producing tables and statistics.

The library is designed for explorative analyses, by linking the information about the samples with the variant data.


### Features

* Import VCF data along with sample information
* Automatically extract basic information about the variants and optionally assign variant IDs
* Extract INFO fields for all variants
* Extract selected FORMAT fields for all samples
* Split the variant data based on sample information (i.e. breed, population, etc)
* Get variant statistics (n called, call rate, allele frequency, nucleotide diversity, variant types and subtypes )


## How to install ##

Install with `pip`:

```shell
pip install vcforge
```

## Dependencies ##
* pandas
* numpy
* cyvcf2


## Why vcforge? ##

I started working on this library for myself, seeing how I found myself increasingly frustrated with the lack of python VCF parsing libraries that would allow me to easily explore VCF files in a readable and workable format, such as pandas DataFrames, and allow me to connect sample-level information with the variants. For example, splitting the variants by sample group is something that allows me to quickly and easily explore variants without having to rely on command-line tools that produce files that need to be examined.


## How to use ##

!!! Documentation is still a work in progress !!!

Import the library and read a VCF with its sample information
```shell
from vcforge import VCFClass

dataset=VCFClass(vcf_path='variants.vcf.gz', sample_info='samples.tsv',sample_id_column="sample", add_info=True, create_ids_if_none=True, threads=4)


```
Show the variant information

```shell

var_info=dataset.variants

```
Split the dataset by sample column(s) and get all genotypes of the variants in one of the populations as a pandas dataframe

```shell

split_dataset_by_breed=dataset.split_by_sample_column(column='population')

split_dataset_by_breed['POP_1'].show_genotypes()
```

Calculate some variant statistics about the selected population, that will be added to the variant information in the VCFClass

```shell

split_dataset_by_breed['POP_1'].get_var_stats(add_to_info=True)

```

Filter out variants with call rate < 0.9

```shell

var_info=split_dataset_by_breed['POP_1'].variants

selected_variants=var_info[var_info['CALL_RATE']>0.9]

```

Save the filtered subset as a VCF, with the ID column previously created.

```shell

split_dataset_by_breed['POP_1'].save_vcf(save_path='pop_1.vcf.gz', add_ids=True, var_ids=selected_variants.index, samples=None)
```
