# 54gene workflow: 54gene phasing pipeline

This pipeline conducts internal population phasing for assorted datasets using pre-existing software.

## Authors

* Arjun Biddanda (@aabiddanda54gene)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Clone this repository to your local system, into the place where you want to perform the data analysis.
```
    git clone https://gitlab.com/data-analysis5/imputation/54gene-phasing.git
```

### Step 2: Configuration definitions & workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. The two key files are the `manifest` and `config.yaml` files. The manifest contains two columns:




#### Global Configuration

The global manifest file must have the following columns:

```
analysis_name	analysis_config
test1	config/analyses/test1.yaml
test2	config/analyses/test2.yaml
```

where each analysis is shown as a row.


#### Analysis-Specific Configuration

The manifest file for a specific analysis is:

```
vcf_file
-------------------
testdata/chr22.vcf.gz
```


The `vcf_file` indications can be relative paths to the top-level directory or absolute paths on your current system. We note that the files are automatically organized into "per-chromosome" inputs so this can be specified regardless if you have a merged VCF file or not.

##### Specifying Recombination Maps

One of the additional/optional specifications available is to define per-chromosome recombination maps:

```
recombination_maps: "config/recombination_maps/gwd.tsv"
```

, where the tab separated file specifies [HapMap-formatted](https://ftp.ncbi.nlm.nih.gov/hapmap/recombination/latest/) genetic maps:
```
chroms	recombination_map
chr22	resources/hg38/GWD/GWD_recombination_map_hapmap_format_hg38_chr_22.txt
```

NOTE: if this is not specified, then the default recombination map for each algorithm (largely the original HapMap Maps) will be used.

##### Specifying Reference Panels


We can also specify reference panels within the pipeline, using the current block within analysis-specific YAML file:

```
reference_panel: "config/ref_panels/hgdp.tsv"
```

where each file contains a specified reference file in per-chromosome format (in VCF/BCF format):

```
chroms	ref_panel
chr22	/home/ec2-user/local_data/hgdp_data/hgdp_wgs.20190516.statphase.autosomes.chr22.bcf
```

NOTE: that the header in this case must be provided and the file must be tab-separated.
NOTE: without this specified the default is to run phasing in "non-reference" mode


##### Phasing Algorithm Options

The various phasing algorithms are organized under the `tools` field in the YAML, with the `enabled` keyword defining whether the algorithm is run. The currently supported algorithms are:

* [SHAPEIT4](https://odelaneau.github.io/shapeit4/)
* [EAGLE2](https://alkesgroup.broadinstitute.org/Eagle/)

If you are interested in specific parameters for each algorithm please look them up on the manual pages. The relevant parameterizations we have included for the algorithms here:

###### Shapeit4

* `mcmc-iterations`: specifies the three different kinds of MCMC iteration (documented [here] (https://odelaneau.github.io/shapeit4/#documentation))
* `pbwt-depth`: number of conditioning haplotypes (default 4)
* `pbwt-mdr`: missing data rate (default 0.5)
* `pbwt-mac`: minimum minor allele count for phasing
* `seed`: Random number seed for MCMC initialization
* `sequencing`: sequencing mode for shapeit4


###### EAGLE2

* `kpbwt`: Number of PBWT conditioning haplotypes
* `pbwt-iters`: Number of PBWT iterations
* `expect-ibd`: Expected IBD length in centiMorgans
* `no-impute-missing`: Option to not impute missing data (default no)
* `geno-err-prob`: Error probability for HMM
* `hist-factor`: history factor for copying model


### Step 3: Install Environment

Install Snakemake and the baseline environment using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```
conda env create -f environment.yaml
```

You may also use [mamba](https://github.com/mamba-org/mamba) for faster dependency management.

### Step 4: Execute workflow

Activate the conda environment:
```
    conda activate 54gene-phasing
```
Test your configuration by performing a dry-run via
```
    snakemake --use-conda -n
```
Execute the workflow locally via
```
    snakemake --use-conda --cores $N
```
using `$N` cores or run it in a cluster environment via
```
    snakemake --use-conda --cluster qsub --jobs 100
```
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

All results will be under the `results/` directory and organized by algorithm (e.g. `results/shapeit4/`).

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/54gene-phasing.git` or `git remote add -f upstream https://github.com/snakemake-workflows/54gene-phasing.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.
