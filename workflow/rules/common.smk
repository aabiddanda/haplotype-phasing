import numpy as np


rule tabix_index:
    """Tabix-index any file."""
    input:
        "{pattern}",
    output:
        "{pattern}.tbi",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "tabix -f {input}"


rule index_bcf:
    """Generate the index for a BCF file."""
    input:
        bcf="{prefix}.bcf",
    output:
        csi="{prefix}.bcf.csi",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index -f {input.bcf}"


checkpoint list_chromosomes:
    """List the chromosomes available in any of the files."""
    input:
        vcfs=expand("results/per_chrom_inputs/{{outfix}}/{chrom}.vcf.gz", chrom=CHROM),
        tbis=expand(
            "results/per_chrom_inputs/{{outfix}}/{chrom}.vcf.gz.tbi", chrom=CHROM
        ),
    output:
        "checkpoints/{outfix}_chrom.list",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "for i in {input.vcfs}; do tabix -l $i; done > {output}"


checkpoint list_samples:
    """List the samples per vcf file."""
    input:
        vcf="{vcf_file}",
        tbi="{vcf_file}.tbi",
    output:
        "checkpoints/samples/{vcf_file}.samples.txt",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools query -l {input.vcf} > {output}"


def check_sample_order(wildcards):
    """Error checking to catch errors in sample naming."""
    manifest = pd.read_csv(analysis_configs[wildcards.outfix]["manifest"], sep="\t")
    samples = []
    for c in manifest.vcf_files:
        filename = checkpoints.list_samples.get(vcf_file=c).output[0]
        samples.append(np.loadtxt(filename, dtype=str))
    if len(samples) == 1:
        return manifest.vcf_files
    elif len(samples) > 1:
        for i in range(1, len(samples)):
            if samples[0].size != samples[i].size:
                raise ValueError(
                    "Each VCF file must contain the same number of samples!"
                )
            if ~np.all(samples[0] == samples[i]):
                raise ValueError("All VCF files must have the same sample order!")
        return manifest.vcf_files


rule split_per_chrom:
    """Split a VCF file into specific chromosomes."""
    input:
        vcf="{vcf_file}",
        tbi="{vcf_file}.tbi",
    output:
        temp("results/per_chrom/{vcf_file}/{chrom}.vcf.gz"),
    threads: 4
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view -r {wildcards.chrom} --threads {threads} {input} -Oz -o {output}"


rule combine_per_chrom:
    """Combine variants per chromosome."""
    input:
        vcfs=lambda wildcards: expand(
            "results/per_chrom/{pattern}/{{chrom}}.vcf.gz",
            pattern=check_sample_order(wildcards),
        ),
        tbis=lambda wildcards: expand(
            "results/per_chrom/{pattern}/{{chrom}}.vcf.gz.tbi",
            pattern=check_sample_order(wildcards),
        ),
    output:
        "results/per_chrom_inputs/{outfix}/{chrom}.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools concat -a -D -Ou {input.vcfs} | bcftools sort -Oz -o {output}"


rule convert_vcf2bcf:
    """Convert VCF formatted files to BCF for faster I/O when running phasing."""
    input:
        unphased_vcf="results/per_chrom_inputs/{outfix}/{chrom}.vcf.gz",
        unphased_vcf_tbi="results/per_chrom_inputs/{outfix}/{chrom}.vcf.gz.tbi",
    output:
        bcf="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf",
    threads: 4
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view -r {wildcards.chrom} {input.unphased_vcf} --threads {threads} -Ob -o {output.bcf}"
