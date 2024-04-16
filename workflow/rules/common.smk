#!python3
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


checkpoint list_chromosomes_full:
    input:
        vcf="{vcf_file}",
        tbi="{vcf_file}.tbi",
    output:
        "checkpoints/list_chromosomes/{vcf_file}.chrom.full.list",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "tabix -l {input.vcf} > {output}"


def uniq_chroms(wildcards):
    """Error checking to catch errors in sample naming."""
    manifest = pd.read_csv(analysis_configs[wildcards.outfix]["manifest"], sep="\t")
    chromosomes = []
    for c in manifest.vcf_files:
        filename = checkpoints.list_chromosomes_full.get(vcf_file=c).output[0]
        chromosomes.append(np.loadtxt(filename, dtype=str))
    uniq_chroms = np.unique(chromosomes)
    filt_chroms = [c for c in uniq_chroms if ("chrY" not in c) and ("chrX" not in c)]
    return filt_chroms


checkpoint list_chromosomes:
    """List the chromosomes available in any of the files."""
    input:
        vcfs=lambda wildcards: expand(
            "results/per_chrom_inputs/{{outfix}}/{chrom}.vcf.gz",
            chrom=uniq_chroms(wildcards),
        ),
        tbis=lambda wildcards: expand(
            "results/per_chrom_inputs/{{outfix}}/{chrom}.vcf.gz.tbi",
            chrom=uniq_chroms(wildcards),
        ),
    output:
        "checkpoints/list_chromosomes/{outfix}_chrom.list",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "for i in {input.vcfs}; do tabix -l $i; done > {output}"


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
        "results/per_chrom/{vcf_file}/{chrom}.vcf.gz",
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


checkpoint valid_reference_panel:
    """Checkpoint for validating the reference panel."""
    output:
        "checkpoints/ref_panel/{outfix}.ref_panel",
    run:
        if analysis_configs[wildcards.outfix]["reference_panel"] == "":
            shell("touch {output}")
        else:
            ref_panel_manifest = pd.read_csv(
                analysis_configs[wildcards.outfix]["reference_panel"], sep="\t"
            )
            for c in ref_panel_manifest.chroms:
                assert c in CHROM
            assert (
                ref_panel_manifest.chroms.size
                == np.unique(ref_panel_manifest.chroms.values).size
            )
            index_files = [
                pc.determine_index(r) for r in ref_panel_manifest.ref_panel.values
            ]
            ref_panel_manifest["file_index"] = index_files
            ref_panel_manifest.to_csv(str(output), sep="\t", index=False)


def extract_ref_panel(wildcards):
    """Extract the information about the reference panel."""
    checkpoint_file = checkpoints.valid_reference_panel.get(
        outfix=wildcards.outfix
    ).output[0]
    if Path(checkpoint_file).stat().st_size == 0:
        return [], []
    else:
        ref_panel_manifest = pd.read_csv(checkpoint_file, sep="\t").set_index("chroms")
        ref_panel = ref_panel_manifest.filter(
            [wildcards.chrom], axis=0
        ).ref_panel.values[0]
        ref_panel_idx = ref_panel_manifest.filter(
            [wildcards.chrom], axis=0
        ).file_index.values[0]
        return ref_panel, ref_panel_idx


rule convert_hapmap_to_formats:
    """Convert hapmap formatted genetic maps to shapeit4/eagle format."""
    output:
        temp("results/recomb_maps/{algo}/{outfix}/{chrom}.gmap.gz"),
    wildcard_constraints:
        algo="shapeit4|eagle",
    resources:
        mem="1G",
        time="0:30:00",
    run:
        if analysis_configs[wildcards.outfix]["recombination_maps"] == "":
            raise ValueError("Cannot convert a recombination map if none are provided!")
        else:
            recomb_map_manifest = pd.read_csv(
                analysis_configs[wildcards.outfix]["recombination_maps"],
                sep="\t",
                dtype=str,
            )
            for c in recomb_map_manifest.chroms:
                assert c in CHROM
            assert (
                recomb_map_manifest.chroms.size
                == np.unique(recomb_map_manifest.chroms.values).size
            )
            assert wildcards.chrom in recomb_map_manifest.chroms.values
            filename = recomb_map_manifest[
                recomb_map_manifest.chroms == wildcards.chrom
            ].recombination_map.values[0]
            transformed_df = pc.convert_hapmap_genmap(filename, wildcards.algo)
            transformed_df.to_csv(str(output), index=False, sep="\t")
