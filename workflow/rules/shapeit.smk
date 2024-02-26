#!python3
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import re

HTTP = HTTPRemoteProvider()


localrules:
    download_shapeit4_genetic_maps,


rule download_shapeit4_genetic_maps:
    """Download genetic maps for Shapeit4."""
    input:
        shapeit4_tar=HTTP.remote(
            "github.com/odelaneau/shapeit4/archive/refs/tags/v4.2.2.tar.gz"
        ),
    output:
        genetic_maps_maps=expand(
            "resources/shapeit4-4.2.2/maps/{chrom}.{x}.gmap.gz",
            chrom=CHROM,
            x="b38",
        ),
    params:
        build="b38",
    shell:
        """
        mv {input} resources/shapeit4.tar.gz
        cd resources/
        tar -zxvf shapeit4.tar.gz shapeit4-4.2.2/maps/
        cd shapeit4-4.2.2/maps/
        tar -xvf genetic_maps.{params.build}.tar.gz
        """


def get_shapeit4_recmap(wildcards):
    """Get the recombination-map set for shapeit4."""
    if analysis_configs[wildcards.outfix]["recombination_maps"] == "":
        return f"resources/shapeit4-4.2.2/maps/{wildcards.chrom}.{wildcards.genome_build}.gmap.gz"
    if analysis_configs[wildcards.outfix]["recombination_maps"] == "constant":
        return f"results/recomb_maps/shapeit4_constant/{wildcards.outfix}/{wildcards.chrom}.{wildcards.genome_build}.gmap.gz"
    else:
        return f"results/recomb_maps/shapeit4/{wildcards.outfix}/{wildcards.chrom}.{wildcards.genome_build}.gmap.gz"


rule create_constant_recombination_map:
    """Create a constant recombination map with 1cM equal to 1 Megabase."""
    input:
        unphased_bcf="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf",
        unphased_bcf_idx="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf.csi",
    output:
        gmap="results/recomb_maps/shapeit4_constant/{outfix}/{chrom}.{genome_build}.gmap.gz",
    shell:
        'bcftools query -f "%POS\t%CHROM" {input.unphased_bcf} | awk \'BEGIN{{c = 0; print "pos\tchr\tcM"}}; {{print $1"\t"$2"\t"c*1e-6; c+=$1;}}\' | gzip > {output.gmap}'


rule phase_shapeit4:
    """Non-reference phasing using shapeit4."""
    input:
        unphased_bcf="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf",
        unphased_bcf_idx="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf.csi",
        reference_panel=lambda wildcards: extract_ref_panel(wildcards)[0],
        referance_panel_idx=lambda wildcards: extract_ref_panel(wildcards)[1],
        genetic_map=lambda wildcards: get_shapeit4_recmap(wildcards),
    output:
        phased_vcf="results/shapeit4/{outfix}.{genome_build}.{chrom}.vcf.gz",
    log:
        "logs/shapeit4/{outfix}.{genome_build}.{chrom}.output.log",
    benchmark:
        "benchmark/shapeit4/{outfix}.{genome_build}.{chrom}.tsv"
    params:
        seed=lambda wildcards: analysis_configs[wildcards.outfix]["tools"]["shapeit4"][
            "seed"
        ],
        sequencing=lambda wildcards: "--sequencing"
        if analysis_configs[wildcards.outfix]["tools"]["shapeit4"]["sequencing"]
        else "",
        pbwt_depth=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "shapeit4"
        ]["pbwt-depth"],
        pbwt_mdr=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "shapeit4"
        ]["pbwt-mdr"],
        pbwt_mac=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "shapeit4"
        ]["pbwt-mac"],
        mcmc_iterations=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "shapeit4"
        ]["mcmc-iterations"],
        ref_panel=lambda wildcards, input: f"--reference {input.reference_panel}"
        if input.reference_panel != []
        else "",
    threads: 8
    resources:
        time="4:00:00",
        mem="8G",
    conda:
        "../envs/shapeit4.yaml"
    shell:
        """
        shapeit4 --input {input.unphased_bcf}\
            --map {input.genetic_map}\
            --seed {params.seed}\
            --region {wildcards.chrom}\
            --mcmc-iterations {params.mcmc_iterations}\
            {params.sequencing}\
            {params.ref_panel}\
            --thread {threads}\
            --output {output.phased_vcf} 2>&1 | tee {log}
        """


rule phase_shapeit5:
    """Non-reference phasing using shapeit4."""
    input:
        unphased_bcf="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf",
        unphased_bcf_idx="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf.csi",
        reference_panel=lambda wildcards: extract_ref_panel(wildcards)[0],
        referance_panel_idx=lambda wildcards: extract_ref_panel(wildcards)[1],
        genetic_map=lambda wildcards: get_shapeit4_recmap(wildcards),
    output:
        phased_vcf="results/shapeit5/{outfix}.{genome_build}.{chrom}.vcf.gz",
    log:
        "logs/shapeit5/{outfix}.{genome_build}.{chrom}.output.log",
    benchmark:
        "benchmark/shapeit5/{outfix}.{genome_build}.{chrom}.tsv"
    params:
        seed=lambda wildcards: analysis_configs[wildcards.outfix]["tools"]["shapeit5"][
            "seed"
        ],
        sequencing=lambda wildcards: "--sequencing"
        if analysis_configs[wildcards.outfix]["tools"]["shapeit5"]["sequencing"]
        else "",
        pbwt_depth=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "shapeit5"
        ]["pbwt-depth"],
        pbwt_mdr=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "shapeit5"
        ]["pbwt-mdr"],
        pbwt_mac=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "shapeit5"
        ]["pbwt-mac"],
        mcmc_iterations=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "shapeit5"
        ]["mcmc-iterations"],
        ref_panel=lambda wildcards, input: f"--reference {input.reference_panel}"
        if input.reference_panel != []
        else "",
    threads: 8
    resources:
        time="4:00:00",
        mem="8G",
    conda:
        "../envs/shapeit5.yaml"
    shell:
        """
        SHAPEIT5_phase_common --input {input.unphased_bcf}\
            --map {input.genetic_map}\
            --seed {params.seed}\
            --region {wildcards.chrom}\
            --mcmc-iterations {params.mcmc_iterations}\
            {params.sequencing}\
            {params.ref_panel}\
            --thread {threads}\
            --output {output.phased_vcf} 2>&1 | tee {log}
        """
