#!python3


rule phase_beagle4:
    """Non-reference phasing using shapeit4."""
    input:
        unphased_vcf="results/per_chrom_inputs/{outfix}/{chrom}.vcf.gz",
        unphased_vcf_tbi="results/per_chrom_inputs/{outfix}/{chrom}.vcf.gz.tbi",
        reference_panel=lambda wildcards: extract_ref_panel(wildcards)[0],
        referance_panel_idx=lambda wildcards: extract_ref_panel(wildcards)[1],
    output:
        phased_vcf="results/beagle/{outfix}.{genome_build}.{chrom}.vcf.gz",
    log:
        "logs/beagle/{outfix}.{genome_build}.{chrom}.output.log",
    benchmark:
        "benchmark/beagle/{outfix}.{genome_build}.{chrom}.tsv"
    params:
        seed=lambda wildcards: analysis_configs[wildcards.outfix]["tools"]["beagle"][
            "seed"
        ],
        burnin=lambda wildcards: analysis_configs[wildcards.outfix]["tools"]["beagle"][
            "burnin"
        ],
        iterations=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "beagle"
        ]["iterations"],
        ref_panel=lambda wildcards, input: f"ref={input.reference_panel}"
        if input.reference_panel != []
        else "",
        outfix=lambda wildcards: f"results/beagle/{wildcards.outfix}.{wildcards.genome_build}.{wildcards.chrom}",
    threads: 8
    resources:
        time="4:00:00",
        mem="8G",
    conda:
        "../envs/beagle.yaml"
    shell:
        """
        beagle -Xmx8g gt={input.unphased_vcf} impute=false nthreads={threads} {params.ref_panel} burnin={params.burnin} iterations={params.iterations} seed={params.seed} out={params.outfix}  2>&1 | tee {log}
        """
