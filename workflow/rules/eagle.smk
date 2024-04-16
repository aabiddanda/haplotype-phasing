#!python3


rule copy_eagle_genetic_maps:
    """Copy the EAGLE genetic maps."""
    output:
        temp(
            "resources/eagle/genetic_map_{}_withX.txt.gz".format(
                pc.convert_build_to_hg("b38")
            )
        ),
    params:
        hg_notation=pc.convert_build_to_hg("b38"),
    conda:
        "../envs/eagle.yaml"
    shell:
        """
        cp ${{CONDA_PREFIX}}/share/eagle/tables/genetic_map_{params.hg_notation}_withX.txt.gz {output}
        """


def get_eagle_recmaps(wildcards):
    if analysis_configs[wildcards.outfix]["recombination_maps"] == "":
        return "resources/eagle/genetic_map_{}_withX.txt.gz".format(
            pc.convert_build_to_hg(wildcards.genome_build)
        )
    else:
        return f"results/recomb_maps/eagle/{wildcards.outfix}/{wildcards.chrom}.gmap.gz"


rule phase_eagle:
    """Non-reference phasing using EAGLE2."""
    input:
        bcf_file="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf",
        bcf_file_idx="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf.csi",
        reference_panel=lambda wildcards: extract_ref_panel(wildcards)[0],
        referance_panel_idx=lambda wildcards: extract_ref_panel(wildcards)[1],
        genetic_map=get_eagle_recmaps,
    output:
        phased_vcf="results/eagle/{outfix}.{genome_build}.{chrom}.vcf.gz",
    log:
        "logs/eagle/{outfix}.{genome_build}.{chrom}.output.log",
    wildcard_constraints:
        chrom="chr[[0-9]+|X]",
    benchmark:
        "benchmark/eagle2/{outfix}.{genome_build}.{chrom}.tsv"
    params:
        outprefix="results/eagle/{outfix}.{genome_build}.{chrom}",
        kpbwt=lambda wildcards: analysis_configs[wildcards.outfix]["tools"]["eagle"][
            "kpbwt"
        ],
        pbwt_iters=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "eagle"
        ]["pbwt-iters"],
        expect_ibd=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "eagle"
        ]["expect-ibd"],
        imp_missing=lambda wildcards: "--noImpMissing"
        if analysis_configs[wildcards.outfix]["tools"]["eagle"]["no-impute-missing"]
        else "",
        hist_factor=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "eagle"
        ]["hist-factor"],
        geno_err_prob=lambda wildcards: analysis_configs[wildcards.outfix]["tools"][
            "eagle"
        ]["geno-err-prob"],
        ref_panel=lambda wildcards, input: f"--vcfRef={input.reference_panel}"
        if input.reference_panel != []
        else "",
        vcf_target=lambda wildcards, input: f"--vcfTarget={input.bcf_file}"
        if input.reference_panel != []
        else "--vcf={input.bcf_file}",
    threads: 8
    resources:
        time="4:00:00",
        mem="8G",
    conda:
        "../envs/eagle.yaml"
    shell:
        """
        eagle --chrom={wildcards.chrom}\
            --Kpbwt={params.kpbwt}\
            --pbwtIters={params.pbwt_iters}\
            --histFactor={params.hist_factor}\
            --genoErrProb={params.geno_err_prob}\
            --expectIBDcM={params.expect_ibd}\
            --numThreads={threads}\
            {params.imp_missing}\
            {params.vcf_target}\
            {params.ref_panel}\
            --geneticMapFile={input.genetic_map}\
            --outPrefix={params.outprefix} 2>&1 | tee {log}
        """
