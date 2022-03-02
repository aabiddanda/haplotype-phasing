rule copy_eagle_genetic_maps:
    output:
        temp(
            "resources/eagle/genetic_map_{}_withX.txt.gz".format(
                pc.convert_build_to_hg(config["genome_build"])
            )
        ),
    params:
        hg_notation=pc.convert_build_to_hg(config["genome_build"]),
    conda:
        "../envs/eagle.yaml"
    shell:
        """
        cp ${{CONDA_PREFIX}}/share/eagle/tables/genetic_map_{params.hg_notation}_withX.txt.gz {output}
        """


rule phase_eagle:
    """Non-reference phasing using EAGLE2."""
    input:
        bcf_file="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf",
        bcf_file_idx="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf.csi",
        genetic_map="resources/eagle/genetic_map_{}_withX.txt.gz".format(
            pc.convert_build_to_hg(config["genome_build"])
        ),
    output:
        phased_vcf="results/eagle/{outfix}.{genome_build}.{chrom}.vcf.gz",
    log:
        "logs/eagle/{outfix}.{genome_build}.{chrom}.output.log",
    params:
        outprefix="results/eagle/{outfix}.{genome_build}.{chrom}",
        kpbwt=config["tools"]["eagle"]["kpbwt"],
        pbwt_iters=config["tools"]["eagle"]["pbwt-iters"],
        expect_ibd=config["tools"]["eagle"]["expect-ibd"],
        imp_missing="--noImpMissing" if config["tools"]["eagle"]["no-impute-missing"] else "",
        hist_factor=config["tools"]["eagle"]["hist-factor"],
        geno_err_prob=config["tools"]["eagle"]["geno-err-prob"],
    threads: 16
    resources:
        time="4:00:00",
        mem="8G",
    conda:
        "../envs/eagle.yaml"
    shell:
        """
        eagle --chrom={wildcards.chrom} --Kpbwt={params.kpbwt} --pbwtIters={params.pbwt_iters} --histFactor={params.hist_factor} --genoErrProb={params.geno_err_prob} --expectIBDcM={params.expect_ibd} --numThreads={threads} {params.imp_missing}  --vcf={input.bcf_file}   --geneticMapFile={input.genetic_map} --outPrefix={params.outprefix} 2>&1 | tee {log}
        """
