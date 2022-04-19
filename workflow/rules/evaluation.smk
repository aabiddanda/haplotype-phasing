from lib.validators import FamFileValidator


rule validate_fam:
    """Validate the FAM file used as input."""
    input:
        sample_list="checkpoints/samples/results/{algo}/{outfix}.{genome_build}.{chrom}.vcf.gz.samples.txt",
    params:
        fam=lambda wildcards: analysis_configs[wildcards.outfix]["evaluation"]["fam"],
    output:
        temp(
            "checkpoints/fam_file/results/{algo}/{outfix}.{genome_build}.{chrom}.vcf.gz.fam"
        ),
    run:
        sample_ids = [x.rstrip() for x in open(input.sample_list).readlines()]
        fam_file_validator = FamFileValidator(params.fam)
        res = fam_file_validator.validate_fam(sample_ids)
        if res:
            fam_file_validator.fam_df.to_csv(
                output[0], sep=" ", index=False, header=False
            )  # noqa
        else:
            shell("touch {output}")


rule evaluate_switch_errors:
    """Rule to collect and evaluate phasing switch error per individual."""
    input:
        gen="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf",
        gen_index="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf.csi",
        hap="results/{algo}/{outfix}.{genome_build}.{chrom}.vcf.gz",
        hap_index="results/{algo}/{outfix}.{genome_build}.{chrom}.vcf.gz.tbi",
        fam_file="checkpoints/fam_file/results/{algo}/{outfix}.{genome_build}.{chrom}.vcf.gz.fam",
    output:
        iser="results/switch_error_estimation/{algo}/switch_error_rate.{outfix}.{genome_build}.{chrom}.iser",
        vser="results/switch_error_estimation/{algo}/switch_error_rate.{outfix}.{genome_build}.{chrom}.vser",
        mser=temp(
            "results/switch_error_estimation/{algo}/switch_error_rate.{outfix}.{genome_build}.{chrom}.mser.gz"
        ),
    log:
        "logs/switchError/{algo}.{outfix}.{genome_build}.{chrom}.output.log",
    params:
        outprefix="results/switch_error_estimation/{algo}/switch_error_rate.{outfix}.{genome_build}.{chrom}",
        maf=lambda wildcards: analysis_configs[wildcards.outfix]["evaluation"]["maf"],
    resources:
        time="1:00:00",
    conda:
        "../envs/switch_error.yaml"
    shell:
        "switchError --gen {input.gen} --hap {input.hap} --reg {wildcards.chrom} --fam {input.fam_file} --maf {params.maf} --out {params.outprefix} 2>&1 | tee {log}"
