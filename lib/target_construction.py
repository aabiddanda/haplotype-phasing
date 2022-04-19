"""Module for target construction for phasing."""


def construct_phasing_targets(config, chroms, outfix):
    """Construct the targets for the phasing rules."""
    targets = []
    for algo in config[outfix]["tools"]:
        if config[outfix]["tools"][algo]["enabled"]:
            for chrom in chroms:
                res = "results/{algo}/{outfix}.{genome_build}.{chrom}.vcf.gz".format(
                    outfix=outfix,
                    genome_build=config[outfix]["genome_build"],
                    chrom=chrom,
                    algo=algo,
                )
                targets.append(res)
                targets.append(f"{res}.tbi")
    return targets
