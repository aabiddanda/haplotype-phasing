from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


localrules:
    download_shapeit4_genetic_maps,


rule download_shapeit4_genetic_maps:
    """Download  genetic maps for Shapeit4."""
    input:
        shapeit4_tar=HTTP.remote("github.com/odelaneau/shapeit4/archive/refs/tags/v4.2.2.tar.gz"),
    output:
        genetic_maps_maps=expand(
            "resources/shapeit4-4.2.2/maps/{chrom}.{x}.gmap.gz",
            chrom=["chr%d" % (i + 1) for i in range(22)] + ["chrX"],
            x=config["genome_build"],
        ),
    params:
        build=config["genome_build"],
    shell:
        """
        mv {input} resources/shapeit4.tar.gz
        cd resources/
        tar -zxvf shapeit4.tar.gz shapeit4-4.2.2/maps/
        cd shapeit4-4.2.2/maps/
        tar -xvf genetic_maps.{params.build}.tar.gz
        """


rule phase_shapeit4:
    """Non-reference phasing using shapeit4."""
    input:
        unphased_bcf="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf",
        unphased_bcf_idx="results/per_chrom_inputs/{outfix}.{genome_build}.{chrom}.bcf.csi",
        genetic_map="resources/shapeit4-4.2.2/maps/{chrom}.{genome_build}.gmap.gz",
    output:
        phased_vcf="results/shapeit4/{outfix}.{genome_build}.{chrom}.vcf.gz",
    log:
        "logs/shapeit4/{outfix}.{genome_build}.{chrom}.output.log",
    params:
        seed=config["tools"]["shapeit4"]["seed"],
        sequencing="--sequencing" if config["tools"]["shapeit4"]["sequencing"] else "",
        pbwt_depth=config["tools"]["shapeit4"]["pbwt-depth"],
        pbwt_mdr=config["tools"]["shapeit4"]["pbwt-mdr"],
        pbwt_mac=config["tools"]["shapeit4"]["pbwt-mac"],
        mcmc_iterations=config["tools"]["shapeit4"]["mcmc-iterations"],
    threads: 16
    resources:
        time="4:00:00",
        mem="8G",
    conda:
        "../envs/shapeit4.yaml"
    shell:
        """
        shapeit4 --input {input.unphased_bcf} --map {input.genetic_map} --region {wildcards.chrom} --mcmc-iterations {params.mcmc_iterations} {params.sequencing} --thread {threads} --output {output.phased_vcf} 2>&1 | tee {log}
        """
