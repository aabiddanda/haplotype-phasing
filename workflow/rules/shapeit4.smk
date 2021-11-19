#!python3

rule download_shapeit4_genetic_maps:
    """ Download appropriate genetic maps for Shapeit4. """
    input:
       shapeit4_tar = 'https://github.com/odelaneau/shapeit4/archive/refs/tags/v4.2.2.tar.gz'
    output:
        b37_maps = expand('resources/shapeit4-4.2.2/maps/chr{chrom}.b37.gmap.gz', chrom=[i+1 for i in range(22)] + ['X'])
        b38_maps = expand('resources/shapeit4-4.2.2/maps/chr{chrom}.b38.gmap.gz', chrom=[i+1 for i in range(22)] + ['X'])
    shell:
    """
    wget {input.shapeit4_tar} -O resources/shapeit4.tar.gz
    tar -zxvf resources/shapeit4.tar.gz shapeit4-4.2.2/maps/
    cd resources/shapeit4-4.2.2/maps/
    tar -xvf genetic_maps.b37.tar.gz
    tar -xvf genetic_maps.b38.tar.gz
    """

# localrules: download_shapeit4_genetic_maps

rule phase_shapeit4:
  input:
    unphased_vcf = "{prefix}.vcf.gz"
    genetic_map = "resources/shapeit4-4.2.2/maps/chr{chrom,\d+}.{genome_build}.gmap.gz",
  output:
    phased_vcf = 'results/'
  log: 'logs/'
  params:
    output = "",
    pbwt_depth=config["tools"]["shapeit4"]["pbwt-depth"],
    pbwt_mdr= config["tools"]["shapeit4"]["pbwt-mdr"],
    pbwt_mac= config["tools"]["shapeit4"]["pbwt-mac"],
    mcmc_iterations=config["tools"]["shapeit4"]["mcmc-iterations"]
  threads: 8
  resources:
    time="4:00:00",
    mem="8G",
    nodes=1
  conda:
    "../envs/shapeit4.yaml"
  shell:
    """
    shapeit4 --input {input.unphased_vcf} --map {input.genetic_map} --region {wildcards.chrom} --pbwt-mdr {params.pbwt_mdr} --pbwt-mac {params.pbwt_mac} --mcmc-iterations {params.mcmc_iterations} --threads {threads} --output {output.vcf}
    """
