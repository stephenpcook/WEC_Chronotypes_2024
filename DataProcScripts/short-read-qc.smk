def get_short_reads(wildcards):
    row = samples.loc[wildcards.sample]
    rtnValue = [f'{config["short_read_dir"]}/{row.short_fwd}',
            f'{config["short_read_dir"]}/{row.short_rev}']
    return rtnValue

rule setup_human:
    output: directory("scratch/ref")
    threads: 16
    log:
        "/lustre/projects/Research_Project-T108347/WoTG/logs"
    benchmark:
        "benchmarks/setup_human.tsv"
    params:
        scratch=config['scratch_dir']
    shell:
        """
            cd {params.scratch}
            wget -O in.gz https://zenodo.org/record/1208052/files/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
            bbmap.sh -Xmx23g ref=in.gz 2>&1 | tee {log};
        """

rule run_bbmap_clumpify:
    input:
        get_short_reads
    output:
        "scratch/{sample}.clumped.fq.gz"
    threads: 16
    benchmark:
        "benchmarks/{sample}/run_bbmap_clumpify.tsv"
    log:
        "logs/{sample}/run_bbmap_clumpify.log"
    shell:
        """
            clumpify.sh -eoom -da in1={input[0]} in2={input[1]} out={output} dedupe optical 2>&1 | tee {log}
        """
rule run_bbmap_filter_by_tile:
    input:
        rules.run_bbmap_clumpify.output
    output:
        "scratch/{sample}.filtered_by_tile.fq.gz"
    threads: 16
    log:
        "logs/{sample}/run_bbmap_filter_by_tile.log"
    benchmark:
        "benchmarks/{sample}/run_bbmap_filter_by_tile.tsv"
    shell:
        """
            filterbytile.sh -eoom -Xmx120g -da in={input} out={output} 2>&1 | tee {log}
        """
rule run_bbmap_bbduk_remove_adapters:
    input:
        rules.run_bbmap_filter_by_tile.output
    output:
        "scratch/{sample}.trimmed.fq.gz"
    threads: 16
    log:
        "logs/{sample}/run_bbmap_bbduk_remove_adapters.log"
    benchmark:
        "benchmarks/{sample}/run_bbmap_bbduk_remove_adapters.tsv"
    shell:
        """
            bbduk.sh -Xmx120g -eoom -da in={input} out={output} ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=50 ref=adapters ftm=5 ordered 2>&1 | tee {log}
        """

rule run_bbmap_bbduk_remove_artefacts:
    input:
        rules.run_bbmap_bbduk_remove_adapters.output
    output:
        "scratch/{sample}.filtered.fq.gz"
    threads: 16
    log:
        "logs/{sample}/run_bbmap_bbduk_remove_artefacts.log"
    benchmark:
        "benchmarks/{sample}/run_bbmap_bbduk_remove_artefacts.tsv"
    shell:
        """
            bbduk.sh -Xmx120g -eoom -da in={input} out={output} k=31 ref=artifacts,phix ordered cardinality 2>&1 | tee {log}
        """

rule tadpole:
    input:
        rules.run_bbmap_bbduk_remove_artefacts.output
    output:
        "scratch/{sample}.filtered.ec.fq.gz"
    threads: 16
    log:
        "logs/{sample}/tadpole.log"
    benchmark:
        "benchmarks/{sample}/tadpole.tsv"
    shell:
        """
            tadpole.sh -Xmx120g -eoom -da in={input} out={output} mode=correct k=62 ecc ecco merge=t prealloc=t 2>&1 | tee {log}
        """
rule remove_human:
    input:
        rules.tadpole.output,
        rules.setup_human.output
    output:
        "scratch/{sample}.filtered.ec.no.hcd.fq.gz"
    threads: 16
    log:
        "logs/{sample}/remove_human.log"
    benchmark:
        "benchmarks/{sample}/remove_human.tsv"
    params:
        scratch=config['scratch_dir']
    shell:
        """
            bbmap.sh -Xmx120g -eoom -da minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
            path={params.scratch}/ qtrim=rl trimq=10 untrim in={input[0]} outu={output} 2>&1 | tee {log}
        """

rule split_reads:
    input:
        rules.remove_human.output
    output:
        fwd="output/{sample}/fwd.ec.no-human.hq.fq.gz",
        rev="output/{sample}/rev.ec.no-human.hq.fq.gz"
    threads: 16
    log:
        "logs/{sample}/split_reads.log"
    benchmark:
        "benchmarks/{sample}/split_reads.tsv"
    params:
        sample_seed=42,
        subsample_count=0
    shell:
        """
            reformat.sh -Xmx120g -eoom -da in={input} out1={output[0]} out2={output[1]} samplereadstarget={params.subsample_count} 2>&1 | tee {log}
        """
