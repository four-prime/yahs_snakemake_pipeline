import glob

configfile: "config.json"

sample_names = config["sample_names"]


def get_sample_input_paths():
    sample_fa_paths = {}
    sample_bam_paths = {}
    
    for sample in sample_names:
        for ext in ["fa", "fasta"]:
            fasta_path_list = glob.glob(f"samples/{sample}/input/*.{ext}")
            if fasta_path_list:
                sample_fa_paths[sample] = fasta_path_list
                break
        
        bam_path = glob.glob(f"samples/{sample}/input/*.bam")
        sample_bam_paths[sample] = bam_path

        if sample not in sample_fa_paths and sample not in sample_bam_paths:
            raise FileNotFoundError(f"No .fa, .fasta or .bam file found for sample {sample}")
        elif sample not in sample_fa_paths:
            raise FileNotFoundError(f"No .fa or .fasta file found for sample {sample}")
        elif sample not in sample_bam_paths:
            raise FileNotFoundError(f"No .bam file found for sample {sample}")

        if len(sample_fa_paths[sample]) > 1:
            raise ValueError(f"""
                            Multiple FASTA files found for sample {sample}: {sample_fa_paths[sample]}
                            Please either remove any accidentaly included files or concatenate your input fastas.
                             """)
    
    return sample_fa_paths, sample_bam_paths


SAMPLE_FA_PATHS, SAMPLE_BAM_PATHS = get_sample_input_paths()


def get_fasta_path(wildcards):
    return SAMPLE_FA_PATHS[wildcards.sample_name][0]


def get_fasta_index_path(wildcards):
    return f"{SAMPLE_FA_PATHS[wildcards.sample_name][0]}.fai"


def get_bam_path(wildcards):
    return SAMPLE_BAM_PATHS[wildcards.sample_name][0]

rule all:
    input:
        [f"samples/{sample_name}/output/out_JBAT.hic" for sample_name in sample_names]


rule samtools_faidx_source:
    container:
        "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"
    threads: 8
    resources:
        mem_mb=20000
    input:
        fa=get_fasta_path
    output:
        fai="samples/{sample_name}/input/filename.{ext}.fai"
    log:
        "samples/{sample_name}/logs/samtools_index_{ext}.log"
    wildcard_constraints:
        ext="(fa|fasta)"
    shell:
        "samtools faidx {input.fa} -o {output.fai} 2> {log}"


rule sambamba_deduplicate:
    container:
        "oras://community.wave.seqera.io/library/sambamba:1.0.1--69782c710ee49140"
    threads: 8
    resources:
        mem_mb=60000
    input:
        bam=get_bam_path,
    output:
        bam="samples/{sample_name}/work/sambamba_deduplicate/out.bam",
    log:
        "samples/{sample_name}/logs/sambamba_deduplicate.log"
    shell:
        "sambamba markdup -t {threads} --show-progress {input.bam} {output.bam}"


rule yahs:
    container:
        "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
    threads: 16
    resources:
        mem_mb=200000
    input:
        bam="samples/{sample_name}/work/sambamba_deduplicate/out.bam",
        fasta=get_fasta_path,
        fai=get_fasta_index_path,
    output:
        agp="samples/{sample_name}/work/yahs/out_scaffolds_final.agp",
        fa="samples/{sample_name}/work/yahs/out_scaffolds_final.fa",
        bin="samples/{sample_name}/work/yahs/out.bin",
    log:
        "samples/{sample_name}/logs/yahs.log"
    shell:
        "yahs -o {wildcards.sample_name}/work/yahs/out {input.fasta} {input.bam} 2> {log}"

rule output_scaffolds:
    input:
        fa="samples/{sample_name}/work/yahs/out_scaffolds_final.fa",
    output:
        fa="samples/{sample_name}/output/out_scaffolds_final.fa",
    log:
        "samples/{sample_name}/logs/output_scaffolds.log"
    shell:
        "mv {input.fa} {output.fa} 2> {log}"

rule juicer_generate_JBAT:
    container:
        "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
    threads: 8
    resources:
        mem_mb=32000
    input:
        agp="samples/{sample_name}/work/yahs/out_scaffolds_final.agp",
        bam="samples/{sample_name}/work/sambamba_deduplicate/out.bam",
        fai=get_fasta_index_path,
    output:
        txt="samples/{sample_name}/work/out_JBAT.txt",
        assembly="samples/{sample_name}/work/out_JBAT.assembly",
        jbat_log="samples/{sample_name}/work/out_JBAT.log"
    log:
        "samples/{sample_name}/logs/juicer_generate_JBAT.log"
    shell:
        "juicer pre -a -o {wildcards.sample_name}/work/out_JBAT {input.bam} {input.agp} {input.fai} > {log} 2> {output.jbat_log}"


rule download_juicer_tools:
    output:
        juicer_tools="downloads/juicer_tools.1.9.9_jcuda.0.8.jar",
    log:
        "downloads/logs/get_juicer_tools.log"
    shell:
        "wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar -O {output.juicer_tools} 2> {log}"


rule juicer_tools_generate_hic:
    container:
        "oras://community.wave.seqera.io/library/java-jdk:8.0.112--c971247410d7cba8"
    threads: 8
    resources:
        mem_mb=33000
    input:
        juicer_tools="utils/juicer_tools.1.9.9_jcuda.0.8.jar",
        txt="samples/{sample_name}/work/out_JBAT.txt",
        jbat_log="samples/{sample_name}/work/out_JBAT.log"
    output:
        hic="samples/{sample_name}/output/out_JBAT.hic"
    log:
        "samples/{sample_name}/logs/juicer_tools_generate_hic.log"
    shell:
        "(java -jar -Xmx{resources.mem_mb}m {input.juicer_tools} pre {input.txt} {output.hic} <(cat {input.jbat_log} | grep PRE_C_SIZE | awk '{{print $2\" \"$3}}')) 2> {log}"
