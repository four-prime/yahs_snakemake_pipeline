import glob

configfile: "config.json"

sample_names = config["sample_names"]


def get_sample_input_paths():
    sample_hifi_paths = {}
    sample_hic_r1_paths = {}
    sample_hic_r2_paths = {}
    
    for sample in sample_names:

        valid_extensions = ('.fa', '.fasta', '.fq', '.fastq', '.fa.gz', '.fasta.gz', '.fq.gz', '.fastq.gz')

        hifi_path_list = glob.glob(f"samples/{sample}/hifi/*")
        hifi_path_list = [f for f in hifi_path_list if f.endswith(valid_extensions)]

        hic_r1_list = glob.glob(f"samples/{sample}/hic/r1/*")
        hic_r1_list  = [f for f in hic_r1_list  if f.endswith(valid_extensions)]

        hic_r2_list = glob.glob(f"samples/{sample}/hic/r2/*")
        hic_r2_list = [f for f in hic_r2_list if f.endswith(valid_extensions)]

        hifi_count = len(hifi_path_list)
        r1_count = len(hic_r1_list)
        r2_count = len(hic_r2_list)

        # check for surplus Hi-C inputs
        if r1_count > 1:
            raise ValueError(
                f"Multiple R1 Hi-C inputs found for sample {sample}: {hic_r1_list}\n"
                f"Please provide only a single input."
            )
        if r2_count > 1:
            raise ValueError(
                f"Multiple R2 Hi-C inputs found for sample {sample}: {hic_r2_list}\n"
                f"Please provide only a single input."
            )

        # check for missing inputs
        missing = []
        if hifi_count == 0:
            missing.append("HiFi")
        if r1_count == 0:
            missing.append("Hi-C R1")
        if r2_count == 0:
            missing.append("Hi-C R2")

        if missing:
            raise FileNotFoundError(
                f"Missing required inputs for sample {sample}: {', '.join(missing)}\n"
                f"Please provide all required inputs."
            )

        sample_hifi_paths[sample] = hifi_path_list
        sample_hic_r1_paths[sample] = hic_r1_list[0]
        sample_hic_r2_paths[sample] = hic_r2_list[0]
   
    return sample_hifi_paths, sample_hic_r1_paths, sample_hic_r2_paths


SAMPLE_FA_PATHS, SAMPLE_HIC_R1_PATHS, SAMPLE_HIC_R2_PATHS = get_sample_input_paths()


def get_hifi_path(wildcards):
    return SAMPLE_FA_PATHS[wildcards.sample_name]

def get_hic_r1_path(wildcards):
    return SAMPLE_HIC_R1_PATHS[wildcards.sample_name]

def get_hic_r2_path(wildcards):
    return SAMPLE_HIC_R2_PATHS[wildcards.sample_name]


rule all:
    input:
        [f"output/{sample_name}/{sample_name}_JBAT.hic" for sample_name in sample_names],
        [f"output/{sample_name}/{sample_name}_JBAT.txt" for sample_name in sample_names],
        [f"output/{sample_name}/{sample_name}_JBAT.assembly" for sample_name in sample_names],
        [f"output/{sample_name}/{sample_name}_scaffolds_final.fa" for sample_name in sample_names]
        

rule generate_samplesheet:
    input:
        hic_1=get_hic_r1_path,
        hic_2=get_hic_r2_path,
    output:
        samplesheet="work/{sample_name}/hic_samplesheet/samplesheet.csv",
    shell:
        "echo -e 'sample,fastq_1,fastq_2\\n{wildcards.sample_name},{input.hic_1},{input.hic_2}' > {output.samplesheet}"


rule hifiasm:
  container:
    "oras://community.wave.seqera.io/library/hifiasm:0.25.0--bcfc60a944a26aaa"
  threads: 24
  resources:
    mem_mb=1_000_000
  input:
    hifi=get_hifi_path,
    hic_1=get_hic_r1_path,
    hic_2=get_hic_r2_path,
  output:
    unitigs="work/{sample_name}/hifiasm/{sample_name}.hic.p_utg.gfa",
    contigs="work/{sample_name}/hifiasm/{sample_name}.hic.p_ctg.gfa",
  log:
    "logs/{sample_name}/hifiasm.log"
  shell:
    # "module load snakemake && "
    "hifiasm -t {threads} -o work/{wildcards.sample_name}/hifiasm/{wildcards.sample_name} -l0 --h1 {input.hic_1} --h2 {input.hic_2} {input.hifi} 2> {log}"


rule gfa2fa:
  threads: 8
  resources:
    mem_mb=50_000
  container:
    "oras://community.wave.seqera.io/library/gfatools_tabix:b5cb1ec01cf54783"
  input:
    fa="work/{sample_name}/hifiasm/{sample_name}.hic.p_ctg.gfa"
  output:
    fa="work/{sample_name}/gfa2fa/{sample_name}.hic.p_ctg.fa"
  log:
    "logs/{sample_name}/gfa2fa.log"
  shell:
    "gfatools gfa2fa {input.fa} > {output.fa} 2> {log}"


rule hic_pro:
    handover: True
    input:
        samplesheet="work/{sample_name}/hic_samplesheet/samplesheet.csv",
        fa="work/{sample_name}/gfa2fa/{sample_name}.hic.p_ctg.fa",
        config="hic_pro.config"
    output:
        bam="work/{sample_name}/hicpro/hicpro/mapping/{sample_name}_0_bwt2pairs.bam",
    shell:
        "nextflow run nf-core/hic -r 2.1.0 "
        "-work-dir work/{wildcards.sample_name}/hicpro_work "
        "--outdir work/{wildcards.sample_name}/hicpro "
        "-c {input.config} "
        "-profile singularity "
        "--input {input.samplesheet} "
        "--fasta {input.fa} "
        "--dnase "
        "--min_cis_dist 1000 "
        "--bin_size 10000,25000,100000 "
        "--skip_maps "
        "--skip_dist_decay "
        "--skip_tads "
        "--skip_compartments "
        "--skip_balancing "
        "--skip_mcool "
        "--skip_multiqc "
        "--save_aligned_intermediates"


rule samtools_faidx_source:
    container:
        "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"
    threads: 8
    resources:
        mem_mb=40_000
    input:
        fa="work/{sample_name}/gfa2fa/{sample_name}.hic.p_ctg.fa",
    output:
        fai="work/{sample_name}/gfa2fa/{sample_name}.hic.p_ctg.fa.fai",
    log:
        "log/{sample_name}/samtools_index.log"
    shell:
        "samtools faidx --threads {threads} {input.fa} -o {output.fai} 2> {log}"


rule yahs:
    container:
        "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
    threads: 16
    resources:
        mem_mb=400_000
    input:
        bam="work/{sample_name}/hicpro/hicpro/mapping/{sample_name}_0_bwt2pairs.bam",
        fa="work/{sample_name}/gfa2fa/{sample_name}.hic.p_ctg.fa",
        fai="work/{sample_name}/gfa2fa/{sample_name}.hic.p_ctg.fa.fai",
    output:
        agp="work/{sample_name}/yahs/{sample_name}_scaffolds_final.agp",
        fa="work/{sample_name}/yahs/{sample_name}_scaffolds_final.fa",
        bin="work/{sample_name}/yahs/{sample_name}.bin",
    log:
        "logs/{sample_name}/yahs.log"
    shell:
        "yahs -o work/{wildcards.sample_name}/yahs/{wildcards.sample_name} {input.fa} {input.bam} 2> {log}"


rule output_scaffolds:
    input:
        fa="work/{sample_name}/yahs/{sample_name}_scaffolds_final.fa",
    output:
        fa="output/{sample_name}/{sample_name}_scaffolds_final.fa",
    log:
        "logs/{sample_name}/output_scaffolds.log"
    shell:
        "mv {input.fa} {output.fa} 2> {log}"


rule juicer_generate_JBAT:
    container:
        "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
    threads: 8
    resources:
        mem_mb=42_000
    input:
        agp="work/{sample_name}/yahs/{sample_name}_scaffolds_final.agp",
        bam="work/{sample_name}/hicpro/hicpro/mapping/{sample_name}_0_bwt2pairs.bam",
        fai="work/{sample_name}/gfa2fa/{sample_name}.hic.p_ctg.fa.fai",
    output:
        txt="output/{sample_name}/{sample_name}_JBAT.txt",
        assembly="output/{sample_name}/{sample_name}_JBAT.assembly",
        jbat_log="work/{sample_name}/JBAT/{sample_name}_JBAT.log"
    log:
        "logs/{sample_name}/juicer_generate_JBAT.log"
    shell:
        "juicer pre -a -o output/{wildcards.sample_name}/{wildcards.sample_name}_JBAT {input.bam} {input.agp} {input.fai} > {log} 2> {output.jbat_log}"


rule download_juicer_tools:
    output:
        juicer_tools="downloads/juicer_tools.1.9.9_jcuda.0.8.jar",
    log:
        "logs/downloads/get_juicer_tools.log"
    shell:
        "wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar -O {output.juicer_tools} 2> {log}"


rule juicer_tools_generate_hic:
    container:
        "oras://community.wave.seqera.io/library/java-jdk:8.0.112--c971247410d7cba8"
    threads: 8
    resources:
        mem_mb=33_000
    input:
        juicer_tools="downloads/juicer_tools.1.9.9_jcuda.0.8.jar",
        txt="output/{sample_name}/{sample_name}_JBAT.txt",
        jbat_log="work/{sample_name}/JBAT/{sample_name}_JBAT.log",
    output:
        hic="output/{sample_name}/{sample_name}_JBAT.hic"
    log:
        "logs/{sample_name}/juicer_tools_generate_hic.log"
    shell:
        "(java -jar -Xmx{resources.mem_mb}m {input.juicer_tools} pre {input.txt} {output.hic} <(cat {input.jbat_log} | grep PRE_C_SIZE | awk '{{print $2\" \"$3}}')) 2> {log}"
