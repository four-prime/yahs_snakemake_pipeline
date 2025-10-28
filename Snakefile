import glob

configfile: "config.json"

sample_names = config["sample_names"]


def get_sample_input_paths():
    sample_hifi_paths = {}
    sample_hic_r1_paths = {}
    sample_hic_r2_paths = {}
    
    for sample in sample_names:

        valid_extensions = ('.fa', '.fasta', '.fq', '.fastq', '.fa.gz', '.fasta.gz', '.fq.gz', '.fastq.gz')

        hifi_path_list = glob.glob(f"samples/{sample}/input/hifi/*")
        hifi_path_list = [f for f in hifi_path_list if f.endswith(valid_extensions)]

        hic_r1_list = glob.glob(f"samples/{sample}/input/hic/r1/*")
        hic_r1_list  = [f for f in hic_r1_list  if f.endswith(valid_extensions)]

        hic_r2_list = glob.glob(f"samples/{sample}/input/hic/r2/*")
        hic_r2_list = [f for f in hic_r2_list if f.endswith(valid_extensions)]

        hifi_count = len(hifi_path_list)
        r1_count = len(hic_r1_list)
        r2_count = len(hic_r2_list)

        # Check for surplus Hi-C inputs
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

        # Check for missing inputs
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
   
    print(sample_hifi_paths, sample_hic_r1_paths, sample_hic_r2_paths)

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
        [f"samples/{sample_name}/output/out_JBAT.hic" for sample_name in sample_names],
        [f"samples/{sample_name}/output/out_JBAT.txt" for sample_name in sample_names],
        [f"samples/{sample_name}/output/out_JBAT.assembly" for sample_name in sample_names],
        [f"samples/{sample_name}/output/out_scaffolds_final.fa" for sample_name in sample_names]
        


rule generate_samplesheet:
    input:
        hic_1=get_hic_r1_path,
        hic_2=get_hic_r2_path,
    output:
        samplesheet="samples/{sample_name}/work/hic_samplesheet/samplesheet.csv",
    shell:
        "echo -e 'sample,fastq_1,fastq_2\\n{wildcards.sample_name},{input.hic_1},{input.hic_2}' > {output.samplesheet}"


rule hifiasm:
  container:
    "oras://community.wave.seqera.io/library/hifiasm:0.25.0--bcfc60a944a26aaa"
  threads: 32
  resources:
    mem_mb=1_200_000
  input:
    hifi=get_hifi_path,
    hic_1=get_hic_r1_path,
    hic_2=get_hic_r2_path,
  output:
    unitigs="samples/{sample_name}/work/hifiasm/out.hic.p_utg.gfa",
    contigs="samples/{sample_name}/work/hifiasm/out.hic.p_ctg.gfa",
  log:
    "samples/{sample_name}/logs/hifiasm.log"
  shell:
    "hifiasm -t {threads} -o samples/{wildcards.sample_name}/work/hifiasm/out -l0 -f0 --h1 {input.hic_1} --h2 {input.hic_2} {input.hifi} 2> {log}"


rule gfa2fa:
  threads: 8
  resources:
    mem_mb=50_000
  container:
    "oras://community.wave.seqera.io/library/gfatools_tabix:b5cb1ec01cf54783"
  input:
    fa="samples/{sample_name}/work/hifiasm/out.hic.p_ctg.gfa"
  output:
    fa="samples/{sample_name}/work/hifiasm/out.hic.p_ctg.fa"
  log:
    "samples/{sample_name}/logs/gfa2fa.log"
  shell:
    "gfatools gfa2fa {input.fa} > {output.fa} 2> {log}"


rule hic_pro:
    handover: True
    input:
        samplesheet="samples/{sample_name}/work/hic_samplesheet/samplesheet.csv",
        fa="samples/{sample_name}/work/hifiasm/out.hic.p_ctg.fa",
        config="hic_pro.config"
    output:
        bam="samples/{sample_name}/work/hicpro/hicpro/mapping/{sample_name}_0_bwt2pairs.bam",
    shell:
        "nextflow run nf-core/hic -r 2.1.0 "
        "-work-dir samples/{wildcards.sample_name}/work/hicpro_work "
        "--outdir samples/{wildcards.sample_name}/work/hicpro "
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
        fa="samples/{sample_name}/work/hifiasm/out.hic.p_ctg.fa",
    output:
        fai="samples/{sample_name}/work/hifiasm/out.hic.p_ctg.fa.fai",
    log:
        "samples/{sample_name}/logs/samtools_index.log"
    shell:
        "samtools faidx {input.fa} -o {output.fai} 2> {log}"


rule yahs:
    container:
        "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
    threads: 16
    resources:
        mem_mb=400_000
    input:
        bam="samples/{sample_name}/work/hicpro/hicpro/mapping/{sample_name}_0_bwt2pairs.bam",
        fa="samples/{sample_name}/work/hifiasm/out.hic.p_ctg.fa",
        fai="samples/{sample_name}/work/hifiasm/out.hic.p_ctg.fa.fai",
    output:
        agp="samples/{sample_name}/work/yahs/out_scaffolds_final.agp",
        fa="samples/{sample_name}/work/yahs/out_scaffolds_final.fa",
        bin="samples/{sample_name}/work/yahs/out.bin",
    log:
        "samples/{sample_name}/logs/yahs.log"
    shell:
        "yahs -o samples/{wildcards.sample_name}/work/yahs/out {input.fa} {input.bam} 2> {log}"


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
        mem_mb=32_000
    input:
        agp="samples/{sample_name}/work/yahs/out_scaffolds_final.agp",
        bam="samples/{sample_name}/work/hicpro/hicpro/mapping/{sample_name}_0_bwt2pairs.bam",
        fai="samples/{sample_name}/work/hifiasm/out.hic.p_ctg.fa.fai",
    output:
        txt="samples/{sample_name}/output/out_JBAT.txt",
        assembly="samples/{sample_name}/output/out_JBAT.assembly",
        jbat_log="samples/{sample_name}/work/out_JBAT.log"
    log:
        "samples/{sample_name}/logs/juicer_generate_JBAT.log"
    shell:
        "juicer pre -a -o samples/{wildcards.sample_name}/work/out_JBAT {input.bam} {input.agp} {input.fai} > {log} 2> {output.jbat_log}"


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
        mem_mb=33_000
    input:
        juicer_tools="downloads/juicer_tools.1.9.9_jcuda.0.8.jar",
        txt="samples/{sample_name}/work/out_JBAT.txt",
        jbat_log="samples/{sample_name}/work/out_JBAT.log",
    output:
        hic="samples/{sample_name}/output/out_JBAT.hic"
    log:
        "samples/{sample_name}/logs/juicer_tools_generate_hic.log"
    shell:
        "(java -jar -Xmx{resources.mem_mb}m {input.juicer_tools} pre {input.txt} {output.hic} <(cat {input.jbat_log} | grep PRE_C_SIZE | awk '{{print $2\" \"$3}}')) 2> {log}"
