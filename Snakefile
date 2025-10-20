rule all:
  input:
    "stipa/output/out_JBAT.hic"

rule samtools_faidx_source:
  container:
    "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"
  threads: 8
  resources:
    mem_mb=20000
  input:
    "{sample_name}/input/filename.fa", 
  output:
    "{sample_name}/input/filename.fa.fai", 
  log:
    "{sample_name}/logs/samtools_index_input.log"
  shell:
    "samtools faidx {input} -o {output} 2> {log}"

rule sambamba_deduplicate:
  container:
    "oras://community.wave.seqera.io/library/sambamba:1.0.1--69782c710ee49140"
  threads: 8
  resources:
    mem_mb=60000
  input:
     "{sample_name}/input/filename.bam",
  output:
     "{sample_name}/work/sambamba_deduplicate/filename.bam",
  log:
    "{sample_name}/logs/sambamba_deduplicate.log"
  shell:
    "sambamba markdup -t {threads} --show-progress {input} {output}"

rule yahs:
  container:
    "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
  threads: 16
  resources:
    mem_mb=200000
  input:
    bam="{sample_name}/work/sambamba_deduplicate/filename.bam",
    fasta="{sample_name}/input/filename.fa",
    fai="{sample_name}/input/filename.fa.fai", 
  output:
    agp="{sample_name}/work/yahs/out_scaffolds_final.agp",
    fa="{sample_name}/work/yahs/out_scaffolds_final.fa",
    bin="{sample_name}/work/yahs/out.bin",
  log:
    "{sample_name}/logs/yahs.log"
  shell:
    "yahs -o {wildcards.sample_name}/work/yahs/out {input.fasta} {input.bam} 2> {log}"

rule juicer_generate_JBAT:
  container:
    "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
  threads: 8
  resources:
    mem_mb=32000
  input:
    agp="{sample_name}/work/yahs/out_scaffolds_final.agp",
    bam="{sample_name}/work/sambamba_deduplicate/filename.bam",
    fai="{sample_name}/input/filename.fa.fai", 
  output:
    txt="{sample_name}/work/out_JBAT.txt",
    assembly="{sample_name}/work/out_JBAT.assembly",
    jbat_log="{sample_name}/work/out_JBAT.log"
  log:
    "{sample_name}/logs/juicer_generate_JBAT.log"
  shell:
    "juicer pre -a -o {wildcards.sample_name}/work/out_JBAT {input.bam} {input.agp} {input.fai} > {log} 2> {output.jbat_log}"

rule download_juicer_tools:
  output:
    juicer_tools="utils/juicer_tools.1.9.9_jcuda.0.8.jar",
  log:
    "utils/logs/get_juicer_tools.log"
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
    txt="{sample_name}/work/out_JBAT.txt",
    jbat_log="{sample_name}/work/out_JBAT.log"
  output:
    "{sample_name}/output/out_JBAT.hic"
  log:
    "{sample_name}/logs/juicer_tools_generate_hic.log"
  shell:
    "(java -jar -Xmx{resources.mem_mb}m {input.juicer_tools} pre {input.txt} {output} <(cat {input.jbat_log} | grep PRE_C_SIZE | awk '{{print $2\" \"$3}}')) 2> {log}"

