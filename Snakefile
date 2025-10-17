rule samtools_faidx_source:
  container:
    "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"
  threads: 8
  resources:
    mem_mb=20000
  input:
    "{sample_name}/{source}.fa", 
  output:
    "{sample_name}/idx_{source}.fa.fai", 
  log:
    "{sample_name}/{source}_samtools_faidx.log"
  shell:
    "samtools faidx {input} -o {output} 2> {log}"

rule sambamba_deduplicate:
  container:
    "oras://community.wave.seqera.io/library/sambamba:1.0.1--69782c710ee49140"
  threads: 8
  resources:
    mem_mb=60000
  input:
     "{sample_name}/{source}.bam",
  output:
     "{sample_name}/deduplicated_{source}.bam",
  log:
    "{sample_name}/{source}_sambamba.log"
  shell:
    "sambamba markdup -t {threads} --show-progress {input} {output}"

rule yahs:
  container:
    "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
  threads: 16
  resources:
    mem_mb=200000
  input:
     bam="{sample_name}/deduplicated_{source}.bam",
     fasta="{sample_name}/{source}.fa",
     fai="{sample_name}/idx_{source}.fa.fai"
  output:
    agp="{sample_name}/{source}/out_scaffolds_final.agp",
    fa="{sample_name}/{source}/out_scaffolds_final.fa",
    bin="{sample_name}/{source}/out.bin",
  log:
    "{sample_name}/{source}_yahs.log"
  shell:
    "yahs -o {wildcards.sample_name}/{wildcards.source}/out {input.fasta} {input.bam} 2> {log}"

rule get_chrom_sizes:
  container:
    "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"
  threads: 8
  resources:
    mem_mb=20000
  input:
    "{sample_name}/{source}/out_scaffolds_final.fa",
  output:
    # "{sample_name}/{source}/out_scaffolds_final.fa.fai",
    "{sample_name}/{source}_scaffolds_final.chrom.sizes"
  log:
    "{sample_name}/{source}_samtools_faidx.log"
  shell:
    "samtools faidx {input} | cut -f 1,2 > {output}"

# rule juicer_pre:
#   container:
#     "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
#   threads: 8
#   resources:
#     mem_mb=32000
#   input:
#     agp="{sample_name}/{source}/out_scaffolds_final.agp",
#     bin="{sample_name}/{source}/out.bin",
#     fai="{sample_name}/idx_{source}.fa.fai",
#     # fai="{sample_name}/{source}/out_scaffolds_final.fa.fai",
#   output:
#     alignments="{sample_name}/{source}_alignments_sorted.txt",
#   log:
#     "{sample_name}/{source}_juicer_pre.log"
#   shell:
#     "(juicer pre {input.bin} {input.agp} {input.fai} | "
#     "sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | "
#     "awk 'NF' > {wildcards.sample_name}/alignments_sorted.txt.part) && "
#     "(mv {wildcards.sample_name}/alignments_sorted.txt.part {output.alignments})"

# rule juicer_tools_generate_hic_maps:
#   container:
#     # NOTE: this could powwibly be replaced with the sequera container for juicer tools directly, once we determine why this isn't working 
#     "oras://community.wave.seqera.io/library/java-jdk:8.0.112--c971247410d7cba8"
#   threads: 8
#   resources:
#     mem_mb=33000
#   input:
#     juicer_tools="utils/juicer_tools.1.9.9_jcuda.0.8.jar",
#     alignments="{sample_name}/{source}_alignments_sorted.txt",
#     chrom_sizes="{sample_name}/{source}_scaffolds_final.chrom.sizes"
#   output:
#     "{sample_name}/{source}.hic"
#   log:
#     "{sample_name}/{source}_juicer_generate_hic_maps.log"
#   shell:
#     "java -Xmx32G -jar {input.juicer_tools} pre {input.alignments} {output} {input.chrom_sizes} 2> {log}"
#     # "juicer_tools

rule juicer_generate_JBAT:
  container:
    "oras://community.wave.seqera.io/library/yahs:1.2a.2--bd14e0d1b929fa78"
  threads: 8
  resources:
    mem_mb=32000
  input:
    agp="{sample_name}/{source}/out_scaffolds_final.agp",
    # Can replace the bin with the original bed, should try this
    # bin="{sample_name}/{source}/out.bin",
    bam="{sample_name}/deduplicated_{source}.bam",
    fai="{sample_name}/idx_{source}.fa.fai",
    # fai="{sample_name}/{source}/out_scaffolds_final.fa.fai",
  output:
    txt="{sample_name}/{source}_JBAT.txt",
    assembly="{sample_name}/{source}_JBAT.assembly",
    jbat_log="{sample_name}/{source}_JBAT.log"
  log:
    "{sample_name}/{source}_juicer_pre.log"
  shell:
    "juicer pre -a -o {wildcards.sample_name}/{wildcards.source}_JBAT {input.bam} {input.agp} {input.fai} > {log} 2> {output.jbat_log}"

rule juicer_tools_generate_hic:
  container:
    # "oras://community.wave.seqera.io/library/java-jdk:8.0.112--c971247410d7cba8"
    "oras://community.wave.seqera.io/library/juicertools:2.20.00--ad76f0f8995fd430"
  threads: 8
  resources:
    mem_mb=33000
  input:
    # juicer_tools="utils/juicer_tools.1.9.9_jcuda.0.8.jar",
    txt="{sample_name}/{source}_JBAT.txt",
    jbat_log="{sample_name}/{source}_JBAT.log"
  output:
    "{sample_name}/{source}_JBAT.hic"
  log:
    "{sample_name}/{source}_juicer_tools.log"
  shell:
    # "(java -jar -Xmx32G {input.juicer_tools} pre {input.txt} {output} <(cat {input.jbat_log} | grep PRE_C_SIZE | awk '{{print $2\" \"$3}}')) 2> {log}"
    "(juicer_tools pre {input.txt} {output} <(cat {input.jbat_log} | grep PRE_C_SIZE | awk '{{print $2\" \"$3}}')) 2> {log}"
    # The above line results in the following error:
    # Exception in thread "main" java.lang.OutOfMemoryError: Java heap space
  	# at java.base/java.io.FilterOutputStream.<init>(FilterOutputStream.java:60)
  	# at htsjdk.tribble.util.LittleEndianOutputStream.<init>(LittleEndianOutputStream.java:25)
  	# at juicebox.tools.utils.original.MatrixZoomDataPP.dumpBlocks(MatrixZoomDataPP.java:361)
  	# at juicebox.tools.utils.original.MatrixZoomDataPP.incrementCount(MatrixZoomDataPP.java:245)
  	# at juicebox.tools.utils.original.MatrixPP.incrementCount(MatrixPP.java:148)
  	# at juicebox.tools.utils.original.Preprocessor.writeBody(Preprocessor.java:759)
  	# at juicebox.tools.utils.original.Preprocessor.preprocess(Preprocessor.java:452)
  	# at juicebox.tools.clt.old.PreProcessing.run(PreProcessing.java:176)
  	# at juicebox.tools.HiCTools.main(HiCTools.java:97)
  	#
  	# Should look into allocating more threads and memory

    # "juicer_tools pre {input.alignments} out.hic.part {input.chrom_sizes}) && (mv out.hic.part {output}) 2> {log}"
