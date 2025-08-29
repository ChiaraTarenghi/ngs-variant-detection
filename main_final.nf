nextflow.enable.dsl=2

// ---- Defaults to silence WARNs ----
params.reads      = params.reads      ?: 'data/*_{R1,R2}.fastq.gz'
params.fasta      = params.fasta      ?: 'ref_demo.fasta'
params.outdir     = params.outdir     ?: 'results'
params.use_snpeff = params.use_snpeff ?: false
params.use_vep    = params.use_vep    ?: false
params.snpeff_db  = params.snpeff_db  ?: 'GRCh38.99'
params.snpeff_data= params.snpeff_data?: 'snpeff_data'
params.vep_cache_dir = params.vep_cache_dir ?: 'vep_cache'
params.vep_species   = params.vep_species   ?: 'homo_sapiens'
params.vep_assembly  = params.vep_assembly  ?: 'GRCh38'
params.min_qual   = (params.min_qual ?: 20)    as Integer
params.min_dp     = (params.min_dp   ?: 1)     as Integer
params.max_dp     = (params.max_dp   ?: 10000) as Integer
params.min_af     = (params.min_af   ?: 0.0)   as Float
params.bwa_opts   = params.bwa_opts  ?: '-k 5 -T 0'


log.info "use_snpeff=${params.use_snpeff} use_vep=${params.use_vep} bwa_opts='${params.bwa_opts}'"

// -------------------- PROCESSES -----------------

process p_fastqc {
  tag { sample }
  container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
  publishDir "${params.outdir}/fastqc", mode: 'copy'
  input:
    tuple val(sample), path(r1), path(r2)
  output:
    path("*_fastqc.zip")
  script:
  """
  set -euo pipefail
  fastqc -q -o . "${r1}" "${r2}"
  """
}

process p_bwa_mem2 {
  tag { sample }
  container 'quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_2'
  publishDir "${params.outdir}", mode: 'copy'
  input:
    tuple val(sample), path(r1), path(r2)
    path fasta
  output:
    tuple val(sample), path("${sample}.sam")
  script:
  """
  set -euo pipefail
  bwa-mem2 index "${fasta}" || true
  RG=\$(printf '@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA\\tLB:lib1\\tPU:%s.1' "${sample}" "${sample}" "${sample}")
  bwa-mem2 mem -t ${task.cpus ?: 2} ${params.bwa_opts} -R "\$RG" "${fasta}" "${r1}" "${r2}" > "${sample}.sam"
  """
}


process p_samtools_sort {
  tag { sample }
  container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'
  publishDir "${params.outdir}", mode: 'copy'
  input:
    tuple val(sample), path(sam)
  output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")
  script:
  """
  set -euo pipefail
  samtools sort -@ ${task.cpus ?: 2} -o "${sample}.sorted.bam" "${sam}"
  samtools index "${sample}.sorted.bam"
  """
}

process p_markdup {
  tag { sample }
  container 'quay.io/biocontainers/picard:3.2.0--hdfd78af_0'
  publishDir "${params.outdir}", mode: 'copy'
  input:
    tuple val(sample), path(bam_sorted), path(bai)
  output:
    tuple val(sample), path("${sample}.marked.bam"), path("${sample}.marked.bam.bai")
  script:
  """
  set -euo pipefail
  picard MarkDuplicates \
    I="${bam_sorted}" \
    O="${sample}.marked.bam" \
    M="${sample}.marked.metrics.txt" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

  # Build index if missing
  if [ ! -s "${sample}.marked.bam.bai" ] && [ ! -s "${sample}.marked.bai" ]; then
    picard BuildBamIndex I="${sample}.marked.bam"
  fi

  # Normalize Picard's .bai to .bam.bai if needed
  if [ -s "${sample}.marked.bai" ] && [ ! -s "${sample}.marked.bam.bai" ]; then
    mv "${sample}.marked.bai" "${sample}.marked.bam.bai"
  fi
  """
}

process p_faidx {
  container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'
  publishDir "${params.outdir}", mode: 'copy'
  input:
    path fasta
  output:
    tuple path(fasta), path("${fasta}.fai")
  script:
  """
  set -euo pipefail
  samtools faidx "${fasta}"
  """
}

process p_bcftools_call {
  tag { sample }
  container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  publishDir "${params.outdir}/vcf", mode: 'copy'
  input:
    tuple val(sample), path(bam), path(bai)
    path fasta
    path fai
  output:
  tuple val(sample), path("${sample}.raw.vcf")

  script:
  """
  set -euo pipefail
  bcftools mpileup \
    -f "${fasta}" \
    -q 0 -Q 0 \
    -B -A \
    -a DP,AD \
    -Ou "${bam}" \
  | bcftools call \
      -c -v \
      --ploidy 2 \
      --pval-threshold 1.0 \
      -Ov -o "${sample}.raw.vcf"
  """
}


process p_bcftools_norm {
  tag { sample }
  container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  publishDir "${params.outdir}/vcf", mode: 'copy'
  input:
    tuple val(sample), path(vcf_in)
    path fasta
    path fai
  output:
    tuple val(sample), path("${sample}.norm.vcf")
  script:
  """
  set -euo pipefail
  # Left-align and split multiallelic records
  bcftools norm -f "${fasta}" -m -both "${vcf_in}" -Ov -o "${sample}.norm.vcf"
  """
}

process p_bcftools_filter {
  tag { sample }
  container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  publishDir "${params.outdir}/vcf", mode: 'copy'
  input:
    tuple val(sample), path(vcf_in)
  output:
    tuple val(sample), path("${sample}.filtered.vcf")
  script:
  """
  set -euo pipefail
  bcftools filter -i "QUAL>=${params.min_qual} && INFO/DP>=${params.min_dp} && INFO/DP<=${params.max_dp} && (INFO/AF>=${params.min_af} || AF>=${params.min_af})" \
    "${vcf_in}" -Ov -o "${sample}.filtered.vcf"
  """
}


process p_bcftools_annotate {
  tag { sample }
  container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  publishDir "${params.outdir}/vcf", mode: 'copy'
  input:
    tuple val(sample), path(vcf_in)
  output:
    tuple val(sample), path("${sample}.annotated.vcf")
  script:
  """
  set -euo pipefail
  bcftools +fill-tags "${vcf_in}" -Ov -o "${sample}.annotated.vcf" -- -t AC,AF,AN
  """
}

process p_bcftools_stats {
  tag { sample }
  container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  publishDir "${params.outdir}/vcf", mode: 'copy'
  input:
    tuple val(sample), path(vcf_in)
  output:
    tuple val(sample), path("${sample}.bcfstats.txt")
  script:
  """
  set -euo pipefail
  bcftools stats "${vcf_in}" > "${sample}.bcfstats.txt"
  """
}

process p_snpeff {
  when:
    params.use_snpeff && !params.use_vep
  tag { sample }
  container 'quay.io/biocontainers/snpeff:5.2--hdfd78af_0'
  publishDir "${params.outdir}/vcf", mode: 'copy'
  input:
    tuple val(sample), path(vcf_in)
  output:
    tuple val(sample), path("${sample}.snpeff.vcf")
  script:
  """
  set -euo pipefail
  SnpEff -dataDir "${params.snpeff_data}" -v "${params.snpeff_db}" "${vcf_in}" > "${sample}.snpeff.vcf"
  """
}

process p_vep {
  when:
    params.use_vep && !params.use_snpeff
  tag { sample }
  container 'quay.io/biocontainers/ensembl-vep:110.1--pl5321hdfd78af_1'
  publishDir "${params.outdir}/vcf", mode: 'copy'
  input:
    tuple val(sample), path(vcf_in)
  output:
    tuple val(sample), path("${sample}.vep.vcf")
  script:
  """
  set -euo pipefail
  vep \
    --offline --cache \
    --dir_cache "${params.vep_cache_dir}" \
    --species "${params.vep_species}" \
    --assembly "${params.vep_assembly}" \
    --vcf --no_stats \
    --input_file "${vcf_in}" \
    --output_file "${sample}.vep.vcf"
  """
}

process p_multiqc {
  container 'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0'
  publishDir "${params.outdir}/multiqc", mode: 'copy'
  input:
    path(anything)
  output:
    path("multiqc_report.html")
  script:
  """
  set -euo pipefail
  multiqc -q -o . .
  """
}

// -------------------- WORKFLOW ------------------

workflow {

  // Pair R1/R2
  reads_ch = Channel
              .fromFilePairs(params.reads, flat: true)
              .map { id, r1, r2 -> tuple(id as String, r1, r2) }

  fasta_ch = Channel.value( file(params.fasta) )

  // QC
  p_fastqc( reads_ch )

  // Map -> Sort -> MarkDuplicates
  p_bwa_mem2( reads_ch, fasta_ch )
  p_samtools_sort( p_bwa_mem2.out )
  p_markdup( p_samtools_sort.out )

  // Reference index
  p_faidx( fasta_ch )

  // Calling -> Norm -> Filter -> Annotate
  p_bcftools_call(
    p_markdup.out,                 // (sample, marked.bam, marked.bam.bai)
    p_faidx.out.map{ it[0] },      // fasta
    p_faidx.out.map{ it[1] }       // fai
  )
  p_bcftools_norm(
    p_bcftools_call.out,
    p_faidx.out.map{ it[0] },
    p_faidx.out.map{ it[1] }
  )
  p_bcftools_filter(   p_bcftools_norm.out )
  p_bcftools_annotate( p_bcftools_filter.out )

  // (Opzionali) snpEff / VEP - abilita esattamente uno
  if( params.use_snpeff && !params.use_vep )
    p_snpeff( p_bcftools_filter.out )
  if( params.use_vep && !params.use_snpeff )
    p_vep( p_bcftools_filter.out )

  // Variant stats -> una sola chiamata
  p_bcftools_stats( p_bcftools_filter.out )

  // MultiQC -> una sola chiamata, consuma FastQC + bcftools stats
  def mqc_inputs = p_fastqc.out.mix( p_bcftools_stats.out.map{ it[1] } )
  p_multiqc( mqc_inputs )
}

