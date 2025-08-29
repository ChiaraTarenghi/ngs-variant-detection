nextflow.enable.dsl=2

params.fasta = 'ref.fasta'

process BWA_MEM2 {
  container "quay.io/biocontainers/bwa-mem2:2.2.1--he4a0461_2"
  input:
    tuple val(sample), path r1, path r2, path fasta
  output:
    tuple val(sample), path("${sample}.sam")
  script:
  """
  bwa-mem2 index ${fasta} || true
  bwa-mem2 mem -t 4 ${fasta} ${r1} ${r2} > ${sample}.sam
  """
}

workflow {
  Channel
    .of( tuple('S1', file('R1.fq.gz'), file('R2.fq.gz'), file(params.fasta)) )
    .set { one }
  BWA_MEM2(one)
}
