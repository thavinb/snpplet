/*
 * Process definitions
 */


/*
 * Step 1a: Create a FASTA genome index (.fai) with samtools for GATK
 */
process PREPARE_GENOME_SAMTOOLS { 
  tag "$genome.baseName"
 
  input: 
    path genome
 
  output: 
    path "${genome}.fai"
  
  script:
  """
  samtools faidx $genome
  """
}


/*
 * Step 1b: Create a FASTA genome index with bwa for bwa
 */
process PREPARE_GENOME_BWA { 
  tag "$genome.baseName"
 
  input: 
    path genome
 
  output: 
    path "${genome}.*"
  
  script:
  """
  bwa index -a is $genome
  """
}


/*
 * Step 1c: Create a FASTA genome sequence dictionary with Picard for GATK
 */
process PREPARE_GENOME_PICARD {
  tag "$genome.baseName"

  input:
    path genome

  output:
    path "${genome.baseName}.dict"

  script:
  """
  gatk CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
  """
}


/*
 * Step 2a: Read QC before trimming using fastqc
 */
process FASTQC_BEFORE_TRIM {
  tag "$id"
  publishDir "$params.results/fastqc_before", mode: 'copy'

  input:
    tuple val(id), path(reads)

  output:
    path "*.zip"

  """
  fastqc $reads --noextract --quiet
  """
}


/*
 * Step 2b: Summarise fastqc results using multiqc
 */
process MULTIQC_FASTQC_BEFORE_TRIM {
  publishDir "$params.results/multiqc_fastqc", mode: 'copy'
  errorStrategy 'ignore'

  input:
    path fastqc_before_all

  output:
    path "before.html"
    path "before_data"

  """
  multiqc ${fastqc_before_all} --interactive --filename before
  """
}


/*
 * Step 2c. Trim adapters and low quality reads
 */
process TRIM {
  tag "$id"
  publishDir "$params.results/trim", mode: 'copy'
  errorStrategy 'ignore'

  input:
    tuple val(id), path(reads)
    path adapter
    val trim_option

  output:
    tuple val(id), path("*_paired.fastq.gz")

  """
  trimmomatic PE -phred33 $reads \
    ${id}_1_paired.fastq.gz ${id}_1_unpaired.fastq.gz \
    ${id}_2_paired.fastq.gz ${id}_2_unpaired.fastq.gz \
    ILLUMINACLIP:$adapter/TruSeq3-PE-2.fa:2:30:10 \
    ILLUMINACLIP:$adapter/NexteraPE-PE.fa:2:30:10 \
    ${trim_option}
  """
}


/*
 * Step 2d: Read QC after trimming using fastqc
 */
process FASTQC_AFTER_TRIM {
  tag "$id"
  publishDir "$params.results/fastqc_after", mode: 'copy'

  input:
    tuple val(id), path(reads)

  output:
    path "*.zip"

  """
  fastqc $reads --noextract --quiet
  """
}


/*
 * Step 2e: Summarise fastqc results using multiqc
 */
process MULTIQC_FASTQC_AFTER_TRIM {
  publishDir "$params.results/multiqc_fastqc", mode: 'copy'
  errorStrategy 'ignore'

  input:
    path fastqc_after_all

  output:
    path "after.html"
    path "after_data"

  """
  multiqc ${fastqc_after_all} --interactive --filename after
  """
}


/*
 * Step 3a: Align reads to the reference genome
 */
process READ_MAPPING_BWA {
  tag "$id"
  publishDir "$params.results/bam", mode: 'copy', \
    saveAs: { filename -> "${id}_$filename" }

  input: 
    path genome
    val genome_name
    path index
    tuple val(id), path(reads)
    val bwa_option

  output: 
    tuple \
      val(id), \
      path("markdup.bam"), \
      path("markdup.bam.bai")
    path "${id}_coverage.txt"

  """
  # read mapping
  bwa mem -R "@RG\\tID:${id}\\tSM:${id}\\tPL:Illumina" \
    ${bwa_option} $genome $reads > sam

  # sort and BAM file
  samtools fixmate -O bam sam bam_fixmate
  samtools sort -O bam -o bam_sort bam_fixmate
  samtools index bam_sort
 
  # mask duplicate
  gatk --java-options "-Xmx1g" MarkDuplicates -I bam_sort \
    -O markdup.bam -M bam_markdup_metrics.txt
  
  # index BAM file
  samtools index markdup.bam

  # depth
  samtools coverage markdup.bam > ${id}_coverage.txt

  sed -i '' "s/${genome_name}/$id/" ${id}_coverage.txt  # for BSD sed
  # sed -i "s/${genome_name}/$id/" ${id}_coverage.txt  # for gnu sed
  
  # free up some disk space
  rm sam bam_fixmate bam_sort
  """
}


/*
 * Step 3b: Combine per-sample coverage info into a single file
 */
process COVERAGE_OUTPUT {
  publishDir "$params.results/bam", mode: 'copy'

  input:
    path coverage

  output:
    path "coverage.tsv"

  """
  awk '(NR == 1) || (FNR > 1)' $coverage > coverage.tsv
  """
}


/*
 * Step 4a: Call variants for each sample
 */
process CALL_VARIANTS {
  tag "$id"
  publishDir "$params.results/vcf", mode: 'copy'

  input:
    path genome
    path index
    path dict
    tuple val(id), path(bam), path(bai)
    val haplotypecaller_option
 
  output: 
    path "${id}.g.vcf.gz", emit: vcf
    path "${id}.g.vcf.gz.tbi", emit: vcf_tbi

  """ 
  gatk HaplotypeCaller \
    -R $genome -I $bam \
    -O ${id}.g.vcf.gz \
    -ERC GVCF \
    ${haplotypecaller_option}
  """
}


process CREATE_SAMPLE_MAP {
  publishDir "$baseDir/data", mode: 'copy'

  input:
    path vcf_all

  output:
    path "sample_map.txt", emit: sample_map

  script:
  def vcf_all_list = vcf_all.collect{ "$it\t$it" }.join('\n')
  """
  echo "${vcf_all_list}" > out1.txt
  sed 's/.g.vcf.gz\t/\t/g' out1.txt > sample_map.txt
  """
}


/*
 * Step 4b: Call variants for all samples
 */
process JOINT_GENOTYPING {
  echo true
  publishDir "$params.results/vcf_all", mode: 'copy'

  input:
    path genome
    path index
    path dict
    val genome_name
    path vcf_all
    path vcf_tbi_all
    path sample_map_usr
    val genomicsdbimport_option
    val genotypegvcfs_option

  output:
    tuple path("merge.vcf.gz"), path("merge.vcf.gz.tbi")
  
  when:
    !params.stop

  script:
  def vcf_all_list = vcf_all.collect{ "-V $it" }.join(' ')
  def is_sample_map_usr = sample_map_usr.exists()

  if (is_sample_map_usr) {
    """
    gatk GenomicsDBImport \
      -R $genome \
      --sample-name-map ${sample_map_usr} \
      --genomicsdb-workspace-path dbcohort \
      -L $genome_name \
      ${genomicsdbimport_option}

    gatk GenotypeGVCFs -R $genome -V gendb://dbcohort \
      -O merge.vcf.gz ${genotypegvcfs_option}
    """
  } else {
    """
    gatk GenomicsDBImport \
      -R $genome \
      ${vcf_all_list} \
      --genomicsdb-workspace-path dbcohort \
      -L $genome_name \
      ${genomicsdbimport_option}

    gatk GenotypeGVCFs -R $genome -V gendb://dbcohort \
      -O merge.vcf.gz ${genotypegvcfs_option}
    """
  }
}


/*
 * Step 4c: Filter variants
 */
process FILTER_VARIANTS {
  publishDir "$params.results/vcf_all", mode: 'copy'

  input:
    path genome
    path index
    path dict
    val variant_filter
    val variant_filter_name
    path mask_positions
    tuple path("merge.vcf.gz"), path("merge.vcf.gz.tbi")
    val selectvariant_option
  
  output:
    tuple path("merge_filtered.vcf.gz"), path("merge_filtered.vcf.gz.tbi")
    path "merge_filtered_stats.txt"

  when:
    !params.stop

  """
  gatk VariantFiltration \
    -R $genome \
    -V merge.vcf.gz \
    -O merge_filt.vcf.gz \
    -filter "${variant_filter}" \
    --filter-name ${variant_filter_name}
  
  gatk SelectVariants \
    -V merge_filt.vcf.gz \
    -O merge_filtered.vcf.gz \
    ${selectvariant_option}

  bcftools stats -F $genome merge_filtered.vcf.gz > merge_filtered_stats.txt
  """
}


/*
 * Step 5: Convert filtered vcf to multiple sequence alignment
 */
process VCF_TO_FASTA {
  publishDir "$params.results/aln", mode: 'copy'

  input:
    tuple path(vcf), path(tbi)

  output:
    path "pos.txt"
    path "aln.fasta"

  when:
    !params.stop
  
  script:
  $/
  gatk VariantsToTable -V $vcf -O gt -F POS -GF GT

  # transpose genotype table from gatk
  datamash transpose --output-delimiter=, < gt > tmp

  # create alignment file
  sed -i '' '1d' tmp  # for BSD sed
  # sed -i '1d' tmp  # for gnu sed

  sed $'s/^/>/;s/.GT,/\\\n/g;s/,//g;s/[\.\*]/-/g' tmp > aln.fasta

  # cleanup
  rm tmp gt

  # optional: positions with ref variants
  gatk VariantsToTable -V $vcf -O pos.txt -F POS -F REF
  /$
}
