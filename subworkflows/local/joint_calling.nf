//
// Joint calling
//

include { GATK4_GENOMICSDBIMPORT      } from '../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS         } from '../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_VARIANTFILTRATION     } from '../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_SELECTVARIANTS        } from '../../modules/nf-core/gatk4/selectvariants/main'
include { BCFTOOLS_STATS              } from '../../modules/nf-core/bcftools/stats/main'

workflow JOINT_CALLING {
    take:
    vcf              // channel: [ path(vcf.collect()) ]
    tbi              // channel: [ path(tbi.collect()) ]
    fasta            // channel: [ val( meta ), path(fasta) ]
    fai              // channel: [ val( meta ), path(fai) ]
    dict             // channel: [ val( meta ), path(dict) ]
    interval         // channel: [ val( meta ), path(interval) ]
    sample_map       // channel: [ path sample_map ]


    main:
    ch_versions = Channel.empty()

    //
    // Switch for sample_map in GATK4_GENOMICSDBIMPORT
    //
    if (params.sample_map) {
        vcf.map {
            it -> { it.plus(0, [ sample_map ]) }
        }.set { ch_vcf }
        sample_map_flag = true
    } else {
        ch_vcf = vcf
        sample_map_flag = false
    }

    //
    // PREPARE INPUT: FOR JOINT_CALLING SUBWORKFLOW.
    // channel: [ val( meta ), path(vcf), path(tbi), path(interval_file), val(interval_value), path(wspace)) ]
    // or
    // channel: [ val( meta ), path(sample_map), path(tbi), path(interval_file), val(interval_value), path(wspace)) ]
    //

    ch_vcf
        .concat(tbi)
        .toList()
        .map { it -> tuple(id:"joint", it[0], it[1], interval, [], []) }
        .set { ch_dbimport }

    ch_dbimport.view()

    //
    // MODULE: Create GENOMICSDC from reference fasta and a collection of vcf file.
    // emit: genomicsdb
    // emit: updatedb
    // emit: intervallist
    //
    GATK4_GENOMICSDBIMPORT (
        ch_dbimport,
        false,
        false,
        sample_map_flag
    )
    ch_versions = ch_versions.mix( GATK4_GENOMICSDBIMPORT.out.versions )

    //
    // PREPARE INPUT FOR GENOTYPEGCVFS
    // channel: [ val(meta), path(gvcf), path(gvcf_index), path(intervals), path(intervals_index) ]
    //
    GATK4_GENOMICSDBIMPORT.out.genomicsdb
        .map { tuple(it[0], it[1], [], interval, []) }
        .set { ch_gatk4_genotypegvcfs }

    //
    // MODULE: Joint Genotype
    // emit: vcf
    // emit: tbi
    //
    GATK4_GENOTYPEGVCFS (
        ch_gatk4_genotypegvcfs,
        fasta[1],
        fai.map { it[1] },
        dict.map { it[1] },
        [],
        []
    )
    ch_versions = ch_versions.mix( GATK4_GENOTYPEGVCFS.out.versions )

    //
    // PREPARE INPUT FOR VARIANTFILTRATION
    // channel: [ val(meta), path(vcf), path(tbi) ]
    //

    GATK4_GENOTYPEGVCFS.out.vcf
        .join(GATK4_GENOTYPEGVCFS.out.tbi)
        .set { gvcf }

    //
    // MODULE: Mask filtered variant
    // emit: vcf
    // emit: tbi
    //
    GATK4_VARIANTFILTRATION (
       gvcf,
       fasta,
       fai,
       dict
    )
    ch_versions = ch_versions.mix( GATK4_VARIANTFILTRATION.out.versions )

    //
    // PREPARE INPUT FOR SELECTFILTRATION
    // channel: [ val(meta), path(vcf), path(tbi) ]
    //
    GATK4_VARIANTFILTRATION.out.vcf
        .join(GATK4_VARIANTFILTRATION.out.tbi)
        .map { tuple(it[0], it[1], it[2], interval) }
        .set { gvcf_filter }

    //
    // MODULE: Filter variant
    // emit: vcf
    // emit: tbi
    //
    GATK4_SELECTVARIANTS (
        gvcf_filter
    )
    ch_versions = ch_versions.mix( GATK4_SELECTVARIANTS.out.versions )

    //
    // PREPARE INPUT FOR BCFTOOLS
    // channel: [ val(meta), path(vcf), path(tbi) ]
    //

    GATK4_SELECTVARIANTS.out.vcf
        .join( GATK4_SELECTVARIANTS.out.tbi )
        .set { selected_gvcf }

    //
    // MODULE: Extract GVCF stats
    // emit: stats
    //
    BCFTOOLS_STATS (
        selected_gvcf,
        fasta,
        fai,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix( BCFTOOLS_STATS.out.versions )

    emit:
    vcf = GATK4_SELECTVARIANTS.out.vcf                              // channel: [ val(meta), vcf ]
    tbi = GATK4_SELECTVARIANTS.out.tbi                              // channel: [ val(meta), tbi ]
    stats = BCFTOOLS_STATS.out.stats                                // channel: [ val(meta), *stats.txt ]
    versions = ch_versions                                          // channel: [ versions.yml ]
}


