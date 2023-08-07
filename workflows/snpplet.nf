/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSnpplet.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.adapter, params.interval ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// TODO: get correct id name for igenome fasta
if (params.fasta) {
    ch_fasta = [  [ id:file(params.fasta).getBaseName() ] , file(params.fasta) ]
    } else if (params.genome.fasta) {
        ch_fasta = [ [ id:file(params.genome.fasta).getBaseName() ] , file(params.genome.fasta) ]
    } else {
        exit 1, 'Reference must be specified!! either locally (--fasta) or through aws igenomes (--genome)'
    }

if (params.adapter) { ch_adapter = file(params.adapter) } else { ch_adapter = [] }

def create_default_interval ( fasta , interval_file ) {
    interval_file.text = ''
    def is_header = { it[0].contains(">") }
    fasta_reader = file(fasta).newReader()
    while ( line = fasta_reader.readLine() ) {
        if ( line.any(is_header) ) {
        header_len = line.length()
        interval_list.append("${line[1..-1]}\n")
        }
    }
    print "[INFO]: Interval file does not specified. Default interval of entire fasta will be used"
    print "[INFO]: Default GATK-style `.list` interval file is created at `assets/default_interval.list`"
    print "[INFO]: To resume the pipeline correctly, please specify `--interval ./assets/default_interval.list` when using `-resume` flag in your next run."

    return interval_file
}

if (params.interval) {
    ch_interval = file(params.interval)
    extension = ch_interval.getExtension()
    if ( extension != "list" && extension != "intervals" && extension != "bed" && extension != "vcf" ) {
        exit 1 , "Interval file must have one of the supported file extensions ([.list, .intervals, .bed, .vcf])"
    }
} else if (params.fasta) {
    interval_list = file("$projectDir/assets/default_interval.list")
    ch_interval = create_default_interval( params.fasta, interval_list )
} else if (params.genome) {
    interval_list = file("$projectDir/assets/default_interval.list")
    ch_interval = create_default_interval( params.genome, interval_list )
}

if (params.sample_map) { sample_map = file(params.sample_map) } else { sample_map = [] }

ch_drgstrmodel = []
ch_dbsnp = []
ch_dbsnp_tbi = []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK   } from '../subworkflows/local/input_check'
include { READ_MAPPING  } from '../subworkflows/local/read_mapping'
include { JOINT_CALLING } from '../subworkflows/local/joint_calling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTP                          } from '../modules/nf-core/fastp/main'
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_FAIDX                 } from '../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_HAPLOTYPECALLER          } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_VARIANTSTOTABLE          } from '../modules/local/gatk4/variantstotable/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS    } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CREATE_MSA                     } from '../modules/local/create_msa'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SNPPLET {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: QC reads
    //
    FASTP (
        INPUT_CHECK.out.reads,
        ch_adapter,
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    //
    // SUBWORKFLOW: Read mapping to reference fasta
    //
    READ_MAPPING (
        FASTP.out.reads,
        ch_fasta
    )
    ch_versions = ch_versions.mix(READ_MAPPING.out.versions)

    //
    // MODULE: Create GATK dictionary from reference
    //
    GATK4_CREATESEQUENCEDICTIONARY( ch_fasta )
    ch_versions = ch_versions.mix( GATK4_CREATESEQUENCEDICTIONARY.out.versions )

    //
    // PREPARE INPUT: FOR VARIANT_CALLING MODULE
    // channel: [ val(meta), bam, bai, interval, drgstrmodel ]
    //
    READ_MAPPING.out.bam
        .map {
            [ it[0], it[1], it[2], ch_interval, ch_drgstrmodel ]
        }.set { ch_variantcalling_input }

    //
    // MODULE: Index reference fasta
    //
    SAMTOOLS_FAIDX ( ch_fasta )
    ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )

    //
    // MOUDLE: Variant Calling
    // TODO: make script to make sample-map of all vcf output.
    //
    GATK4_HAPLOTYPECALLER (
        ch_variantcalling_input,
        ch_fasta[1],
        SAMTOOLS_FAIDX.out.fai
            .map { it[1] },
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
            .map { it[1] },
        ch_dbsnp,
        ch_dbsnp_tbi
    )
    ch_versions = ch_versions.mix( GATK4_HAPLOTYPECALLER.out.versions )

    //
    // PREPARE INPUT: collect all vcf into a list to make a joint vcf
    //
    ch_vcf = GATK4_HAPLOTYPECALLER.out.vcf.collect { it[1] }
    ch_tbi = GATK4_HAPLOTYPECALLER.out.tbi.collect { it[1] }

    //
    // SUBWORKFLOW: Joint calling
    //
    JOINT_CALLING (
        ch_vcf,
        ch_tbi,
        ch_fasta,
        SAMTOOLS_FAIDX.out.fai,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict,
        ch_interval,
        sample_map
    )
    ch_versions = ch_versions.mix( JOINT_CALLING.out.versions )

    //
    // MODULE: Turn variants to tabular format
    //
    GATK4_VARIANTSTOTABLE (
       JOINT_CALLING.out.vcf,
       JOINT_CALLING.out.tbi,
    )
    ch_versions = ch_versions.mix( GATK4_VARIANTSTOTABLE.out.versions )

    CREATE_MSA (
        GATK4_VARIANTSTOTABLE.out.table
    )
    ch_versions = ch_versions.mix( CREATE_MSA.out.versions )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSnpplet.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSnpplet.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(JOINT_CALLING.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(READ_MAPPING.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(READ_MAPPING.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
