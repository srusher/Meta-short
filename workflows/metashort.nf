/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

//skip the line below to prevent "--fasta not specified" errors
//WorkflowMetashort.initialise(params, log)

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
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW                           } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED                       } from '../modules/nf-core/fastqc/main'
include { TRIMMOMATIC                                    } from '../modules/nf-core/trimmomatic/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_MAIN                } from '../modules/nf-core/kraken2/kraken2/main'
include { KRONA_KRONADB                                  } from '../modules/nf-core/krona/krona_db/main'
include { KRAKENTOOLS_KREPORT2KRONA                      } from '../modules/nf-core/krakentools/kreport2krona/main'
include { KRAKENTOOLS_EXTRACTKRAKENREADS                 } from '../modules/nf-core/krakentools/extractkrakenreads'
include { KRONA_KTIMPORTTEXT                             } from '../modules/nf-core/krona/ktimporttext/main'
include { METAPHLAN                                      } from '../modules/nf-core/metaphlan/main'
include { SPADES                                         } from '../modules/nf-core/spades/main'
include { MEGAHIT                                        } from '../modules/nf-core/megahit/main'
include { MAXBIN2                                        } from '../modules/nf-core/maxbin2/main'
include { BUSCO_BUSCO                                    } from '../modules/nf-core/busco/busco/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// MODULE: custom, local modules
//
include { BBMAP_REFORMAT                               } from '../modules/local/bbmap_reformat'
include { FASTP                                        } from '../modules/local/fastp'
include { FASTQSCREEN                                  } from '../modules/local/fastq_screen'
include { UNZIP                                        } from '../modules/local/unzip'
include { UNZIP as UNZIP_POLISHED                      } from '../modules/local/unzip'
include { QUAST                                        } from '../modules/local/quast'
include { BLAST_BLASTN                                 } from '../modules/local/blastn'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

ch_versions = Channel.empty()
ch_multiqc_files = Channel.empty()

workflow METASHORT {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //

    if (!params.skip_subsample) {

        BBMAP_REFORMAT (
            INPUT_CHECK.out.reads,
            params.num_subsamples
        )

        raw_reads = BBMAP_REFORMAT.out.fastq

    } else {

        raw_reads = INPUT_CHECK.out.reads

    }

    FASTQC_RAW (
        raw_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    if (params.trim_tool == "fastp") {

        if (params.use_adapter_reference) {

            FASTP (
                raw_reads,
                ["${params.adapt_ref}"],
                [],
                [],
                []
            )

        } else {

            FASTP (
                raw_reads,
                [],
                [],
                [],
                []
            )
            
        }

        trimmed_reads = FASTP.out.reads
    
    } else if (params.trim_tool == "trimmomatic") {

        TRIMMOMATIC (

            raw_reads

        )

        trimmed_reads = TRIMMOMATIC.out.trimmed_reads

    }

    FASTQC_TRIMMED (

        trimmed_reads

    )

    ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())

    // taxonmic profiling with kraken2
    if (!params.skip_kraken2) {
        
        KRAKEN2_MAIN (

            filtered_reads,
            params.kraken_db_main,
            true,
            true

        )

        KRAKENTOOLS_KREPORT2KRONA (
            KRAKEN2_MAIN.out.report
        )

        KRONA_KTIMPORTTEXT (
            KRAKENTOOLS_KREPORT2KRONA.out.txt
        )


    }

    // contaminant screening and limited taxonmic profiling with fastqscreen
    if (!params.skip_fastqscreen) {

        FASTQSCREEN (
            filtered_reads,
            params.fastq_screen_conf
        )

    }

    if (!params.skip_methaphlan) {

        METAPHLAN (
            filtered_reads,
            params.methaphlan_db
        )

    }

    // filter reads 
    if (!params.skip_extract_kraken_reads) {

        KRAKENTOOLS_EXTRACTKRAKENREADS (

            ["${params.kraken_tax_ids}"],
            KRAKEN2_MAIN.out.classified_reads_assignment,
            KRAKEN2_MAIN.out.classified_reads_fastq,
            KRAKEN2_MAIN.out.report

        )

        filtered_reads = KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads

    } else {

        filtered_reads = trimmed_reads

    }

    if (!params.skip_assembly) {
        
        if (params.assembler == "spades") {
            
            spades_map = filtered_reads.map { meta, fastq -> [ meta, fastq, [], []] }

            contigs_produced = true

            SPADES (

                spades_map,
                [],
                [],         

            )

            assembly_ch = SPADES.out.contigs

        } else if (params.assembler == "flye") {

            FLYE (

                filtered_reads,
                '--nano-hq'

            )

            assembly_ch = FLYE.out.fasta

        } else if (params.assembler == "megahit") {

            MEGAHIT (

                filtered_reads

            )

            assembly_ch = MEGAHIT.out.contigs
        
        }

        UNZIP (

        assembly_ch

        )

        unzip_channel = UNZIP.out.unzip_contigs


        //assembly qc with quast
        QUAST (

            unzip_channel, // consensus (one or more assemblies)

        )

        if (!params.skip_binning) {

            maxbin_map = unzip_channel.join(filtered_reads)

            MAXBIN2 (

                maxbin_map

            )

            BUSCO_BUSCO (

                MAXBIN2.out.binned_fastas,
                "genome",
                "eukaryota_odb10",
                [],
                []

            )


        }

        if (!params.skip_blast) {

            BLAST_BLASTN (
                unzip_channel,
                params.blast_db
            )

        }

    }
    
    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMetashort.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMetashort.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect { it[1] }.ifEmpty([]))
    if (params.trim_tool == "trimmomatic") {

        ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.trim_log.collect{it[1]}.ifEmpty([]))

    } else if (params.trim_tool == "fastp") {

        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] }.ifEmpty([]))

    }
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect { it[1] }.ifEmpty([]))

    if (!params.skip_kraken2) {
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_MAIN.out.report_mqc.collect().ifEmpty([]))
    }
    
    if (!params.skip_assembly) {
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.report.collect().ifEmpty([]))
    }

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
    NfcoreTemplate.dump_parameters(workflow, params)
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
