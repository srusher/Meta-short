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
include { NONPAREIL_NONPAREIL                            } from '../modules/nf-core/nonpareil/nonpareil/main'
include { NONPAREIL_CURVE                                } from '../modules/nf-core/nonpareil/curve/main'
include { NONPAREIL_NONPAREILCURVESR                     } from '../modules/nf-core/nonpareil/nonpareilcurvesr/main'
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
include { UPDATE_NODES_DB                                 } from '../modules/local/update_nodes_db'
include { BBMAP_REFORMAT as BBMAP_REFORMAT_SUBSAMPLE      } from '../modules/local/bbmap_reformat_subsample'
include { FASTP                                           } from '../modules/local/fastp'
include { FASTQSCREEN                                     } from '../modules/local/fastq_screen'
include { KRAKEN2_KRAKEN2 as KRAKEN2_MAIN                 } from '../modules/local/kraken2'
include { MINIMAP2_ALIGN as ALIGN_READS                   } from '../modules/local/minimap2'
include { ALIGNMENT_CLASSIFY                              } from '../modules/local/alignment_classify'
include { SAMTOOLS_VIEW as SAMTOOLS_QUALITY_FILTER        } from '../modules/local/samtools_view'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED         } from '../modules/local/samtools_fastq'
include { BBMAP_REFORMAT as BBMAP_REFORMAT_CLEAN_MAPPED   } from '../modules/local/bbmap_reformat'
include { KRAKEN_ALIGNMENT_COMPARISON                     } from '../modules/local/alignment_and_kraken_comparison'
include { ALIGNMENT_CLASSIFICATION_GRAPH as ALIGNMENT_CLASSIFICATION_GRAPH_READS } from '../modules/local/alignment_classification_graph'
include { UNZIP                                            } from '../modules/local/unzip'
include { UNZIP as UNZIP_POLISHED                          } from '../modules/local/unzip'
include { MINIMAP2_ALIGN as ALIGN_CONTIGS                  } from '../modules/local/minimap2'
include { ALIGNMENT_CLASSIFY_CONTIGS                       } from '../modules/local/alignment_classify_contigs'
include { ALIGNMENT_CLASSIFICATION_GRAPH as ALIGNMENT_CLASSIFICATION_GRAPH_CONTIGS } from '../modules/local/alignment_classification_graph'
include { QUAST                                            } from '../modules/local/quast'
include { BLAST_BLASTN                                     } from '../modules/local/blastn'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// check to make sure only one filtering technique is being used
if (!params.skip_alignment_based_filtering && !params.skip_extract_kraken_reads) {

    println "Error: You can only select one filtering technique - check the parameters [skip_alignment_based_filtering] and [skip_extract_kraken_reads] to ensure you are skipping at least one!"
    System.exit(1)

}

//clearing out kraken and minimap2 queues if memory_saver mode is enabled (only required for local compute; memory allocation should generally be handled by the job scheduler when using the cluster)
if (params.memory_saver) {

    def kraken_queue = new File("${projectDir}/queue/kraken")

    if (kraken_queue.exists() && kraken_queue.isDirectory()) {
        kraken_queue.eachFile { file ->
            file.delete()
        }
    }

    def minimap2_queue = new File("${projectDir}/queue/minimap2")

    if (minimap2_queue.exists() && minimap2_queue.isDirectory()) {
        minimap2_queue.eachFile { file ->
            file.delete()
        }
    }

}


// Info required for completion email and summary
def multiqc_report = []

ch_versions = Channel.empty()
ch_multiqc_files = Channel.empty()

workflow METASHORT {

    if (!params.skip_alignment_based_filtering && params.filter_alignment_by_id) {

        UPDATE_NODES_DB (

            params.local_nodes_db,
            params.ncbi_taxonomy_nodes,
            params.my_tax_ids

        )

        placeholder = UPDATE_NODES_DB.out.complete

    } else {

        placeholder = []

    }

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input),
        placeholder
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //

    if (!params.skip_subsample) {

        BBMAP_REFORMAT_SUBSAMPLE (
            INPUT_CHECK.out.reads,
            params.num_subsamples
        )

        raw_reads = BBMAP_REFORMAT_SUBSAMPLE.out.fastq

    } else {

        raw_reads = INPUT_CHECK.out.reads

    }

    FASTQC_RAW (
        raw_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    if (!params.skip_trimming) {

        if (params.trim_tool == "fastp") {

            if (!params.adapter_auto_detect) {

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

        FASTQC_CH = FASTQC_TRIMMED
        ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())

    } else {

        FASTQC_CH = FASTQC_RAW
        trimmed_reads = raw_reads

    }

    if (!params.skip_nonpareil) {
        
        NONPAREIL_NONPAREIL (

            trimmed_reads,
            "fastq",
            "alignment"

        )

        NONPAREIL_CURVE (

            NONPAREIL_NONPAREIL.out.npo

        )

        // npo_reports = Channel.empty()
        // npo_reports = npo_reports.mix(NONPAREIL_NONPAREIL.out.npo.collect().ifEmpty([]))
        // npo_reports_tuple = npo_reports
        //                     .map {meta, npo -> [[id: 'all'], npo] }
        //                     .groupTuple()

        NONPAREIL_NONPAREILCURVESR (

            NONPAREIL_NONPAREIL.out.npo,
            NONPAREIL_CURVE.out.png

        )

    }

    // taxonmic profiling with kraken2
    if (!params.skip_kraken2) {
        
        KRAKEN2_MAIN (

            trimmed_reads,
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
            trimmed_reads,
            params.fastq_screen_conf
        )

    }

    if (!params.skip_methaphlan) {

        METAPHLAN (
            trimmed_reads,
            params.methaphlan_db
        )

    }

    if (!params.skip_alignment_based_filtering) {

        ALIGN_READS (

            trimmed_reads,
            [[params.minimap2_meta],[params.minimap2_index]],
            true,
            false,
            false        

        )
        
        ALIGNMENT_CLASSIFY (

            ALIGN_READS.out.bam.join(FASTQC_CH.out.html),
            params.seqid2taxid_map,
            params.filter_alignment_by_id,
            params.my_tax_ids,
            params.include_children

        )

        if (params.filter_alignment_by_id) {

            alignment_classified_bam = ALIGNMENT_CLASSIFY.out.classified_plus_filtered_bam

        } else {

            alignment_classified_bam = ALIGNMENT_CLASSIFY.out.classified_bam

        }

        // filtering for reads with a mapping quality at or above params.mapping_quality
        SAMTOOLS_QUALITY_FILTER (

            alignment_classified_bam,
            false

        )

        // capturing aligned reads and converting to fastq
        SAMTOOLS_FASTQ_MAPPED (

            SAMTOOLS_QUALITY_FILTER.out.bam,
            false

        )

        // using bbmap suite to remove empty reads and deduplicate aligned reads
        //bbmap dedup has been removing reads with unqiue names - which is unexpected. So we're going to skip this by default
        if (!params.skip_bbmap_dedup) {

            BBMAP_REFORMAT_CLEAN_MAPPED (

                SAMTOOLS_FASTQ_MAPPED.out.fastq

            )

            // setting filtered_reads channel equal to aligned reads
            filtered_reads = BBMAP_REFORMAT_CLEAN_MAPPED.out.fastq

        } else {

            filtered_reads = SAMTOOLS_FASTQ_MAPPED.out.fastq

        }
    }

    // filter reads 
    if (!params.skip_extract_kraken_reads) {

        KRAKENTOOLS_EXTRACTKRAKENREADS (

            ["${params.kraken_tax_id}"],
            KRAKEN2_MAIN.out.classified_reads_assignment,
            KRAKEN2_MAIN.out.classified_reads_fastq,
            KRAKEN2_MAIN.out.report

        )

        filtered_reads = KRAKENTOOLS_EXTRACTKRAKENREADS.out.classified_reads_fastq

    }

    //if kraken2 and custom alignment are both used for classification then compare the results
    if (!params.skip_kraken2 && !params.skip_alignment_based_filtering && !params.non_standard_reference) {

        KRAKEN_ALIGNMENT_COMPARISON (

            KRAKEN2_MAIN.out.report.join(ALIGNMENT_CLASSIFY.out.summary_tsv)

        )

    } else if (!params.skip_alignment_based_filtering) {

        ALIGNMENT_CLASSIFICATION_GRAPH_READS (

            ALIGNMENT_CLASSIFY.out.summary_tsv

        )

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

        ALIGN_CONTIGS (

            unzip_channel,
            [[params.minimap2_meta],[params.minimap2_index]],
            true,
            false,
            false        

        )
        
        ALIGNMENT_CLASSIFY_CONTIGS (

            ALIGN_CONTIGS.out.bam.join(QUAST.out.qc),
            params.seqid2taxid_map,
            false,
            params.my_tax_ids,
            false

        )

        ALIGNMENT_CLASSIFICATION_GRAPH_CONTIGS (

            ALIGNMENT_CLASSIFY_CONTIGS.out.summary_tsv

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

    if (!params.skip_trimming) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect { it[1] }.ifEmpty([]))
        
        if (params.trim_tool == "trimmomatic") {

            ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.trim_log.collect{it[1]}.ifEmpty([]))

        } else if (params.trim_tool == "fastp") {

            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] }.ifEmpty([]))

        }
    }

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
