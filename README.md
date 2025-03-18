## Introduction

**metashort** is a bioinformatics workflow that accepts short reads as input and runs them through the following processes/analyses:


1. _OPTIONAL_: Subsampling ([`BBMap`](https://github.com/BioInfoTools/BBMap))

2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))

3. Adapter Trimming and Quality Filtering ([`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic) or [`fastp`](https://github.com/OpenGene/fastp))

4. Taxonomic Classification ([`kraken2`](https://github.com/DerrickWood/kraken2))

5. Taxonomy distribution visualization ([`Krona`](https://github.com/marbl/Krona))

6. _OPTIONAL_: Taxonomic Filtering ([`KrakenTools`](https://github.com/jenniferlu717/KrakenTools))

7. De Novo Assembly [`Spades`](https://github.com/ablab/spades)

8. Assembly QC ([`quast`](https://github.com/ablab/quast))

9. Binning ([`Maxbin2`](https://sourceforge.net/projects/maxbin2/))

10. Contig alignment and identification ([`blast`](https://www.ncbi.nlm.nih.gov/books/NBK279684/))

11. Generate summary report ([`MultiQC`](http://multiqc.info/))

## Setup

This workflow uses assets and depencies native to the CDC's SciComp environment. If you do not have access to the SciComp environment, you can request an account [`here`](https://info.biotech.cdc.gov/info/helpdesk-ticket/?category=Account%20Requests). 


## Usage

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
SAMPLE_1,/data/reads/sample1-R1.fastq.gz,data/reads/sample1-R2.fastq.gz
SAMPLE_2,/data/reads/sample2-R1.fastq.gz,data/reads/sample2-R2.fastq.gz
```

The top row is the header row ("sample,fastq_1,fastq_2") and should never be altered. Each row below the header, represents two paired-end fastq file with a unique identifier in the "sample" column (SAMPLE_1 and SAMPLE_2 in the example above). Each fastq file needs to be gzipped/compressed to prevent validation errors from occuring at the initialization of the pipeline

There is an example samplesheet located under the assets folder (`assets/samplesheet.csv`) that you can view and edit yourself. **NOTE** If you use this samplesheet, please make a back up copy of it as it will be overwritten each time you pull an updated version of this repository. 

Once the samplesheet has been formatted, we can run the workflow using one of the 3 methods methods listed below.


**Method 1: Cluster Submission**:

The `qsub` method allows you to submit the job to SciComp's high memory cluster computing nodes for fast performance and load distribution. This is a good "fire and forget" method for new users who aren't as familiar with SciComp's compute environment

Format:
```bash
bash ./run_qsub.sh --input "/path/to/samplesheet" --outdir "/path/to/output/directory" "<additional-parameters>"
```

Example:
```bash
bash ./run_qsub.sh --input "assets/samplesheet.csv" --outdir "results/test" "--skip_subsample false --num_subsamples 1000 --skip_kraken2 false"
```


**Method 2: Local Execution**:

The `local` method may be a better option if you are experiencing technical issues with the `qsub` method. `qsub` adds additonal layers of complexity to workflow execution, while `local` simply runs the workflow on your local machine or the host that you're connected to, _**provided it has sufficient memory/RAM and CPUs to execute the workflow**_

Format:
```bash
bash ./run_local.sh --input "/path/to/samplesheet" --outdir "/path/to/output/directory" "<additional-parameters>"
```

Example:
```bash
bash ./run_local.sh --input "./assets/samplesheet.csv" --outdir "./results/test" "--skip_subsample false --num_subsamples 1000 --skip_kraken2 false"
```


**Method 3: Native Nextflow Execution**:

If you are familiar with nextflow and Scicomp's computing environment, you can invoke the `nextflow` command straight from the terminal. **NOTE: if you are using this method you will need to load up a nextflow environment via `module load` or `conda`**
 
Format:
```bash
nextflow run main.nf -profile singularity,local --input "/path/to/samplesheet" --outdir "/path/to/output/directory" \<additional flags\>
```

Example:
```bash
nextflow run main.nf -profile singularity,local --input "./assets/samplesheet.csv" --outdir "./results/test" --skip_subsample false --num_subsamples 1000 --skip_kraken2 false
```


## Parameters

See below for all possible input parameters:


**Global Variables**:
| Parameter | Data Type | Default Value |
|:---------:|:---------:|:-------------:|
| `--metagenomic_sample` | boolean | true |

**Workflow processes**:
| Parameter | Data Type | Default Value |
|:---------:|:---------:|:-------------:|
| `--skip_subsample` | boolean | true |
| `--skip_fastq_screen` | boolean | true |
| `--skip_kraken2` | boolean | true |
| `--skip_extract_kraken_reads` | boolean | true |
| `--skip_metaphlan` | boolean | true |
| `--skip_assembly` | boolean | true |
| `--skip_medaka` | boolean | true |
| `--skip_binning` | boolean | true |
| `--skip_blast` | boolean | true |                       


**BBmap subsampling parameters**:
| Parameter | Data Type | Default Value |
|:---------:|:---------:|:-------------:|
| `--num_subsamples` | integer | 1000 |

**Global Trimming parameters**:
| Parameter | Data Type | Default Value |
|:---------:|:---------:|:-------------:|
| `--trim_tool` | string | "trimmomatic" |
| `--adapt_ref` | string | "./assets/sequencing-adapters.fasta" |

**Trimmomatic parameters**:
| Parameter | Data Type | Default Value | Notes |
|:---------:|:---------:|:-------------:|:-------------:|
| `--trimmomatic params` | string | "ILLUMINACLIP:./assets/sequencing-adapters.fasta:2:30:10 SLIDINGWINDOW:3:20 MINLEN:36" | This string can be modified to any known command line arguments for trimmomatic - simply format the string in the same way you would enter it on the command line |

**fastp parameters**:
| Parameter | Data Type | Default Value | Notes |
|:---------:|:---------:|:-------------:|:-------------:|
| `--adapter_auto_detect` | boolean | false | when set to `false` fastp will use the `--adapter_ref` fasta to locate adapter sequences |
| `--fastp params` | string | "" | This string can be modified to any known command line arguments for fastp - simply format the string in the same way you would enter it on the command line |

**Kraken2 parameters**:
| Parameter | Data Type | Default Value |
|:---------:|:---------:|:-------------:|
| `--kraken_db_main` | string | "/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/Projects/Long_Read_Analysis/data/kraken-db/bact_arch_vir_fungi_amoeba-DB_41-mer" |
| `--kraken_custom_params` | string | "" |

**Kraken tools - Extract Kraken Reads**:
| Parameter | Data Type | Default Value | Notes |
|:---------:|:---------:|:-------------:|:-------------:|
| `--kraken_tax_id` | string | "5754" | Taxonomic ID value(s) that you want `krakentools` to pull out of your classified reads |
| `--include-children` | boolean | true | filter for all child taxonomic IDs of the parent tax ID declared in `--kraken_tax_id` |

**FastqScreen parameters**:
| Parameter | Data Type | Default Value |
|:---------:|:---------:|:-------------:|
| `--fastq_screen_conf` | string | "./assets/fastq_screen.conf" |

**Metaphlan parameters**:
| Parameter | Data Type | Default Value |
|:---------:|:---------:|:-------------:|
| `--methaphlan_db` | string | "/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/Projects/Long_Read_Analysis/data/metaphlan/metaphlan_databases" |

**Assembler paramaters**:
| Parameter | Data Type | Default Value |
|:---------:|:---------:|:-------------:|
| `--assembler` | string | 'spades' |

**BLAST parameters**:
| Parameter | Data Type | Default Value |
|:---------:|:---------:|:-------------:|
| `--blast_db` | string | "/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/Projects/Long_Read_Analysis/data/blast/arch-bact-fung-hum-amoeba_refseq/arch-bact-fung-hum-amoeba_refseq" |
| `--blast_evalue` | string | "1e-10" |
| `--blast_perc_identity` | string | "90" |
| `--blast_target_seqs` | string | "5" |

## Credits

Meta-short was originally written by Sam Rusher (rtq0@cdc.gov)..

## Citations

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
