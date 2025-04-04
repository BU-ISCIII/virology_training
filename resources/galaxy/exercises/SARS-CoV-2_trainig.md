# SARS-CoV-2 training with galaxy

In this report you will find all the information necessary to follow the steps to analyze SARS-CoV-2 data with Galaxy.

## Training 

During this training we will following these steps:

* [Register/Login](#registerlogin): Register or login into Galaxy website
* [Data](#data): Upload data for the analysis
* [Quality](#quality): Analysis of the quality of the raw reads
* [Trimming](#trimming): Quality trimming using fastp
* [Mapping](#mapping): Mapping reads to reference genome with Bowtie2
* [Stats](#stats): Mapping statistics with picard
* [Amplicons](#amplicons): Preprocessing steps mandatory for amplicon sequencing data
* [Variants](#variants): Variant calling and annotation
* [Consensus](#consensus): Consensus genome generation
* [Additional info](#additional-info): Websites for additional info.

## Register/Login

First of all, you have to create a user (or login into your user if you already have one) in the [European Galaxy website](https://usegalaxy.eu/) by clicking in the button indicated in the image:

![galaxy_login](../docs/images/galaxy_login.png)

You can do this analysis without creating a user, but your history will be saved into your user so you can access later at home.

Once you are registered we can start with the analysis workflow.

**_From now on, each job we run in Galaxy will have a unique number for identifying each process. This numbers can differ depending on the number of samples and the times you run or delete any process. This training's snapshots were taken using two samples and some process were deleted for any reason, so numbers MAY DIFFER. Also, this exercise's snapshots were taken with other dataset, so sample names may differ too._**

## Data

Before starting with any analysis we have to upload de data we want to analyze into galaxy. The are three different ways to upload data to Galaxy, from which we will explain you some of them along this training.

### Uploading data from URL

In order to do upload files from URL we have to follow these steps:

1. In the left side panel, select **Upload**
2. In the new panel select **Paste/Fetch Data**
3. Then copy the following block of text for the samples:

```bash
https://zenodo.org/record/5724464/files/SARSCOV2-1_R1.fastq.gz?download=1
https://zenodo.org/record/5724464/files/SARSCOV2-1_R2.fastq.gz?download=1
```

4. Now, in the **Download data from the web by entering URLs (one per line) or directly paste content.** square, paste the text you copied before
5. Select **Start**
6. When everything is green in the screen, select *Close*

<img src="../docs/images/Upload_1.png" alt="Upload 1" width="700"/>
<img src="../docs/images/Upload_2.png" alt="Upload 2" width="700"/>

Once files are uploaded we can close the window clicking on "_close_" button. Both files will appear in the right panel.

If we select the "eye" in the datasets from the history, we can see the content of the data:

<img src="../docs/images/fastq.png" alt="fastq" width="200"/>

## Quality

### Quality Analysis (FastQC)

Once we have the raw data, an important step is to analyze the quality of the reads, to know if the reads are worth it. To do this, we have to:

1. Search for the **fastqc** tool
2. Select **FastQC Read Quality reports** and set the following parameters:
3. Select _Multiple datasets_ in _Raw read data from your current history_
4. then select in the bar the two datasets from the history
5. Then go down and select **Run tool**

<p align="center"><img src="../docs/images/fastqc_run.png" alt="fastqc_run" width="900"></p>

This program will generate a message like this one, were we can read that, each .fastq file is going to generate two different jobs, one for the Raw Data and another one for the .HTML report.

![fastqc_message](../docs/images/fastqc_message.png)

#### FastQC results visualization

To visualize the information coming from FastQC we just have to select the job of interest. In this case we are interested in the "_Web page results_" so for the sample we want to see the results we have to click in the _eye_ to visualize galaxy results:

![fastqc_results](../docs/images/fastqc_results.png)

This report gives different type of information about the sequencing quality of the reads. Clicking in the arrows in the bottom right/left corners of the page we can hide/show the side panels.

![bottom_arrows_1](../docs/images/bottom_arrows_1.png)

So the central panel with the results we want to visualize will bee better seen.

![bottom_arrows_2](../docs/images/bottom_arrows_2.png)

**_For more information about FastQC output visit [FastQC website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)_**

## Trimming

### Quality trimming (Fastp)

Once we have check the quality of our reads, it's important to trim low quality nucleotides from those reads, for which we will use _Fastp_. So, in the search bar you look for fastp and then select "_fastp - fast all-in-one preprocessing for FASTQ files_". There, we will have to change some parameters ensure the trimming accuracy for this amplicon data. First of all we are going to do the analysis for the sample we gave to you (201569). These are the field we will have to change:

3.Single-end or paired reads > Paired
4.Input 1 > Select `SARSCOV2-1_R1.fastq.gz?download=1`

   Input 2 > Select `SARSCOV2-1_R2.fastq.gz?download=1`
  
5.Display Filter Options
   - Quality Filtering options
     
     6.Qualified Quality Phred = 20

     7.Unqualified percent limit = 10

   - Length Filtering Options

     8.Length required = 50

9.Read modification options

  10.PoliX tail trimming > Enable polyX tail trimming

- Per read cutting by quality options
  
  11.Cut by quality in front (5') > Yes
  
  12.Cut by quality in tail (3') > Yes

13."_Execute_"

![fastp1](../docs/images/fastp1.png)
![fastp2](../docs/images/fastp2.png)
![fastp3](../docs/images/fastp3.png)

A message like this one will appear, which means that 4 results will be generated:

  1. One with the R1 trimmed reads: `fastp on data 2 and data 1: Read 1 output`
  2. Another one with the R2 trimmed reads: `fastp on data 2 and data 1: Read 2 output`
  3. Another one with the HTML results: `fastp on data 2 and data 1: HTML report`
  4. A last one which is not relevant for this training

![fastp_message](../docs/images/fastp_message.png)

#### Fastp results

Once fastp analysis is done, you can see the results by clicking in the eye ("_View Data_") in the fatp HTML results. We will see a report like this one:

![fastp_results](../docs/images/fastp_results.png)

This report can be downloaded by clicking in the results name ("_fastp on data 2 and data 1: HTML report_") and then in "_Download_".

![fastp_download](../docs/images/fastp_download.png)

This will download a .zip folder that contains the .html report, that can be visualize in your computer using any browser (such as Google Chrome) and interactively scroll around the results.

Among the most relevant results, you have the:

- Summary: Stats summary
  - After filtering: Statistics of the reads after quality filtering
    - reads passed filters: Reads remaining after quality filter trimming
    - reads with low quality: Reads that were remove due to low quality
    - reads too short: Reads that didn't pass the minimum length filter.
  - After filtering: Plots after filtering
    - After filtering: read1: quality: Plot with the evolution of R1 quality over read position. Usually it decays in the last nucleotides.
    - After filtering: read2: quality: Same plot for R2.

**_For more information about FastQC output visit [Fastp github](https://github.com/OpenGene/fastp)_**

*:warning: If you are working with two samples, don't forget to run [trimming steps](#quality-trimming-fastp) with the other sample.*

## Mapping

In order to call for variants between the samples and the reference, it's mandatory to map the sample reads to the reference genome. To do this we need the fasta file of the reference and the Bowtie2 index of that fasta file.

### Mapping reads with reference genome (Bowtie2)

Now we can start with the main mapping process. The first thing we have to do is look for the program "_Bowtie2_" in the search bar and then select "_Bowtie2 - map reads against reference genome_". Here we will have to set the following parameters, for the first sample:

3.Is this single or paired library > Paired-end

4.`Fasta/Q file #1`: **fastp on data 2 and data 1: Read 1 output**

5.`Fasta/Q file #2`: **fastp on data 2 and data 1: Read 2 output**

-Will you select a reference genome from your history or use a built-in index?

    6.Select reference genome > `SARS-CoV-2 isolate Wuhan-Hu-1, complete genome (NC_045512.2)`

7.Do you want to use presets? > Very sensitive local

8.Save the bowtie2 mapping statistics to the history > Yes

9.Run tool

![bowtie1](../docs/images/bowtie1.png)
![bowtie2](../docs/images/bowtie2.png)

We will see a message like this one:

![bowtie_message](../docs/images/bowtie_message.png)

#### Mapping results

Now we can see the mapping results for the samples. The bowtie2 resulting file is a .bam file, which is not easy to read by humans. This .bam file can be downloaded by clicking in the alignment file and then into download. Then, the resulting .gz file will contain the alignment .bam file that can be introduced in a software such as [IGV](http://software.broadinstitute.org/software/igv/) with the reference genome fasta file.

![bowtie2_bam](../docs/images/bowtie2_bam.png)

In our case, the file that can be visualize is the statistics file, which contains information such as the percentage of reads that aligned.

![bowtie2_results](../docs/images/bowtie2_results.png)

_:warning: If you are working with two samples, don't forget to run [trimming steps](#quality-trimming-fastp) with the other sample._

## Stats

The previously shown files give few human readable information, because mapping files are supposed to be used by other programs. In this sense, we can use some programs to extract relevant statistical information about the mapping process.

### Picard CollectWgsMetrics

Another program that gives statistical information about the mapping process is Picard. To run this program you just have to search "_Collect Wgs Metrics_" and then select "_CollectWgsMetrics compute metrics for evaluating of whole genome sequencing experiments_".

3. Using reference genome > `SARS-CoV-2 isolate Wuhan-Hu-1, complete genome (NC_045512.2)`
4. Run tool

![picard_wgsmetrics1](../docs/images/picard_wgsmetrics1.png)

This process will generate one output file per .bam alignment file selected as input.

![picard_wgsmetrics_message](../docs/images/picard_wgsmetrics_message.png)

#### Picard results

Picard results consist in quite long files, so the best is to download those results and visualize them in your computer. Yo you have to click in the CollectWgsMetrics job you want to download, and then click in the save button:

![download_picard](../docs/images/download_picard.png)

Then you just have to open the file with Excell in your computer, and you will see a file with different columns with information about the percentage of the reference genome that is covered by the reads at a specific depth or the mean depth of coverage of the reference genome.

![picard_results](../docs/images/picard_results.png)

So in the results table you can see that the "Mean Coverage" is around 115, which means that each position in the reference genome is supported by 115 reads, by mean. The "PCT 10X" represents the percentage of the reference genome that is sequence with a 10X depth, which represents the percentage of the reference genome that is at least sequenced 10 times (79,57%).

## Amplicons

After mapping the reads to the reference genome, we are interested in removing the sequences of the amplicon primers. To do that you will use a program called iVar, and you will need a bed file with the positions of those amplicon primers.

### Trim amplicon sequences

You have to search for "_ivar trim_" in the search bar and select "_ivar trim Trim reads in aligned BAM_". Then follow these steps:

3. Source of primer information > `Built-in`
4. Primer scheme name > `SARS-CoV-2 ARTICv3`
5. Include reads not ending in any primer binding sites? > `Yes.`
6. Require a minimum length for reads to retain them after any trimming? > `Yes, and provide a custom threshold`
7. Minimum trimmed length threshold = `20`
8. Run tool

![ivar_trim1](../docs/images/ivar_trim1.png)

![ivar_trim2](../docs/images/ivar_trim2.png)

![ivar_trim3](../docs/images/ivartrim_message.png)

#### iVar trimming results

The resulting file from iVar will be a new BAM file where amplicon primer positions will be removed, so there's no result to visualize.

## Variants

Once we have the alingment statistics and files with amplicon primers trimmed, we can start with the variant calling process.

### iVar

#### Reference genome for calling variants

We have to upload the fasta reference genome, and we are going to do it using the following URL:

```bash
https://zenodo.org/record/5724970/files/GCF_009858895.2_ASM985889v3_genomic.200409.fna.gz?download=1
```

So you have to copy the URL direction, and do then select "_Download from web or upload from disk_", the select "_Paste/Fetch data_" and in the new window that appears you should paste the URL above and finally select "_Start_". Now the fasta file of the reference is download.

![reference_genome_download](../docs/images/reference_genome_download.png)

To call for variants between the sample and the reference we are going to use iVar variants, so you have to search for "_ivar_" in the search bar and select "_ivar variants Call variants from aligned BAM file_". Then select the following parameters:

3. In "Bam file" make sure it is selected ivar Trimmed bam files
4. Reference > `GCF_009858895.2_ASM985889v3_genomic.200409.fna.gz?download=1`
5. Minimum frequency threshold > `0.25`
6. Output format > `Both Tabular and VCF`
7. Run tool

![ivar_variants1](../docs/images/ivar_variants1.png)

iVar will create two output files per each bam file uploaded, one with the tabular results, and another one with the vcf files.

![ivar_variants_message](../docs/images/ivar_variants_message.png)

#### iVar results

To display iVar variants results, select the :eye: icon in the right pannel of the Tabular results and then extend the central panel by clicking in the arrows in the bottom.

![ivar_variants_results1](../docs/images/ivar_variants_results1.png)

iVar results consist in a Tab separated file containing all the variants found between the reference and the sample with an Alle Frequency higher than threshold (0.25). Each line represents a variant and the columns give information about that variant, such as the position in the reference genome, the reference allele, the alternate allele, if that variant passed the filters, and so on.

This variants have passed a filter for the minimum quality of the variant, which we set as 20, and the allele frequency threshold. However, we will have a closer look only to those variants present in an allele frequency higher than 0.75 which are the ones that are going to be included in the consensus.

### Annotation with SnpEff

Once we have the variants called, it's interesting to annotate those variants, for which you will use SnpEff. Search for "_snpeff_" in the searh bar and select "_SnpEff eff: annotate variants for SARS-CoV-2_", then change the following parameters:

3. In "Sequence changes (SNPs, MNPs, InDels)" make sure it is selected iVar variants results
4. Run tool

![snpeff](../docs/images/snpeff.png)
![snpeff_message](../docs/images/snpeff_message.png)

### SnpEff results

The SnpEff gives three different results, from which the most interesting ones are:

1. Snpeff eff: Which is a VCF file with the annotation results. It is a very similar file to the ones we saw before for VarScan and Bcftools but with the last column different, containing relevant information about that variant.

![snpeff_results1](../docs/images/snpeff_results1.png)

2. Snpeff eff HTML stats: This file is a CSV file that contains statistics about the variant annotation process, such as the percentage of variants annotated, the percentage of variants that are MISSENSE or SILENT, the percentage that have HIGH, MODERATE or LOW effect, and so on.

![snpeff_results2](../docs/images/snpeff_results2.png)

## Consensus

Once we have the most relevant variants that can be considered to include in the consensus genome, you can start with the consensus genome generation.

### iVar consensus

Now we are going to generate the consensus genome using iVar. We are going to search for _ivar_ and the select "_ivar consensus Call consensus from aligned BAM file_". Now follow these steps:

3. In "BAM fle" make sure iVar trim bam is selected
4. Minimum frequency threshold > `0.75`
5. Run tool

![ivar_consensus](../docs/images/ivar_consensus.png)

![ivar_message](../docs/images/ivar_message.png)

This will just generate a fasta file identical to the reference one, except for those nucleotides that are variants from the VCF file.

![ivar_consensus_results](../docs/images/ivar_consensus_results.png)

## Lineage

Now we are going to determine the lineage of the samples. We will use a software called pangolin. We are going to use the masked consensus genomes generated in the previous steps as follows:

3. Search for the **pangolin** tool
4. Select **Pangolin Phylogenetic Assignment of Outbreak Lineages** and set the following parameters:
5. In "Input FASTA File(s)" make sure iVar consensus fasta is selected
6. Version of pangolin-data to use > `Download latest pangolin-data version from web`
7. Version of constellations to use > `Download latest pangolin-data version from web`
8. **Run tool**

<p align="center"><img src="../docs/images/pangolin.png" alt="pangolin" width="900"></p>

<p align="center"><img src="../docs/images/pangolin_message.png" alt="pangolin_messasge" width="500"></p>

Now we are going to have a look to the results from pangolin:

<p align="center"><img src="../docs/images/pangolin_results.png" alt="pangolin_results" width="900"></p>

As you can see, results are in table format, where you have in first place the reference genome, then de lineage assingned, which seems to be [B.1.1.7 lineage](https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html). This lineage corresponds as you can see in the 5th column to the Alpha lineage.

## All results

If you have any problem following this training and you want to visualize the resulting file you can access them through [this URL](https://usegalaxy.eu/u/svarona/h/sars-cov-2-1)

And viralrecon workflfow in [this URL](https://usegalaxy.eu/u/svarona/w/viralreconupdated)

## Additional info

* [GiSaid](https://gisaid.org/): The GISAID Initiative promotes the rapid sharing of data from all influenza viruses and the coronavirus causing COVID-19. This includes genetic sequence and related clinical and epidemiological data associated with human viruses, and geographical as well as species-specific data associated with avian and other animal viruses, to help researchers understand how viruses evolve and spread during epidemics and pandemics.
* [ENA](https://www.ebi.ac.uk/ena/browser/): The European Nucleotide Archive (ENA) provides a comprehensive record of the world’s nucleotide sequencing information, covering raw sequencing data, sequence assembly information and functional annotation.
* [Nextstrain](https://nextstrain.org/ncov/gisaid/global/6m): Nextstrain is a project to harness the scientific and public health potential of pathogen genome data. Our goal is to aid epidemiological understanding of pathogen spread and evolution and improve outbreak response.
* [outbreak.info](https://outbreak.info/): Outbreak.info aggregates data across scientific sources, providing tools to meet three major aims: (1) Track daily developments in SARS-CoV-2 variants. (2) Integrate publications, preprints, clinical trials, datasets, protocols, and other resources into one searchable library of COVID-19 research. (3) Track trends in COVID-19 cases and deaths.