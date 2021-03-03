# SARS-CoV-2 training with Nextflow

In this report you will find all the information necessary to follow the steps to analyze SARS-CoV-2 data with Galaxy.

## Training overview
During this training we will following these steps:
* [Conda](#conda): Installation
* [Nextflow](#nexflow): Overview and installation
* [Viralrecon](#viralrecon): Overview and download
	* [Github](#github): Glance to main commands
* [Pipeline](#pipeline): Running the pipeline with command line and main results.
	* [Results](#results): Most important results.
* [Statistics](#statistics): Parse resulting files to obtain statistics.
* [Lineage](#Lineage): Lineage clasification with Pangolin

## Conda
In case you haven't install conda yet, follow the intructions in [this link](https://docs.conda.io/en/latest/miniconda.html)

## Nextflow
[Nextflow](https://www.nextflow.io/) is a bioinformatics workflow manager that enables the development of portable and reproducible workflows. It supports deploying workflows on a variety of execution platforms including local, HPC schedulers, AWS Batch, Google Cloud Life Sciences, and Kubernetes. Additionally, it provides support for manage your workflow dependencies through built-in support for Conda, Docker, Singularity, and Modules.

### Installation
We will install Nextflow using conda running the following commands:

Create an environment
```
conda create --name nextflow
```

proceed ([y]/n)?
```
y
```

Activate the environment:
```
conda activate nextflow
```

Install nextflow in the environment through bioconda channel:
```
conda install -c bioconda nextflow
```

proceed ([y]/n)?
```
y
```

Test if it's correctly installed:
```
nextflow -v
```

## Viralrecon

[viralrecon](https://github.com/nf-core/viralrecon) is a bioinformatics analysis pipeline used to perform assembly and intra-host/low-frequency variant calling for viral samples. The pipeline supports short-read Illumina sequencing data from both shotgun (e.g. sequencing directly from clinical samples) and enrichment-based library preparation methods (e.g. amplicon-based: ARTIC SARS-CoV-2 enrichment protocol; or probe-capture-based).

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with Docker containers making installation trivial and results highly reproducible. Furthermore, automated continuous integration tests that run the pipeline on a full-sized data set using AWS cloud ensure that the code is stable.


### GitHub
To download and run viralrecon locally first you have to learn a little bit about github in command line (git). If you need more information about how to manage GitHub repositories, read our [Wiki documentation](https://github.com/BU-ISCIII/BU-ISCIII/wiki/Github--gitflow). You need to follow these steps:

In case you don't have git installed run:
```
sudo apt install git-all
```

Download the repo from github:
:warning: Copy the viralrecon link below in the code, the one above won't work.
```
git clone git@github.com:BU-ISCIII/viralrecon.git
```

Moving inside the repo:
```
cd viralrecon
```

Moving to another branch:
```
git checkout dev
```

Leave repo's folder:
```
cd ..
```

Now you have your repo (viralrecon) installed and the program to run it (Nextflow), the next step is to run the pipeline.

## Pipeline
To run the pipeline you have to type in your terminal:

```
nextflow run Repositories/viralrecon/main.nf -profile conda,test -resume
```
With this command line you are going to run the default test data set with all the default parameters. This data is a subset of amplicon based sequenced SARS-CoV2 samples. This means that during the pipeline the sequences/positions corresponding to the amplicon's primers are going to be removed.

This is the pipeline overview:
1. FastQC: Quality control
2. Fastp: Quality+size trimming
3. Mapping approach:
	1. bowtie2: Mapping to the reference genome
	2. ivar: Amplicon's adapter trimming by position
	3. VarScan2: Variant calling
	4. bcftools: Consensus generation
4. Assembly approach:
	1. cutadapt: Adapter trimming by sequence
	2. MetaSpades: De Novo assebmly
	3. ABACAS: Reference-based scaffold ordering

### Results
Now we are going to see the most important results of the pipeline. In case you were not able to run the pipeline, you can see the [results in this folder](../results/).

One of the most interesting results is the [MultiQC report](../results/multiqc/multiqc_report.html). In this report you have an overview of all the steps of the pipeline and their results.

The first results consist in a summary of both approaches performed by viralrecon (mapping and assembly). It gives different statistical measures such as *the number of input reads*, *the number of reads remaining after trimming*, *the percentage of mapped reads* and more.

![multiqc](../docs/images/multiqc_1.png)

You can interactively move across the report to see more results, such as the number of variants called by VarScan2 for each sample in a barplot.

![varscan2](../docs/images/varscan2.png)


Also, another interesting result is the one showed by [PlasmidID](https://github.com/BU-ISCIII/plasmidID) in [a circos plot](../results/assembly/metaspades/plasmidid/SAMPLE1_PE/images/SAMPLE1_PE_NC_045512.2.png). In the *de novo* asssembly approach, the resulting scaffolds assemblies are plotted against the reference SARS-CoV2 genome from Wuhan. In this case we can see that the whole reference genome is covered and represented by the contigs 1 and 2 from the assembly fasta.

![plasmid](../results/assembly/metaspades/plasmidid/SAMPLE1_PE/images/SAMPLE1_PE_NC_045512.2.png)
