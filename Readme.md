# Sentieon BWA-MEM alignment for UMI

## Please read this important information before running the app

DNAnexus has partnered with [Sentieon](http://www.sentieon.com) to bring the Sentieon suite of secondary analysis tools to the DNAnexus platform. In order to run this app, you will need to request access by contacing info@sentieon.com and providing your DNAnexus organization or username and the region(s) you plan on running the app in. If you encounter any issue with this app, please contact support@sentieon.com.

**Software License Agreement**

Please review terms and conditions for using Sentieon software on DNAnexus: https://www.sentieon.com/EULA/eula-dnanexus.html

**Free trial**

After running one of the Sentieon apps for the first time, you will receive a 2 week trial evaluation license for all the Sentieon apps. If you would like to continue using the Sentieon tools after the trial period is completed, please contact info@sentieon.com.

## What does this app do?

This app implements [Sentieon](http://www.sentieon.com)'s pre-processing and alignment pipeline for next-generation sequence data while taking advantage of molecular barcode information (also called unique molecular indices or UMIs that introduce unique tags on the ends of template DNA molecules prior to sequencing to greatly reduce the impact of PCR duplicates and sequencing errors on the variant calling process). It receives raw FASTQ reads as input, extracts the UMI tags from read data and performs barcode-aware consensus generation, and finally maps the consensus reads to a reference genome using an optimized BWA. The app outputs the sorted mappings in BAM or CRAM format.

## What data are required for this app to run?

This app requires:

- An array of either single gzipped FASTQ files (for single-end experiments) or paired FASTQ files (for paired-end experiments).
- The reference genome and a BWA indexed reference genome.
- A read structure that defines the format that the sequence of the read follows.
- An optional read group info file in CSV format can be provided to specify the read group information for the input FASTQ files. 
Each row of the CSV file corresponds to one pair of FASTQ file, and contains four read group attributes separated by comma: 
  1. filename of the first FASTQ mate, 
  2. RG ID (read group identifier) tag, 
  3. RG LB (read group library) tag, 
  4. RG PU (read group platform unit) tag. 

  Note that the filenames in the first column need to include all the intput FASTQ filenames of the first mate. 

  Below is an example:
  ```
  SRR123453_1.fastq.gz,SRR123453_1,Library1,SRR123453_1 
  SRR123454_1.fastq.gz,SRR123454_2,Library1,SRR123454_2 
  SRR123455_1.fastq.gz,SRR123455_3,Library2,SRR123455_3 
  SRR123456_1.fastq.gz,SRR123456_4,Library3,SRR123456_4 
  SRR123457_1.fastq.gz,SRR123457_5,Library3,SRR123457_5 
  SRR123458_1.fastq.gz,SRR123458_6,Library4,SRR123458_6
  ```

## What does this app output?

This app outputs the following files:

- The sorted and mapped mappings file in BAM (`*.bam`) or CRAM (`*.cram`) format along with the associated BAM index file (`*.bai` or `*.crai`).
- Optionally, a collection of (`*_metrics/`) files containing the statistics metrics of the reads (`*.gc_summary.txt`, `*.gc_metric.txt`, `*.mq_metric.txt`, `*.qd_metric.txt`, `*.is_metric.txt`, `*.aln_metric.txt`, and `*.metrics.pdf`).

## How does this app work?

This app performs the following steps:

- Fetches your gzipped FASTQ input file(s) and the reference genome data.
- Runs Sentieon umi extract to extract the barcode sequences from the input reads and tags the reads with the corresponding barcode.
- Runs BWA to map the reads to the reference.
- Runs Sentieon umi consensus to create consensus molecules from aligned barcode-tagged reads.
- Runs BWA to map the consensus reads to the reference and sorts the output.

