# atac-seq-pipeline
Pipeline to analyze ATACseq data
##################################################################################################<br>
# Developed by Tsviya Olender, Molecular Genetics Department, Weizmann Institute #<br>
last update: 26 June 2018<br>
The script is still under development and might contains bugs.<br>
For any problem: please write to tsviya.olender@weizmann.ac.il, I will appreciate your feedback # <br>
##################################################################################################<br>
<br>
######################################################################################################<br>
# INSTALLATION<br>
#####################################################################################################<br>
The script has a hard coded parameter =>$progD<br>
which defines the location of the code.<br>
Should be changed when the package is moved.<br>
The script collect_bedtools_count.pl has a hard coded parameter=> $pipelinePATH<br>
Should contain the location of the pipeline<br>
#####################################################################################################<br>
LOCATION OF FASTQ<br>
Very often, crude reads come in multiple fastq files per sample. In this case, the script expects a file structure with a folder per sample,<br>
all samples in the same folder. The fastq files should be gzipped, with R1 for R1 and R2 for R2.<br>
In this case, the parameter file has to be defined as follows:<br>
1. Define the crude reads location: by the parameter crude_reads_location_for_merging<br>
2. set the option combine_fastq to 1<br>
<br>
Otherwise, if you have allready 2 fastq files per sample (one for R1, and one for R2)- put the files under 1_fastq, in the working directory.<br>
Set the option combine_fastq to 0<br>
The fastq files much be named as: XXX_R1.fastq.gz, XXX_R2.fastq.gz- where XXX is the sample name.<br>
(that means gziped)<br>
########################################################################################################<br>
# USAGE<br>
########################################################################################################<br>
# the script atacpipeline.pl send the queries to the server.<br>
It runs the script run_ATAC_Ts_V2.pl in parallel<br>

The script accepts the following:<br>
1. name of file with names of all samples. the file should contain a row per sample. <br>
2. parameters file, e.g. run_ATAC_Ts_V2_params.txt<br>
3. name of the LSF queue. The default is new-short.<br>
4. The memory (default is 4000). In order to change the default memory, one has to define the queue in the command line, otherwise the script will not work. <br>
Note that the script runs in 4 threds, therefore a memory of 4000G means 16000Gb memory per job.<br>
perl PATH_TO_prog/atacpipeline.pl samples.txt params_file [queue_name memory]<br>
example<br>
perl PATH_TO_prog/atacpipeline.pl samples.txt run_ATAC_Ts_V2_params.txt<br>
submits all samples listed in samples.txt to the analysis. The parameters are defined in the file run_ATAC_Ts_V1_params.txt.<br>
The jobs will be submitted to the queue new-all with memory of 4000G.<br>

perl PATH_TO_prog/atacpipeline.pl samples.txt run_ATAC_Ts_V2_params.txt new-short 8000<br>
here the jobs was submitted with more memory. Because we changed the default memory we add also to define the queue.<br>

perl PATH_TO_prog/atacpipeline.pl samples.txt run_ATAC_Ts_V2_params.txt new-all.q <br>
here we submitted the job to different queue, with the default memory<br>

perl PATH_TO_prog/atacpipeline.pl samples.txt run_ATAC_Ts_V2_params.txt new-all.q 8000<br>
both: queue name and memory are not the default options<br>
########################################################################################################<br>
# Generate a quality report<br>
########################################################################################################<br>

The script collect_qual_params.pl generates a report with the quality parameters of the run (e.g. num of reads in every analysis step).<br>
# run example:<br>
perl PATH_TO_prog/collect_qual_params.pl samples.txt > qual_report.txt<br>
This script also generates a file called TSS_counts_table.txt with the crude read counts, and the counts in FPKM of the promoters.<br>
########################################################################################################<br>
# Explanation on the parameter file<br>
########################################################################################################<br>
this files allows the user to set up the run parameters.<br>
It is build from 2 parts:<br>
[params] and [setup_run]<br>
the parameters in params are self explanatory<br>
in the set up the user decides which parts of the analysis will be processes<br>
0 mean no, 1 means yes<br>
collectReads = 0 => In case the fastq is divided into several part. The script excpet to find a folder per sample, gzipped.<br>
R1 files should be named as *R1* (meaning- they should contain the string R1 in the file name), and R2 for R2 files.<br>
The location of the crude reads is defined by crude_reads_location_for_merging<br>
fastqc = 0  => run fastqc<br>
trim_adapter = 1 => runs cutadapt<br>
make_body = 1 => runs bowtie alignment, and process the alignment to contain only uniquily mapped reads. the opuput is sorted bam<br>
make_plots = 1 => generates ngsplot plots<br>
nucleosome_free = 1 => filters the bam file, to contain only paired-end reads with < 130bp insert<br>
call_peaks = 1 => uses MACS2 to call peaks<br>
countTSS_reads = 1 => counts the reads on the TSS reagions, performs FPKM normalization. The report contains crude read count, and FPKM count<br>



