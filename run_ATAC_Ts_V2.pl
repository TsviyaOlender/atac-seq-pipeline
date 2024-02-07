#!/usr/local/bin/perl
##############################
use Config::Simple;
use File::Path qw(make_path remove_tree);
##############################################################
# NOTE: the script requires loading samtools and bedtools
##############################################################

my($sample,$paramFile) = (@ARGV);

## load modules
load_modules();

## read run params
read_params($paramFile);

###################################################################
## folders definition
my $fastqDir = "1_fastq";
my $fastqcD = "2_fastqc";
my $outD="3_align";
my $plotsD = "4_plots";
my $nucleosomeFreeD = "5_nucleosome_free";
my $MACS_dir = "6_MACS_2";
my $tssD = "7_TSS";
###################################################################
# files for step 0
#fastq files before gzip
#$fastq1_1 = "$fastqDir/"."$sample"."_R1.fastq";
#$fastq2_1 = "$fastqDir/"."$sample"."_R2.fastq";
$fastq1_1 = "fastq/"."$sample"."_R1_001.fastq";
$fastq2_1 = "fastq/"."$sample"."_R2_001.fastq";

#fastq files after gzip
#$fastq1 = "$fastqDir/"."$sample"."_R1.fastq.gz";
#$fastq2 = "$fastqDir/"."$sample"."_R2.fastq.gz";
$fastq1 = "$fastqDir/"."$sample"."_R1_001.fastq.gz";
$fastq2 = "$fastqDir/"."$sample"."_R2_001.fastq.gz";
###################################################################
# files for step 1
#fastq after cutadapt
$fastq1_t = "$fastqDir/"."$sample"."_R1_t.fastq.gz";
$fastq2_t = "$fastqDir/"."$sample"."_R2_t.fastq.gz";
###################################################################
# files for step 2
#initial sam file
my $hits = "$sample".".sam";
# bam file, sorted bam file, bam after remove PCR duplicates
my $mapped = $sample."_mapped\.bam";
my $sorted = $sample."_sorted\.bam";
my $rem_dup = $sample."_rem\.bam";
# for PICARD
my $metrics = $sample."_metrics\.txt";
# flagstat report
my $flagstat=$sample.".flagstat";
# alignment report
my $log ="$outD/$sample".".log";
###################################################################
# files for step 3: plots
my $plotF = $sample."\.pdf";
my $plotTSS = $sample."\_tss.pdf";
my $plotGEN = $sample."\_genbody.pdf";
###################################################################
# files for step 4: nucleosome free
$tempB = "$nucleosomeFreeD/$sample"."_free.bam";
$tempF = "$nucleosomeFreeD/$sample"."_temp.sam";
$tdf = "$nucleosomeFreeD/$sample".".tdf";
$bedF = "$nucleosomeFreeD/$sample".".bed";
$bedFs = "$nucleosomeFreeD/$sample".".s.bed";
$flagstatR = "$nucleosomeFreeD/$sample"."_free.flagstat";
###################################################################
# files for step 5: count of reads in TSS
$TSScountF = "$tssD/"."$sample"."_TSScounts.txt";
###################################################################
 

# step 0- collect fastq
if($collectReads == 1){
  collect_reads($crude_reads_location);
}

#step 0_1- run fastqc
   #run fastqc
if($fastqc == 1){
  run_fastqc($fastq1,"$sample","1",$fastqcD);
  run_fastqc($fastq2,"$sample","2",$fastqcD);
}


# step 1: trim
if($trim_adapter == 1){
	system("cutadapt -q 25 -a $adaptor1 -A $adaptor2 --minimum-length 30 -o $fastq1_t -p $fastq2_t $fastq1 $fastq2"); 
}
# step 2: bowtie alignment,
# filter alignment, mark duplicates, index bam file, generate flagstat report, remove temporary files
if($make_body ==1){ 
	print "bowtie2 -X2000 --local -p4 --mm -x $genome -1 $fastq1_t -2 $fastq2_t -S $outD/$hits >& $log\n";
	system("bowtie2 -X2000 --local -p4 --mm -x $genome -1 $fastq1_t -2 $fastq2_t -S $outD/$hits >& $log");
	process_alignment();
}

#step 3: make plots
if($make_plots==1){
  make_plots();

}

# step 4 - generate a bam file with nucleosome free regions
if($nucleosome_free == 1){
  make_nucleosome_free("$outD/$rem_dup");
}
# step 5 - peak calling
if($call_peaks == 1){
  system("macs2 callpeak -t $tempB --bw 120 -B -f BAMPE --SPMR -B -g mm --nomodel --shift -50  --extsize 100 --broad -n $sample --keep-dup all --outdir $MACS_dir");
}
# step 6 - do read count per TSS
# to overcome memory issues: I generate bed file
# sort the file
# and use the option
# -sorted- to keep the memory usage low
if($countTSS_reads == 1){
  system("bedtools bamtobed -i $tempB >$bedF");
  system("sort -k1,1 -k2,2n $bedF > $bedFs");
  system("bedtools coverage -counts -sorted -a $TSSfile -b $bedFs > $TSScountF");
  system("rm $bedF");
  system("rm $bedFs");
}

# clean
$tempSamtools = $sample."_metrics.txt";
system("rm $tempSamtools");
$tempUD = $sample."_rem.bam.cnt";
system("rm $tempUD");

####################################################################################
sub make_plots{
  system("ngs.plot.r -G mm10 -R tss -C $outD/$rem_dup -O $plotsD/$plotTSS -D refseq -T $sample");
  system("ngs.plot.r -G mm10 -R genebody  -C $outD/$rem_dup -O $plotsD/$plotGEN -D refseq -T $sample");
   system("java -jar /apps/RH7U2/general/picard/2.8.3/picard.jar CollectInsertSizeMetrics I=$outD/$rem_dup MINIMUM_PCT=0.5 O=$sample.log H=$plotsD/$plotF W=1000");
  return();
}
####################################################################################
sub process_alignment{
  # step 3: clean
  system("grep -v 'chrM' $outD/$hits| samtools view -b -h -F 4 -f 0x2 - >$outD/$mapped");
  system("java -jar /apps/RH7U2/general/picard/2.8.3/picard.jar SortSam SO=coordinate I=$outD/$mapped O=$outD/$sorted");
  system("java -jar /apps/RH7U2/general/picard/2.8.3/picard.jar MarkDuplicates INPUT=$outD/$sorted OUTPUT=$outD/$rem_dup M=$metrics REMOVE_DUPLICATES=true");
  system("samtools index $outD/$rem_dup");
  system("samtools flagstat $outD/$rem_dup >$outD/$flagstat");
 
  
  # delete unrequired files
  system("rm $outD/$hits");
  system("rm $outD/$mapped");
  system("rm $outD/$sorted");
  return();
}

####################################################################################
sub run_fastqc{
  ($fastqF,$sample,$side,$fastqcD)=@_;

  # run fastqc
  $fastqcD_r = "$fastqcD/"."$sample"."_"."$side";
   make_path("$fastqcD_r");
   print "fastqc -o $fastqcD_r -f fastq $fastqF\n";
   system("fastqc -o $fastqcD_r -f fastq $fastqF");

  print "Run fastqc on $fastqF=> finished\n";
  return();
}

####################################################################################
sub collect_reads{
  my($crudeFASTQ) = @_;
  my(@files1) = `ls -1 $crudeFASTQ/$sample/*R1*`;
  my(@files2) = `ls -1 $crudeFASTQ/$sample/*R2*`;
  chomp(@files1);
  chomp(@files2);

  system("zcat $files1[0] > $fastq1_1");
  system("zcat $files2[0] > $fastq2_1");
  for($i=1;$i < $#files1+1;$i++){
     system("zcat $files1[$i] >> $fastq1_1");
    system("zcat $files2[$i] >> $fastq2_1");
   
  }
 
  system("gzip $fastq1_1");
  system("gzip $fastq2_1");

 return();
}
####################################################################################
sub make_nucleosome_free{
  my($bamFile) = @_;


  system("samtools view -H $bamFile > $tempF");
  system("samtools view $bamFile | awk -F \"\\t\" \'{if ((\$9>- 120) && (\$9< 120))  print \$_}\' >> $tempF");
  system("samtools view -h -b $tempF > $tempB");
  system("samtools index $tempB");
  
  # bam to bed- required for counting
  #system("bedtools bamtobed -i $tempB >$bedF");
  
  # generate flagstat report
  system("samtools flagstat $tempB > $flagstatR");
  
  # clean
  system("rm $tempF");
  
  # generate tdf file
  system("igvtools count -w 5 $tempB $tdf mm10");
  
 
 return();
}
####################################################################################
sub read_params{
  my($param_file) = @_;
  
  $cfg = new Config::Simple("$param_file");
  # general params
  $genome = $cfg->param("params.genome");
  $adaptor1 = $cfg->param("params.adaptor1");
  $adaptor2 = $cfg->param("params.adaptor2");
  $TSSfile = $cfg->param("params.TSS_file");
  $crude_reads_location= $cfg->param("params.crude_reads_location_for_merging");
  
  # run setup_run
  $collectReads = $cfg->param("setup_run.combine_fastq");
  $fastqc = $cfg->param("setup_run.fastqc");
  $trim_adapter= $cfg->param("setup_run.trim_adapter");
  $make_body=$cfg->param("setup_run.make_body");
  $make_plots=$cfg->param("setup_run.make_plots");
  $nucleosome_free=$cfg->param("setup_run.nucleosome_free");
  $call_peaks = $cfg->param("setup_run.call_peaks");
  $countTSS_reads=$cfg->param("setup_run.countTSS_reads");
  
  
   return();
}
####################################################################################
sub load_modules{
  do ('/apps/RH7U2/Modules/default/init/perl.pm');
  module('load fastqc');
  module('load jdk');
  module('load R');
  module('load ngsplot');
  module('load bowtie2');
  module('load picard-tools/1.119');
  module('load IGVTools');
  module('load samtools');
  module('load bedtools-gnu/2.27.1');
  module('load python');
  
   return();
}
