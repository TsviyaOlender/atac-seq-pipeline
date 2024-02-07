#!
##############################
use Config::Simple;
use File::Path qw(make_path remove_tree);
##############################################################
##### Before INSALLATION: update the location of $pipelinePATH
$pipelinePATH="/home/labs/bioservices/bareket/Yael/ATACscript/ATACpipeline_V2";
##############################################################


my($inSamples) = @ARGV;

$reportsPath = "8_reports";
$alignedPath = "3_align";
$nucleosomeFreePath = "5_nucleosome_free";

open(IN,"$inSamples")|| die "can not locate file with samples names\n";
@samples=<IN>;
chomp(@samples);
close(IN);


# prepare folder with all reports
make_path("$reportsPath");
foreach $sample (@samples){
  # get num of crude reads
   $file = "$sample"."_runlog.txt";
   system("mv $file $reportsPath/");  

}
foreach $sample (@samples){
  # get num of crude reads
   $file = "$sample"."_runlog.txt";
    $lineToParse = `grep \'Total read pairs processed\' $reportsPath/$file`;
    (@data) = split(/\s+/,$lineToParse);
    $data[4] =~s/\,//g;
    
    $qual{$sample}{crude} = $data[4]*2;
    
  # after trimming
      $lineToParse = `grep \'Pairs written (passing filters)\' $reportsPath/$file`;
    (@data) = split(/\s+/,$lineToParse);
    $data[4] =~s/\,//g;
    $qual{$sample}{trim} = $data[4]*2;
    
   # alignment report
   $file = $sample.".flagstat";
    $lineToParse= `grep \'in total (QC-passed reads\' $alignedPath/$file`;
    ($mapped) = $lineToParse =~/^(\d+)/;
    $qual{$sample}{mapped} = $mapped;
    
    # nucleosome free
   $file = $sample."_free.flagstat";
    $lineToParse= `grep \'in total (QC-passed reads\' $nucleosomeFreePath/$file`;
    ($mapped) = $lineToParse =~/^(\d+)/;
    $qual{$sample}{nucleosomefree} = $mapped;
    
    
}


### print summary
print "sample\tcrudeReads\tpassedCutAdapt\tMapped\tNucleosomeFree\n";
foreach $sample (@samples){
  print "$sample\t$qual{$sample}{crude}\t$qual{$sample}{trim}\t$qual{$sample}{mapped}\t$qual{$sample}{nucleosomefree}\n";
}

## save counts for temporary file
open(OUT,">final_counts.txt")|| warn "can not save counts\n";
foreach $sample (@samples){
  print OUT "$qual{$sample}{nucleosomefree}\n";
}
close(OUT);

system("perl $pipelinePATH/collect_bedtools_count.pl $inSamples");