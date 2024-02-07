#!
##############################
my($inSamples) = @ARGV;

$suffix="_TSScounts.txt";
$path = "7_TSS";
# samples
open(IN,"$inSamples")|| warn "no samples file\n";
@samples = <IN>;
chomp(@samples);
close(IN);

#read count
open(IN,"final_counts.txt")|| warn "no file with read count\n";
@readCount=<IN>;
chomp(@readCount);
close(IN);
system("rm final_counts.txt");


$fpkmN = 1;
# calculate FPKM factor
if($fpkmN == 1){
  # calculate the per million factor
  $i=0;
  foreach $readC (@readCount){
    $factor[$i] = $readC/1000000;
    $i++;
  }
}

foreach $file (@samples){
  $file1 = "$path/".$file."$suffix";
  open(IN,"$file1")|| warn "no $file\n";
  while($line = <IN>){
    chomp($line);
    
   ($chr,$start,$end,$name,$score,$strand,$count) = split(/\t/,$line);
    #($chr,$start,$end,$count) = split(/\t/,$line);
    $data{$chr}{$start}{$end}{$file} = $count;
    $data{$chr}{$start}{$end}{name} = $name;
    $data{$chr}{$start}{$end}{strand} = $strand;
    
  }
  close(IN);
}

### print
open(OUT,">TSS_counts_table.txt")|| warn "can not write to outfile\n";
print OUT "\t";
foreach $sample (@samples){
  print OUT "$sample\t";
}
print OUT "\n";

foreach $chr (keys %data){
  foreach $start (keys %{$data{$chr}}){
    foreach $end (keys %{$data{$chr}{$start}}){
      $key = "$chr"."_"."$start"."_"."$end";
       print OUT "$key\t$data{$chr}{$start}{$end}{name}\t$data{$chr}{$start}{$end}{strand}\t";
       foreach $file (@samples){
         print OUT "$data{$chr}{$start}{$end}{$file}\t";
       }
       if($fpkmN == 1){
         $i=0;
         $peakL = ($end-$start+1)/1000;
         foreach $file (@samples){
           $rpm = $data{$chr}{$start}{$end}{$file}/$factor[$i];
           $rpkm = $rpm/$peakL;
           print OUT "$rpkm\t";
         }
       }
       print OUT "\n";
    }
  }
}
close(OUT);