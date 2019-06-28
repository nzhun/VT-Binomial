#!/usr/bin/perl;

use strict;
use warnings;

sub info2str {
  my @sets=@{$_[0]};
  my %info_cols=%{$_[1]};
#  print $var."\n";
#  print join("\t",%info_cols)."\n";
  my $str="";
  foreach my $key(keys %info_cols){
  #  print $key."\t".$info_cols{$key}."\t".$sets[$info_cols{$key}]."\n";
    $str.=$key."=".$sets[$info_cols{$key}].";";
  }
  return $str;
}

sub main {
      my $file=$_[0];
      my $fout=$file.".recode.vcf";
      open IN,"zless $file|";
      open OUT, ">$fout";
      open IDS, $_[1] ; #"/home/local/ARCS/nz2274/Resources/Control/InHouse_Control/src/EUR_control_4597.txt"; ## cases list
      my %inds=();
      my $col=0;
      my $header="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
      while(my $line=<IDS>){
        chomp($line);
	my @sets1=split(/\s+/,$line);
	my $id=$sets1[0];
        $inds{$id}=$col;
    #    print $col."\t".$id."\n";
        $col=$col+1;
        $header=$header."\t".$id;

      }

      print OUT $header."\n";
      my %info_cols=();
      my $lastvar="";
      my $lastinfo="";
      my $info="";
      my $nr=0;
      my @refarr=();
      my $flag=0;
      my $gt_str="";
      push @refarr, '0/0' foreach (1..scalar(keys %inds));
      my @curGT=@refarr;
      while(my $line=<IN>){
	  chomp($line);
          my @sets=split(/\t+/,$line);
          if($nr==0){
            for(my $j=0;$j<@sets;$j++){
              my $item=$sets[$j];
              $info_cols{$item}=$j;
            }
            $nr=$nr+1;
             next;
          }

           my $var=join("\t",@sets[0..4]);
           if($lastvar eq ""){

             $lastvar=$var;
           }
           if($var eq $lastvar){
               my $info_str=info2str(\@sets,\%info_cols);
               $lastvar=$var;
               $lastinfo=$info_str;
           }else{
               ## GT to string
                $gt_str=join("\t",@curGT);
               ## info to string



               if($flag==1){print OUT $lastvar."\t100\t.\t".$lastinfo."\tGT\t".$gt_str."\n";} #print "out\n";}
            #   else{
            #     print $var."\t".$lastvar."\t".$sets[$info_cols{"proband"}]."\tNo\n";
              # }

               @curGT=@refarr;
               $lastvar=$var;
               $lastinfo=info2str(\@sets,\%info_cols);
               $flag=0;
           }
      #     print $info_cols{"Genotype"}."\t".$info_cols{"proband"}."\t".$sets[$info_cols{"Genotype"}]."\t".$inds{$info_cols{"proband"}}."\n";
           if(exists($inds{$sets[$info_cols{"proband"}]})){
             $curGT[$inds{$sets[$info_cols{"proband"}]}]=$sets[$info_cols{"Genotype"}];
#             print "Y\t$lastvar\t".$inds{$sets[$info_cols{"proband"}]}."\t".$sets[$info_cols{"proband"}]."\n";
             $flag=1;
            # exit;
        #  print $sets[$info_cols{"proband"}]."\tY\n";
      }else{
 #            print $sets[$info_cols{"proband"}]."\t"."Not Eur\n";
           }

         }
        #if($nr>1){exit;}


      if($flag==1){print OUT $lastvar."\t100\t.\t".$lastinfo."\tGT\t".$gt_str."\n";}
      close IN;
      close OUT;
      close IDS;
}
main($ARGV[0],$ARGV[1]);
