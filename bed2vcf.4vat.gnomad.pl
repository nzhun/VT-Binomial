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

sub generate {
  my ($ac_nfe,$hom_nfe,$total)=@_;
  my $het=$ac_nfe-2*$hom_nfe;
  my %gts=();
  while(scalar(keys %gts) <$het){
     $gts{int(rand($total))}=1;
  }

  while(scalar(keys %gts) <$hom_nfe){
      my $id=int(rand($total));
      if(!exists($gts{$id})){$gts{$id}=2;}
  }
  return (\%gts);

}

sub main {
      my $file=$_[0];
      my $fout=$file.".vcf";
      open IN,"less $file|" or die "cannot find $file\n";
      open OUT, ">$fout";
      open IDS,$_[1]; #"/home/local/ARCS/nz2274/PAH/PAH_10032017/src/EUR.case.gnomad.list"; #"/home/local/ARCS/nz2274/Resources/Control/InHouse_Control/src/EUR_control_4597.txt"; ## cases list
      my %inds=();
      my $col=0;
      my $header="#CHROM\tPOS\tID\tREF\tALT\tGene\tREVEL\tMCAP\tCADD_phred\tFORMAT";
	  my $total=0;
	  my $start=0;
      while(my $line=<IDS>){
		  
        chomp($line);
		my @sets=split(/\s+/,$line);
		my $id=$sets[0];
		#if($id =~ /CCHMC/){$start=$start+1;}
		if( $id =~ /gnomad/){
			$total=$total+1;
			if($start==0){$start=$col}
		}
        $inds{$id}=$col;
       # print $col."\t".$id."\t$start\t$total\n";
        $col=$col+1;
        $header=$header."\t".$id;

      }
     # print "\t$start\t$total\n";
      print OUT $header."\n";
      my %info_cols=();
      my $lastvar="";
      my $lastinfo="";
      my $lastscore="";
      my $lastgene="";
      my $info="";
      my $nr=0;
      my @refarr=();
      my $flag=0;
      my $ac_nfe=0;
      my $gt_str="";
      push @refarr, '0' foreach (1..scalar(keys %inds));
      my @curGT=@refarr;
      while(my $line=<IN>){
	      chomp($line);
          my @sets=split(/\t/,$line);
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
           if($var ne $lastvar){
               ## GT to string
                $gt_str=join("\t",@curGT);
               ## info to string
               if($flag==1){print OUT $lastvar."\t$lastgene\t".$lastscore."\tGT\t".$gt_str."\n";} #print "out\n";}
               else{
                 print $var."\t".$lastvar."\t".$sets[$info_cols{"proband"}]."\tNo\n";
               }

               @curGT=@refarr;
               $lastvar=$var;
               #$lastinfo=info2str(\@sets,\%info_cols);
               $flag=0;
           }


             {
               my $info_str=info2str(\@sets,\%info_cols);
               $lastvar=$var;
	           $lastgene=$sets[$info_cols{"Gene.refGene"}];
               #$lastinfo=$info_str;
		       my $cadd=$sets[$info_cols{"CADD_phred"}];
		       #if($cadd eq "."){$cadd=70;}
		       #$lastscore="1\t1\t$cadd";
	           if($sets[$info_cols{"ExonicFunc.refGene"}] !~ /non/){
		       #my $cadd=$sets[$info_cols{"CADD_phred"}];
		          if($cadd eq "."){$cadd=70;}
		          $lastscore="1\t1\t$cadd";
	          }elsif($sets[$info_cols{"REVEL"}] eq "."){
		 		    if($cadd eq "."){$cadd=0;}
		 			$lastscore="0\t0\t$cadd";
		      }else{
	                $lastscore=$sets[$info_cols{"REVEL"}]."\t".$sets[$info_cols{"MCAP"}]."\t".$sets[$info_cols{"CADD_phred"}];
	          }
           }      
		#    print $info_cols{"proband"}."\t"."\t".$sets[$info_cols{"proband"}]."\n";
           if(exists($info_cols{"proband"})){		
                if($sets[$info_cols{"proband"}] eq "."){
                 	my $ac_nfe=$sets[$info_cols{"AC_NFE"}]; ## homo_nfe,
                	my $hom_nfe=$sets[$info_cols{"Hom_NFE"}];
				#	print $start."\t".$total."\n";
                	my %pinds=%{generate($ac_nfe,$hom_nfe, $total)}; ### better to know homo or het
                	foreach my $ind(keys %pinds){
						 my $gt=1;
						 if($sets[0] ne "X") {$gt=$pinds{$ind};} 
                 	  	 $curGT[$ind+$start]=$gt;
						# print  "YYY:  ".($ind+$start)." :\t: ".($curGT[$ind+$start])." \t ".$total."\n";
               	  	}
				 	$flag=1;
          		 }  
				 if(exists($inds{$sets[$info_cols{"proband"}]})){
	 	   			my $gt=1;
	 	     		if( $sets[0] ne "X"  && $sets[$info_cols{"Is_homo"}] eq "T"){$gt=2}
             		$curGT[$inds{$sets[$info_cols{"proband"}]}]=$gt; #$sets[$info_cols{"Genotype"}];
       #      print "Y\t$lastvar\t".$inds{$sets[$info_cols{"proband"}]}."\t".$sets[$info_cols{"proband"}]."\n";
	   	  			$flag=1;
            	}
     	  	  
         
      		}

         }
        #if($nr>1){exit;}


      if($flag==1){print OUT $lastvar."\t$lastgene\t".$lastscore."\tGT\t".$gt_str."\n";}
      close IN;
      close OUT;
      close IDS;
}
main($ARGV[0],$ARGV[1]);
