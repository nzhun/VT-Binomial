#!/usr/bin/perl;
use warnings;
use strict;
#use Math::GSL::CDF qw/:binomial/;
#use Math::GSL::Randist qw/:binomial/;
use binomtest;
#use Math::CDF;
#use Statistics::R;
our $keycol=6;
our $start=0.2;
our $step=0.05;
our $NR=(1-$start)/$step+1;
print "Totally bins ".$NR."\n";
if(@ARGV<4){
	print "Input format: fin fped region outpredix\n";
	exit;
}

print "Input: ".join(" ",@ARGV)."\n";
my $fin=$ARGV[0]; #"/home/local/ARCS/nz2274/PAH/PAH_10032017/hg38/VAT/test.vcf.gz";
my $fped=$ARGV[1]; #;"/home/local/ARCS/nz2274/PAH/PAH_10032017/src/EUR.case_control.affected.ped";
my $chr=$ARGV[2];
my $prefix=$ARGV[3];
#our $R = Statistics::R->new();
main($fin,$fped,$chr,$prefix);
#$R->stopR();
sub generate_arr {
  my $cnt=$_[0];
  my $mp=$_[1];
  my @arr=();
  for(my $i=0;$i<$cnt;$i++){
     my $r=rand(1);
	 if($r>$mp){push(@arr,0);}else{push(@arr,1)}
  }

  return (\@arr);
}

#sub linear {
#	my $num=$_[0];
#	my $denum=$_[1];
#	my $z0=$num[$k]/sqrt($denum[$k]);
#	return ($z0); 
#}

sub clean {
	 my ($mat,$phe)=@_;
	 my @matrix2=@$mat;
	 my @phe2=@$phe;
	 my @newmat=();
	 my @newphe=();
	 for (my $i=0;$i<@matrix2;$i++){
		 if($matrix2[$i]==0){next}
		 push(@newmat,$matrix2[$i]);
		 push(@newphe,$phe2[$i]);
	 }
	 return (\@newmat,\@newphe);
}


sub Btest {
	my $num=$_[0];
	my $denum=$_[1];
	my $mp=$_[2];
#	print "ss\t".$num."\t".$denum."\t".$mp."\t\n";
	#print join("\t",@_)."  rece\n";
	my $r=-1;
	if(@_ <3){print "$num:$denum:$mp\nERROR\n";return ($r);}
	$r=binomtest::binomtest($num,$num+$denum,$mp);
#	my $gsl_p=gsl_ran_binomial_pdf($num,$mp,$num+$denum);
 #       my $gsl_cdf=1-gsl_cdf_binomial_P($num,$mp,$denum);
 #       $r=$gsl_p+$gsl_cdf;
	return ($r);
}

sub best_chose {
  my @denum=@{$_[0]};
  my @num=@{$_[1]};
  my $mp=$_[2];
  my $bestZ=1;
  my $bestT=$start;
  my $bestk=0;
  my @r=();
 # print join("\t",@denum)."\n".join("\t",@num)."\n";
  for(my $k=0;$k<$NR;$k++){
   # print "switch ".$k."\t"."\t$denum[$k]\t$mp\n";
    if($denum[$k]==0 && $num[$k]==0){last;}
	#print "1\t$k\t".$num[$k]."\t".$denum[$k]."\t".($start+$k*$step)."\t".($denum[$k]+$num[$k])."\t".$mp."\n";
	#print "1\t$k\t".$num[$k+1]."\t".$denum[$k+1]."\t".($start+$k*$step)."\t".$mp."\n";
	#exit;
	#print "start $k\n";
	my $z0=Btest($num[$k],$denum[$k],$mp);
	#print "end\n";
	#exit;
     #my $z0=$num[$k]/sqrt($denum[$k]);
	# print "swith: For best\t".$num[$k]."\t".$denum[$k]."\t".$k."\t".$z0."\t"."\n";
     if($z0>$bestZ){next;}
     $bestZ=$z0;
     $bestT=$start+$k*$step;
	 $bestk=$k;
    #  print "swith: For best\t".$k."\t".$z0."\t".$bestZ."\t".$bestT."\n";
      ## output the best threshold and best Z for each gene
    #}

   }
   #print "swith: For best\t"."\t".$bestZ."\t".$bestT."\n";
#   if($bestZ==1){return \@r;}
if($num[$bestk]==0){$num[$bestk]+=0.5;}   
@r=($bestZ,$bestT,$denum[$bestk],$num[$bestk]);
   return \@r;
}

sub get_best {
   # print "b1".localtime."\n";
 	my $cutoff=0;
	my $p=1;
	my ($new_mat,$new_phe,$mp)=@_;
	my @mat=@$new_mat;
	my @phe=@$new_phe;
	my @cs=();
	my @ct=();
	push @cs, 0 foreach (0..($NR-1));
	push @ct, 0 foreach (0..($NR-1));
	for(my $i=0;$i<@mat;$i++){
		my $hit=int(($mat[$i]-$start)/$step);
		for(my $h=0;$h<$hit+1;$h++){
				 if($phe[$i]==0){
				 	$ct[$h]=$ct[$h]+1;
				 }else{
				 	$cs[$h]=$cs[$h]+1;
				 }
		 	}

		}
	#	print $mp."rrr\n";
	 my @r=@{best_chose(\@ct,\@cs,$mp)};
	#	  print "b2".localtime."\n";
	# print join("\t",@r)."\n";
     return (\@r);
}

sub deal_gene {
	 my ($matd,$phed,$mp)=@_; 
	 my @matrix=@$matd;
	 my @phe=@$phed;
	 my ($new_mat,$new_phe)=clean(\@matrix,\@phe);
	 my $LEFT=@$new_phe;
	 if($LEFT==0){return 1;}
	# print "reduce to ".@$new_mat."\t".@$new_phe."\n";
	 ## compute the best cutoff and p-value 
	 my ($p,$cutoff,$N1,$N2)=@{get_best($new_mat,$new_phe,$mp)};
	 if($N1==0){$N1=0.05;}
	 my $OR=($N2/$N1)*(1/$mp-1);	
# print $mp."\t".$p."\t".$cutoff."\n";
	
	 my @rs=($cutoff,$p);
	 my $T1=1000;
	 my $Te=$T1;
	 my $Ts=0;
	 my $pcount=0;
	 LOOP:
	 while($Te < 10000000){
		# print $Te."\n";
		 for(my $ccct=$Ts;$ccct<$Te;$ccct++){
			 ##permut new_phe; ## resign phenotype based on the freq, ignore the acutal number of case-control in new_phe
		#	   print "b3".localtime."\n";
			  my $new_phe2=generate_arr($LEFT,$mp);
		#	  print "b5".localtime."\n";
			  my ($p2,$cutoff2,$n1,$n2)=@{get_best($new_mat,$new_phe2,$mp)};
		#	  print $ccct."\t".$p2."\n";
			  if($p2 >$p){
				  next;
			  }
			  $pcount+=1;
		 }
		 my $permup=($pcount+1)/($Te+1);
		 if($permup >5/($Te+1)){
			 push(@rs,$permup);
			 push(@rs,$Te);
			 push(@rs,$OR);
			 push(@rs,$N2);
			 push(@rs,$N1);
			 #print OUT $lastgene."\t".$cutoff."\t".$p."\t".$permup;
			 last LOOP;
	      }else{
			  $Ts=$Te;
			  $Te=$Te*10;
	      }
    }
	if($Te >10000000 ||@rs<4){
		 my $permup=($pcount+1)/($Te+1);
		 push(@rs,$permup);
		 push(@rs,$Te);
		 push(@rs,$OR);
		 push(@rs,$N2);
		 push(@rs,$N1);
	 	#print OUT $lastgene."\t".$cutoff."\t".$p."\t".$permup;
	
	}  
	
	return \@rs;
}
sub gen_best {
        my $fin=$_[0];
        my $chr=$_[1];
        #my $fout=$_[1];
        my @phe=@{$_[2]};
        my $mp=$_[3];
        my $keycol=$_[4];
		my @includes=@{$_[5]};
		my $fout=$_[6];
        my $lastgene="";
		my @matrix=();
		push @matrix, 0 foreach (0..(@phe-1));
		my @init_matrix=();
		push @init_matrix, 0 foreach (0..(@phe-1));
        my %result=();
        my $bestT=0;
        my $bestZ=1;
        my $N=10; ##Gt column
		my $call=0;
        print "score: ".$start."\t".$step."\t$NR\n";
        open IN, "tabix $fin $chr|sort -k6,6d|";
		open OUT, " > $fout";
	
        while(my $line=<IN>){
          chomp($line);
          my @sets=split(/\t+/,$line);
          my @c=@sets[@includes];
          my $site_num=0;
          my $site_denum=0;
          my $revel=$sets[$keycol];
          my $gene=$sets[5];
		  if($gene =~ /;/){next}
          if($revel eq "." || $revel<$start){next;}
		 # print $gene."\t".$revel."\n";
	      if($lastgene eq ""){$lastgene=$gene}
	      if ($gene ne $lastgene){
			#  print $gene."\t".$lastgene."\n";
			  if($call ==1){
				 ## clean matrix,cleam phe 
			      my ($cutoff,$p,$permutp,$Te,$OR,$Nc,$Ns)= @{deal_gene(\@matrix,\@phe,$mp)};
				  print OUT $lastgene."\t".$cutoff."\t".$p."\t".$permutp."\t".$Te."\t".$OR."\t$Nc\t$Ns\n";
				  ## need to clean the array
				  @matrix=@init_matrix;
	        	  $lastgene=$gene;
				  $call=0;
			#  exit;
	      }
	     }
		  for(my $i=0;$i<@c;$i++){
			  if($c[$i]==0){next;}
			  if($matrix[$i] <$revel) {$matrix[$i]=$revel}
			  $call=1;
		  }		 
	 
	  
        
             $lastgene=$gene;
        	 # last;
          
	  }
          close IN;
      #    print @denum."\t".join("\t",@denum)."\n num\t".join("\t",@num)."\n";	  
		  if($call ==1){
		      my ($cutoff,$p,$permutp,$Te,$OR,$Nc,$Ns)= @{deal_gene(\@matrix,\@phe,$mp)};
			  print OUT $lastgene."\t".$cutoff."\t".$p."\t".$permutp."\t".$Te."\t".$OR."\t$Nc\t$Ns\n";
	 	   }
		   close OUT;
         # return \%result;
}
#close OUT;



sub main {
    print "input format: inputFile pedfile(case+control)  region  OutPut_prefix\n ";
	print "input format pedfile: sampleName	Affected_state(1 or 2)\n ";
	if(@_<4){print "check the input parameters please!\n";exit;}
    my ($fin,$fped,$chr,$prefix)=@_;
	 
	my %ped=();
	my @phe=();
	open PED, "$fped";
	while(my $line=<PED>){
	   chomp($line);
	   my @sets=split(/\t+/,$line);
	   $ped{$sets[0]}=$sets[1]-1;
	 }
	close PED;



	open IN,"tabix $fin  -H |" or die "$fin cannot find!\n";
	print "tabix $fin $chr -h\n";#open OUT,">$fout";
	my $line=<IN>;
	chomp($line);
	my @sets=split(/\t+/,$line);
	my $outkey=$sets[$keycol];
	my $fout="$prefix.$chr.$outkey.txt";
	close IN;
	my $N=10;
	my $mp=0;
	my @includes=();
	#my $PC=0;
#	  print @sets."\n";
	for(my $i=$N;$i<@sets;$i++){
	    if(exists($ped{$sets[$i]})){
	      push(@phe,$ped{$sets[$i]});
		  push(@includes,$i);
           
	    }else{
	      next;
	    }
	   #print $sets[$i]."\t".$ped{$sets[$i]}."\n";
	    $mp+=$ped{$sets[$i]};
	  }

	print "Total samples:".@phe.", cases: $mp\n";
	$mp=$mp/(@phe);

	my $permut=1000;#1000; #00;
	#print $mp."\t".join(":",@phe)."\n";
	#exit;
	gen_best($fin,$chr,\@phe,$mp,$keycol,\@includes,$fout);

}
