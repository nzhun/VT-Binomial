#!/usr/bin/perl;
use warnings;
use strict;
#use Math::GSL::CDF qw/:binomial/;
#use Math::GSL::Randist qw/:binomial/;
use binomtest;

#use Math::CDF;
#use Statistics::R;

#$R->stopR();
our %self_dist=();
our $keycol=6;
our $start=0.1;  ## mimum interested variants' revel 
our $step=0.02; ## step of revel to lookup
our $NR=(1-$start)/$step+1; ## based on start and step, the total number of bins can be calculated.

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

sub selectN {
	my ($N, $total)=@_;
#	my %pool=%$admap;
	my %sm=();
	my @arr= (0) x $total;
	while(scalar(keys %sm) <$N ){
		my $k=rand($total);
		#print "$k selected\n";
		 $arr[$k]=1;
		 $sm{$k}=1;
	 }
	 #print $N."\t".(keys %sm)." ".$total."\n";
	return (\@arr);
}


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
	if(exists($self_dist{$num.":".$denum})){
		$r=$self_dist{$num.":".$denum}
	}else{
		$r=binomtest::binomtest($num,$num+$denum,$mp);
		$self_dist{$num.":".$denum}=$r;
	}
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
	#print "aks ".@phe."\n";
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
	# if($p >0.1){return (\@rs);}
	 my $T1=100;
	 my $Te=$T1;
	 my $Ts=0;
	 my $pcount=0;
	 my $MAXT=10000000;
	 LOOP:
	 while($Te < $MAXT+1){
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
		 #my $permup=($pcount+1)/($Te+1);
		 if($pcount >9||$Te==$MAXT){
			 #push(@rs,$permup);
			 #push(@rs,$Te);
			 #push(@rs,$OR);
			 #push(@rs,$N2);
			 #push(@rs,$N1);
			 #print OUT $lastgene."\t".$cutoff."\t".$p."\t".$permup;
			 last LOOP;
	      }else{
			  $Ts=$Te;
			  $Te=$Te*10;
			  print "increase to ".$Te."\n";
	      }
    }
#	if($Te >10000000 ||@rs<4){
		 my $permup=($pcount+1)/($Te+1);
		 push(@rs,$permup);
		 push(@rs,$Te);
		 push(@rs,$OR);
		 push(@rs,$N2);
		 push(@rs,$N1);
	 	#print OUT $lastgene."\t".$cutoff."\t".$p."\t".$permup;

#	}

	return \@rs;
}
sub gen_best {
 #       my $fin=$_[0];
#	my $gene=
    #    my $chr=$_[1];
        #my $fout=$_[1];
   #     my @phe=@{$_[2]};
  #      my $mp=$_[3];
 #       my $keycol=$_[4];
#		my @includes=@{$_[5]};
#		my $fout=$_[6];
	my ($fin,$gene,$chr,$adphe,$mp,$keycol,$adinclude)=@_;
	my @includes=@$adinclude;
	my @phe=@$adphe;
	#print "come ".join(",",@phe).":::".@phe."\n";
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
	my $ostr="";
        #ggprint "$gene\t$chr\n";
        open IN, "tabix $fin $chr |";
	#open OUT, " >> $fout";

        while(my $line=<IN>){
          chomp($line);
	  if($line =~ /^#/){next;}
          my @sets=split(/\t+/,$line);
          my @c=@sets[@includes];
          my $site_num=0;
          my $site_denum=0;
          my $revel=$sets[$keycol];
          my $gene=$sets[5];
	#	  if($gene =~ /;/){next}
          if($revel eq "." || $revel<$start){next;}
		  for(my $i=0;$i<@c;$i++){
			  if($c[$i]==0){next;}
			  if($matrix[$i] <$revel) {$matrix[$i]=$revel}
			  $call=1;
		  }
	  }
          close IN;
      #    print @denum."\t".join("\t",@denum)."\n num\t".join("\t",@num)."\n";
		  if($call ==1){
	#	      print "still : ".@phe."\n";
		      my @ret=@{deal_gene(\@matrix,\@phe,$mp)};
		      if(@ret <5){print "$gene went wrong!\n";return($ostr);}
		      my ($cutoff,$p,$permutp,$Te,$OR,$Nc,$Ns)= @ret; #{deal_gene(\@matrix,\@phe,$mp)};
	#		 if($permutp <= 0.05) {
			 return ($gene."\t".$cutoff."\t".$p."\t".$permutp."\t".$Te."\t".$OR."\t$Nc\t$Ns");
	 #	   	}
		}
		   #close OUT;
		  return ($ostr);
}
#close OUT;



sub main {
    print "input format: inputFile pedfile(case+control)  region  OutPut_prefix\n ";
	print "input format pedfile: sampleName	Affected_state(1 or 2)\n ";
	if(@_<4){print "check the input parameters please!\n";exit;}
    my ($fin,$fped,$trans,$prefix)=@_;

	my %ped=();
	my @phe=();
	open PED, "$fped";
	while(my $line=<PED>){
	   chomp($line);
	   if($line =~ /^#/){next;}
	   my @sets=split(/\t+/,$line);
	   $ped{$sets[0]}=$sets[@sets-1]-1;
	 }
	close PED;



	open IN,"tabix $fin  -H |" or die "$fin cannot find!\n";
	print "tabix $fin -H\n";#open OUT,">$fout";
	my $line=<IN>;
	chomp($line);
	my @sets=split(/\t+/,$line);
	my $outkey=$sets[$keycol];
	my $fout="$prefix.$outkey.txt";
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
	my $Ncase=$mp;
	$mp=$mp/(@phe);

	my $permut=1000;#1000; #00;
	#print $mp."\t".join(":",@phe)."\n";
	open OUT, ">$fout" ;
	open TRIN,$trans or die "$trans cannot find!\n";
	while(my $tl=<TRIN>){
		chomp($tl);
		my @tsets=split(/\s+/,$tl);
		my $gene=$tsets[0];
		my $region=join("\t",@tsets[1..(scalar(@tsets)-1)]);
		my $st=gen_best($fin,$gene,$region,\@phe,$mp,$keycol,\@includes);
		if($st ne "") {print OUT $st."\n";}
	}
	close	TRIN;
	close OUT;
}

## /home/local/ARCS/nz2274/PAH/PAH_10032017/Result/Type1_Error/data/refGene_mRNA_protein_coding_hg37_Ensemble95_4VT.cannoical.bed

print "Warning: Before runing the script, please add the binomtest.pm to your perl library\n";

if(@ARGV<4){
	print "more input is required\n
	Input format: fin fped transcript outpredix\n";
	exit;
}



print "Warning: the input required specific format, it is tab-delimited\n";
print "Warning: it has to include: #CHR,POS,END,REF,ALT,GENENAME,REVEL,MCAP,CADD,FORMTAT,INDIVIDUAL.......\n";
print "Warning: please see the example input in the test file\n";

print "Warning: the input pedigree file must have two columns, tab-delimited,ID\tcoded_ohenotype.  1=healthy, 2=affected\n";
print "Warning: the input transcript file must follow the format: Name\\tregion1 region2......\n";

print "Warning: please see the example in the test folder, before runing the script\n";

my $fin=$ARGV[0]; 
if ( !-f "$fin.tbi"){
	print "Error: the input file has to be tabixed file\n";
	exit;
}
print "Totally bins ".$NR."\n";
print "Input: ".join(" ",@ARGV)."\n";

# required tabixed .gz file, it is required specific format, tab-delimited #CHR,POS,END,REF,ALT,GENENAME,REVEL,MCAP,CADD,FORMTAT,INDIVIDUAL....... 
# Genotype code: 0-> reference, 1-> heterozygous 2-> homozygous 
my $fped=$ARGV[1]; 
## required ID and phenotype 1-> healthy 2-> affected
my $ftrans=$ARGV[2];
## transcript region file: Name\tchr1:pos1-pos2 chr1:pos3-pos4.....
my $prefix=$ARGV[3];
## prefix for the output

main($fin,$fped,$ftrans,$prefix);