#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N VT
#$ -l h_rt=12:00:00
#$ -l h_vmem=12G
#$ -cwd

i=$SGE_TASK_ID
fvcf=$1
fped=$2
prefix=$3
#for i in 2 12
#for i in {1..22} X Y
#do
if [ "$i" == 23  ];
then 
i="X"
fi

if [ "$i" == 24 ]
then 
i="Y"
fi
ftrans=/home/nz2274/Resources/Gene/refGene_mRNA_protein_coding_hg37_Ensemble95_cannonical/refGene_mRNA_protein_coding_hg37_Ensemble95_4VT.$i.cannoical.bed;  #"/home/nz2274/Resources/Gene/cannoical_ensembl95_hg38/refGene_mRNA_protein_coding_hg38_Ensemble95_4VT.cannoical.$i.txt"
#ftrans="/home/nz2274/Resources/Gene/refGene_mRNA_protein_coding_hg38_chr/refGene_mRNA_protein_coding_hg38_Ensemble95_4VT.$i.anno.txt"
 #"/home/local/ARCS/nz2274/Resources/Gene/cannoical_ensembl95/refGene_mRNA_protein_coding_hg38_Ensemble95_4VT.cannoical.$i.txt"
fout=../result/$(basename $fped|sed 's/.ped//g').$i.$prefix.full
echo "perl VAT_binomial_trans.pl $fvcf $fped $ftrans $fout \n "
perl VAT_binomial_transcript.pl $fvcf $fped $ftrans $fout 
echo "done"
a=$(wc -l  $ftrans|awk '{print $1}')
b=$(wc -l $fout.REVEL.txt |awk '{print $1}')

if [ $a  != $b ]
then

echo "incomplete, only $b of $a transcripts were tested.\n"

fi 



#wait
