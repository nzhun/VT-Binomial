for fvcf in ../data/Obesity_SPARK.0.0001.merge.bed.gz.vcf.gz  ../data/Obesity_SPARK.merge.bed.gz.vcf.gz  #../src/US_UK_SPARK.filter0.001.bed.vcf.gz  #../src/US_UK_SPARK.filter0.0001.bed.vcf.gz 

do
for fped in ../src/asso.EUR.ped   #../scr/asso.EUR.ped  #../src/APAH_US_UK_SPARK.ped   ../src/US_UK_SPARK.ped ../src/IPAH_US_UK_SPARK.ped 

do 
pref=$(basename $fvcf|sed 's/.bedi.gz.vcf.gz//g' )
qsub -t 1-24 qsub_VT.sh $fvcf $fped $pref

done

done 
