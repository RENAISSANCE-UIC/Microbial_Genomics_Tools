#!/bin/bash

# VARIABLES & FILES MANAGEMENT #

# Subsampling #
  minlength=100



# File names #
  RAWREADS="raw_reads"
  OUTPUT="processedReads2"
  subsampled="$OUTPUT/subsampled"
  samples="samples"


# Directories #
mkdir -p ${OUTPUT}/subsampled
mkdir -p results/{centrifugeWoL,centrifuge,kraken2,breseq}

# /FILES #

main_ReadProcessing(){

sampling
trimming
#subsampling
#centrifuge_wol
#centrifuge_default
#kraken2_default
#breseq_default
}

sampling(){

cd ${RAWREADS}

   ls *R1*.fastq.gz | sed s/_L001_R1_001.fastq.gz// > ../samples

cd -
}

trimming(){
# you will need to generate a text file called "samples" with the sample names.
# the switches you can adjust as necessary, but usable as is.

#files=$(ls -d $RAWREADS/*R1*.fastq.gz | sed s/_L001_R1_001.fastq.gz// | awk '{print $NF}')
cd $RAWREADS

echo "Processed Reads for " $PWD > processedReads2.txt

while read sampleID; do
echo "Trimming ${x}..."
/environments/builds/bbmap/38.90/bbmap/bbduk.sh -Xmx84g \
  in1=${sampleID}_L001_R1_001.fastq.gz \
  in2=${sampleID}_L001_R2_001.fastq.gz \
  out1=../$OUTPUT/${sampleID}_R1.fastq \
  out2=../$OUTPUT/${sampleID}_R2.fastq \
  t=16 \
  qtrim=rl \
  trimq=20 \
  minlen=$minlength \
  maxns=0 \
  rieb=t |& tee -a ${sampleID}_bbduk_rpt.txt

# Generate processed reads report
  echo -e -n ${sampleID} '\t' | tee >> processedReads2.txt
  sampleReads=$(grep "^@" ../$OUTPUT/${sampleID}_R?.fastq | wc -l)
  echo $sampleReads | tee >>processedReads2.txt


done < ../$samples

#  >> processedReads2/${x}_bbduk.log 2>&1

# copy this info into the results files for R
#cp processedReads2.txt results/centrifuge/
#cp processedReads2.txt results/centrifugeWoL/
#cp processedReads2.txt results/kraken2/
cd -
}

subsampling(){

cd $OUTPUT

   N=500000	#### Set the number of subsampled reads ####

while read y;
 do
echo "Subsampling ${y}..."
        seqtk sample -s100 ${y}_R1.fastq $N > $subsampled/${y}_500K_R1.fastq
        seqtk sample -s100 ${y}_R2.fastq $N > $subsampled/${y}_500K_R2.fastq
        cat $subsampled/${y}_500K_R1.fastq $subsampled/${y}_500K_R2.fastq > $subsampled/${y}_combined_1M.fastq

done <../$samples
cd -

}

centrifuge_wol(){

databaseWoL="/galaxy/biodatabases/web_of_life/centrifuge/WoLr1"
# add line to generate report in folder for database, numreads, etc.
echo "Centrifuge:WoL"

while read sampleID; do

echo "${sampleID}..."

centrifuge \
-x $databaseWoL \
-1 $OUTPUT/${sampleID}_R1.fastq \
-2 $OUTPUT/${sampleID}_R2.fastq \
-p 36 \
-S /dev/null \
--report-file ./results/centrifugeWoL/${sampleID}_centrifuge.txt

REPORT="./results/centrifugeWoL"

# Sort the output from centrifuge using the numReads
  (head -n 1 $REPORT/${sampleID}_centrifuge.txt && \
   tail -n +2 $REPORT/${sampleID}_centrifuge.txt | sort -nr -k5 -t$'\t') \
   > $REPORT/${sampleID}_centrifuge_sorted.txt


done < $samples
}

centrifuge_default(){

databaseDefault="/galaxy/biodatabases/centrifuge/09262018/p+h+v"
#add report in folder for database and numreads

while read sampleID; do

echo "Centrifuge:default ${sampleID}..."

### centrifuge, run subsampled reads

centrifuge \
-x $databaseDefault \
-1 $OUTPUT/${sampleID}_R1.fastq \
-2 $OUTPUT/${sampleID}_R2.fastq \
-p 16 \
-S /dev/null \
--report-file results/centrifuge/${sampleID}_centrifuge.txt

REPORT="results/centrifuge"

# Sort the output from centrifuge using the numReads
(head -n 1 $REPORT/${sampleID}_centrifuge.txt && \
 tail -n +2 $REPORT/${sampleID}_centrifuge.txt | \
 sort -nr -k5 -t$'\t')  > $REPORT/${sampleID}_centrifuge_sorted.txt


done < $samples

}

kraken2_default(){

db_default="/galaxy/biodatabases/kraken2-default"

export PATH="/environments/builds/kraken2/scripts:$PATH"
pwd
while read x; do

echo "kraken2 processing ${x}..."

# confidence changed 0.001-->0.05 on 2021-10-29; should result in lower sens
#  but higher spec. of low-abundance taxa

kraken2 --db $db_default \
 --threads 36 \
 --use-names \
 --confidence 0.05 \
 --report ./results/kraken2/${x}_k2_report.txt \
 --paired \
 $OUTPUT/${x}_R1.fastq \
 $OUTPUT/${x}_R2.fastq \
 --output ./results/kraken2/${x}.txt

cut -f1-3,6-8 ./results/kraken2/${x}_k2_report.txt > ./results/kraken2/${x}_stdK2_report.txt

done < $samples

}

breseq_default() {

eval "$(conda shell.bash hook)"
conda activate breseq_env

cd ${OUTPUT}       #processedReads2
#cd ${RAWREADS}	    #raw_reads

while read x; do
#conda activate breseq_env
breseq \
-r ../Ecoli_NIST0056.gbk \
-p \
-j 30 \
${x}_R1.fastq \
${x}_R2.fastq \
-o ../results/breseq/${x}
done <../$samples
cd -
#${x}_L001_R1_001.fastq.gz \
#${x}_L001_R2_001.fastq.gz \
#-o ../results/breseq_trimmed_r/${x}
#done <../$samples
#cd -

#-o ../results/breseq/${x}
#done <../$samples  	#done < $samples
#cd –
}


### Call the main list to execute ###
main_ReadProcessing
