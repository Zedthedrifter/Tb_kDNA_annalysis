#!/bin/bash
#SNPs calling
export PATH=$PATH:/home/zed/script_collections:/home/zed/.local/bin/:/usr/bin:/bin
now=$(date '+%y-%m-%d-%H-%M')
log_file=log_${now}.txt
stats_file=stats_${now}.txt
function komics_assembly {

log_file=$1
r1=$2
r2=$3
strain=$4
stats_file=$5
kmers=$(echo 99 119 129 149)

  echo -e "start KOMICS assembly of ${strain} with kmer size ${kmers}\nr1= ${r1}\nr2= ${r2}" >> ${log_file}
  #assemble minicircle: higher k
  for k in ${kmers}
  do
    komics assemble --threads 16 --kmin ${k} --kmax ${k} ${strain}_k${k} ${r1} ${r2}
    #circularize
    komics circularize ${strain}_k${k} tmp.${strain}_k${k}.csb3contigs.fasta
    #collect all contigs for polishing
    cat  tmp.${strain}_k${k}.circularized.fasta >> tmp.${strain}.circularized.fasta
  done
  #extract circularized contigs
  extract_circ.py tmp.${strain}.circular.fasta tmp.${strain}.circularized.fasta
  komics polish --minidentity 95 ${strain} tmp.${strain}.circular.fasta
  rm -r tmp* *maxicir*
  minicircles=$(grep -c 'circular' ${strain}.minicircles.fasta)
  printf "Total minicircle found: \n%s \n" ${minicircles} >> ${log_file}

#completeness assessment
fasta_extend.py ${strain}.minicircles.extended.fasta ${strain}.minicircles.fasta 150
bowtie2-build -f ${strain}.minicircles.extended.fasta mini_ref
/usr/bin/bowtie2-align-s --wrapper basic-0 -x mini_ref -p 20 --threads 16 --phred33 --very-sensitive-local --quiet --time --rg-id ${strain} --rg SM:${strain} --rg PL:'ILLUMINA' -S map_mini.bam -1 ${r1} -2 ${r2}
  all_reads=$(samtools view map_mini.bam |wc -l)
  mapped=$(samtools view -F 4 map_mini.bam |wc -l)
  mapped_q=$(samtools view -F 4 -q 10 map_mini.bam|wc -l)
  percentmapped=$(echo "scale=2 ; ${mapped} / ${all_reads} * 100" | bc)
  percentmapped_q=$(echo "scale=2 ; ${mapped_q} / ${all_reads} * 100" | bc)
  printf "Total reads:\n%s\nMapped reads:\n%s %s%%\nMapped reads with quality > 10:\n%s %s%%\n" ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q} >> ${log_file}
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s" ${strain} ${minicircles} ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q}>>  $stats_file
 #CSB3 reads
  all_reads=$(samtools view map_mini.bam |egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
  mapped=$(samtools view -F 4 map_mini.bam |egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
  mapped_q=$(samtools view -F 4 -q 10 map_mini.bam |egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
  percentmapped=$(echo "scale=2 ; ${mapped} / ${all_reads} * 100" | bc)
  percentmapped_q=$(echo "scale=2 ; ${mapped_q} / ${all_reads} * 100" | bc)
  printf "Total CSB3 reads:\n%s\nMapped CSB3 reads:\n%s %s%%\nMapped CSB3 reads with quality > 10:\n%s %s%%\n" ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q} >> ${log_file}
  printf "\t%s\t%s\t%s\t%s\t%s\n" ${all_reads} ${mapped} ${percentmapped} ${mapped_q} ${percentmapped_q} >>  $stats_file
  echo -e "\n\n---------------------------\n\n" >> ${log_file}
  
rm map_mini.bam *bt2 ${strain}.minicircles.extended.fasta
}

function metadata_process {

stats_file=$1
meta=$2

#make duplicates
cp $meta tmp.meta.txt; meta=tmp.meta.txt; cp $stats_file tmp.stats.txt; stats_file=tmp.stats.txt
#process stats file
head -1 $stats_file > tmp1.txt; sed -i 1d $stats_file; cat $stats_file|sort -k1,1 > tmp3.txt
#process meta file
head -1 $meta|cut -d , -f 1,2,3,4,5,6,7,8,9,10,11,12,13 --output-delimiter=$'\t' > tmp2.txt; sed -i 1d $meta; cat $meta|sort -k1,1|cut -d , -f 1,2,3,4,5,6,7,8,9,10,11,12,13 --output-delimiter=$'\t' > tmp4.txt
#paste together
paste tmp1.txt tmp2.txt > stats_meta_${now}.txt; paste tmp3.txt tmp4.txt >> stats_meta_${now}.txt
#clean up
rm tmp*

}

function main {

#list input data
printf "strain\tminicircle_nb\ttotal_reads\tmapped_reads\tmapped_reads_percentage\tmapped_q10\tmapped_q10_percentage\ttotal_csb3\tmapped_csb3\tmapped_csb3_percentage\tmapped_csb3_q10\tmapped_csb3_q10_percentage\n" >> $stats_file
for r1 in $(ls /home/zed/disk3/Fre_WGS/se*/*.reads1.fq.gz)
  do
  strain=$(echo ${r1}|awk -F'/' '{print $NF}'); strain=${strain/.reads1.fq.gz/}
  r2=${r1/.reads1.fq.gz/.reads2.fq.gz}
#run komics assembly and completeness assessment
  komics_assembly ${log_file} ${r1} ${r2} ${strain} $stats_file
done
#combine metadata
metadata_process stats_22-02-11-14-52.txt /home/zed/disk1/Tb_Fre_WGS/common_ref/metadata.csv
}


main
