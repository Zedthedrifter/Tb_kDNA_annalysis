#!/bin/bash
#download files from SRA with SRR accession numbers

export PATH=$PATH:/home/zed/script_collections:/home/zed/.local/bin/:/usr/bin:/bin
now=$(date '+%y-%m-%d-%H-%M')
log_file=log_${now}.txt
stats_file=stats_${now}.txt

function main {

printf "%s\t%s\t%s\t%s\n" 'read_file' 'total_reads' 'CSB3_containing_reads' 'percentage' >>  $stats_file
for acc in SRR7739651 SRR7739655 SRR7721317 SRR7721318	

do

prefetch $acc -O .
fastq-dump ${acc}.sra
gzip ${acc}.fastq
#clean up intermediates
rm ${acc}.fastq ${acc}.sra
#stats
fq=${acc}.fastq.gz
all_reads=$(grep -c '@SRR' ${fq})
csb3_reads=$(egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC' ${fq})
percentage=$(echo "scale=2 ; ${csb3_reads} / ${all_reads} * 100" | bc)
printf "%s\t%s\t%s\t%s\n" ${fq} ${all_reads} ${csb3_reads} ${percentage} >>  $stats_file
done

}


main


#SRR7721317 wild-type R1 [Hi-C I]
#SRR7721318 wild-type R2 [Hi-C I]
#SRR7739651 ATACseq wild-type, 10 million cells
#SRR7739655 ATACseq wild-type, 20 million cells
