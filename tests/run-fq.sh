#!/bin/bash


TIME="/usr/bin/time -v"
FILES=/ccb/salz4-1/fbreitwieser/microbiome-pipeline/staging/CP_PT{[0-9],10}-*.fastq.gz

for J in 1 2; do
	for THREADS in 10 5 2 1; do 

		krakenhll --preload --db ../dbs/refseq-oct2017-k31 --report-file /dev/null --fasta <(printf ">A\nA")
		for B in $*; do 
			echo $B
			BN=`basename $B .fastq.gz` && echo $BN && BN="$BN.t$THREADS.j$J" 
			$TIME -o log/$BN.krakenhll.timing.log krakenhll --db ../dbs/refseq-oct2017-k31 --fastq --gzip --report-file report/$BN.krakenhll.report --threads $THREADS $B > kraken/$BN.krakenhll.kraken 2> log/$BN.krakenhll.log; 
		done

		kraken --preload --db ../dbs/refseq-oct2017-k31 --fasta <(printf ">A\nA")
		for B in $*; do 
			echo $B
			BN=`basename $B .fastq.gz` && echo $BN && BN="$BN.t$THREADS.j$J" 
			$TIME -o log/$BN.kraken.timing.log kraken --db ../dbs/refseq-oct2017-k31 --fastq --gzip  --threads $THREADS $B > kraken/$BN.kraken 2> log/$BN.kraken.log
			$TIME -o log/$BN.kraken-report.timing.log kraken-report --db ../dbs/refseq-oct2017-k31 kraken/$BN.kraken > report/$BN.kraken.report; 
		done

	done
done

