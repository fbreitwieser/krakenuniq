#!/usr/bin/env bash

set -xeu

case "$OSTYPE" in
  darwin*)  TIME="gtime -v" ;; 
  *)        TIME="usr/bin/time -v" ;;
esac

DB=$1
DB_BN=`basename $DB`
shift

mkdir -p log kraken report

for J in 1 2; do
	for THREADS in 4 2 1; do 

		krakenhll --preload --db $DB --report-file /dev/null --fasta <(printf ">A\nA") > /dev/null
		for B in $*; do 
			echo $B
			if [[ "$B" == *.gz ]]; then
				BN=`basename $B .fastq.gz`
				PARAM="--gzip"
			else
				BN=`basename $B .fastq`
				PARAM=
			fi
			echo $BN && BN="$BN.$DB_BN.t$THREADS.j$J" 
			$TIME -o log/$BN.krakenhll.timing.log krakenhll --db $DB --fastq $PARAM --report-file report/$BN.krakenhll.report --threads $THREADS $B > kraken/$BN.krakenhll.kraken 2> log/$BN.krakenhll.log; 
		done

		kraken --preload --db $DB --fasta <(printf ">A\nA") > /dev/null
		for B in $*; do 
			echo $B
			if [[ "$B" == *.gz ]]; then
				BN=`basename $B .fastq.gz`
				PARAM="--gzip"
			else
				BN=`basename $B .fastq`
				PARAM=
			fi
			echo $BN && BN="$BN.$DB_BN.t$THREADS.j$J" 
			$TIME -o log/$BN.kraken.timing.log kraken --db $DB --fastq $PARAM  --threads $THREADS $B > kraken/$BN.kraken 2> log/$BN.kraken.log
			$TIME -o log/$BN.kraken-report.timing.log kraken-report --db $DB kraken/$BN.kraken > report/$BN.kraken.report; 
		done

	done
done

