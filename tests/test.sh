#!/bin/bash

set -xeu

`dirname $0`/../install_krakenuniq.sh `pwd`/install
install/krakenuniq-download --dust --db DB refseq/archaea/Chromosome taxonomy
install/krakenuniq-build --db DB --threads 2 --kmer-len 21 --minimizer-len 12 --taxids-for-genomes --taxids-for-sequences
