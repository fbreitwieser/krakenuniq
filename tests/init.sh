
## Install KrakenU locally into install/
../install_kraken.sh `pwd`/install

## Download taxonomy and genomic data into data/
install/krakenu-download --db data -R --include-viral-neighbors taxonomy refseq/archaea refseq/bacteria refseq/viral/Any

