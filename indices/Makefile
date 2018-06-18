
all: nt-db

THREADS=10
DODUST=1

nt-db: nt/database.kdb

.PHONY: all

mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
TAXONLIST=$(mkfile_path)/microbial-taxa-incl-helminths.tsv

nt:
	mkdir -p nt

nt/nt.gz: nt
	curl -o nt/nt.gz ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz

nt/nt.fna: nt/nt.gz
	gunzip -c nt/nt.gz > nt/nt.fna

nt/taxonomy/accession2taxid:
	mkdir -p nt/taxonomy/accession2taxid

nt/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz: nt/taxonomy/accession2taxid
	curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz > nt/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

nt/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz: nt/taxonomy/accession2taxid
	curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz > nt/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz

nt/taxonomy/taxdump.tar.gz: nt/taxonomy/accession2taxid
	curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz > nt/taxonomy/taxdump.tar.gz

nt/taxonomy/names.dmp: nt/taxonomy/taxdump.tar.gz
	tar --directory nt/taxonomy -mvxf nt/taxonomy/taxdump.tar.gz names.dmp

nt/taxonomy/nodes.dmp: nt/taxonomy/taxdump.tar.gz
	tar --directory nt/taxonomy -mvxf nt/taxonomy/taxdump.tar.gz nodes.dmp

nt/nt-microbial.fna: nt/taxonomy/nodes.dmp nt/taxonomy/names.dmp $(TAXONLIST) nt/nt.fna nt/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz nt/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
	time perl $(mkfile_path)/filter-nt.pl -m nt/nt.fna.map -a nt/nt-filtered-acs.txt --nodes nt/taxonomy/nodes.dmp --names nt/taxonomy/names.dmp -t $(TAXONLIST) -f nt/nt.fna nt/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz nt/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz > nt/nt-microbial.fna.tmp
	mv nt/nt-microbial.fna.tmp nt/nt-microbial.fna

ifeq ($(DODUST), 1)
DB_F=nt/library/nt-microbial.fna
else
DB_F=nt/library/nt-microbial-dustmasked.fna
endif

nt/library/nt-microbial.fna: nt/nt-microbial.fna
	cp -v nt/nt-microbial.fna nt/library/nt-microbial.fna

nt/library/nt-microbial-dustmasked.fna: nt/nt-microbial.fna
	mkdir -p nt/library
	dustmasker -infmt fasta -in nt/nt-microbial.fna -level 20 -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g' > nt/library/nt-microbial-dustmasked.fna.tmp
	mv -v nt/library/nt-microbial-dustmasked.fna.tmp nt/library/nt-microbial-dustmasked.fna

nt/database.kdb: $(DB_F)
	krakenhll-build --db nt --threads $(THREADS)
