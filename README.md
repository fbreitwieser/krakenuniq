KrakenUniq: confident and fast metagenomics classification using unique k-mer counts
===============================================

False-positive identifications are a significant problem in metagenomics classification. KrakenUniq (formerly KrakenHLL) is a novel metagenomics classifier that combines the fast k-mer-based classification of [Kraken](https://github.com/DerrickWood/kraken) with an efficient algorithm for assessing the coverage of unique k-mers found in each species in a dataset. On various test datasets, KrakenUniq gives better recall and precision than other methods and effectively classifies and distinguishes pathogens with low abundance from false positives in infectious disease samples. By using the probabilistic cardinality estimator HyperLogLog, KrakenUniq runs as fast as Kraken and requires little additional memory. NEW in v0.7.0: KrakenUniq can run on standard laptops and desktops with as little as 16GB of RAM using the --preload-size option (see below).

**If you use KrakenUniq in your research, please cite our publication:** [KrakenUniq: confident and fast metagenomics classification using unique k-mer counts. Breitwieser FP, Baker DN, Salzberg SL. Genome Biology, Dec 2018. https://doi.org/10.1186/s13059-018-1568-0](https://doi.org/10.1186/s13059-018-1568-0)
## KrakenUniq databases available for direct download
We now have two standard Kraken1/KrakenUniq databases available for free download from the Amazon cloud. You can find links at [https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2). One is our "standard" database with all RefSeq bacteria, archaea, and viruses, plus common vectors and the human genome. The other is all of the first database plus all available genomes of eukaryotic pathogens. Each DB is over 300GB, and by downloading them you can avoid having to build them yourself.

# Announcements

## New release v1.0.4
In this release we removed the requirement to have "file" command in the Docker or Singularity image (thanks @lskatz and @boulund).  We now force --preload switch when building the database for speed.

## New release v1.0.3
In this release we added a Dockerfile (thanks @Jessime).  There are also few minor updates to documentation. This code of this release has been reviewed in connection with our publication in JOSS.

## New release v1.0.2
This release fixes the issue with possibly incorrect output produced when running multiple krakenuniq processes in the same folder with --paired input files.  

## New release v1.0.1
This release fixes the Issue #116 and #117. Thanks to @boulund for the fix!

## New release v1.0.0
This is the official 1.0.0 release of KrakenUniq. This release fixes the bug with downloading databases (Issue #87). Thanks to @clescoat for the fix!

## New release v0.7.3
This maintenance release provides the following updates: 
(1) fixes issues with building large databases
(2) installs Jellyfish version 1 under KRAKENUNIQ_INSTALL_DIR/jellyfish-install/bin/, if -j switch is used (KrakenUniq requires Jellyfish version 1 to build databases)
(3) fixes --work-on-disk option (#97)

## New release v0.7.2
This maintenance release fixes the --paired option in krakenuniq and the minor problem at the last stage of building a new database (report).

## New release v0.7.1
This minor release fixes a bug in the Makefile that resulted in installation of unusable executables count_unique and set_lcas. The bug resulted in fatal error in building a new krakenuniq database. Classification with an existing database was not affected.

## New Release v0.7
New option for low-memory computers: --preload-size.

By default, KrakenUniq performs memory mapping to load the database; i.e., it does not load the entire database into main memory. (Kraken 1 employs the same strategy.) This makes classification of larger read datasets much slower, but it allows KrakenUniq to run on machines with low available main memory. If enough free RAM is available to hold the entire database in main memory, users are recommended to explicitly load the entire database prior to classification using the flag --preload, which dramatically speeds up the classification, often by a factor of 20 or more.

To improve the performance when not enough main memory is available to load the entire database into RAM, we added a new capability to KrakenUniq. When using this new feature, only a chunk of the database is loaded into memory at a time, after which the algorithm iterates over the reads and looks up all k-mers in those reads that are matching in this database chunk. This process is repeated until the entire database has been processed. The k-mer lookups are then merged, and reads are classified based on the results of the full database. This new feature makes it feasible to run KrakenUniq on very large datasets and huge databases on virtually any computer, even a laptop, while providing exact classifications that are identical to those of KrakenUniq in its other modes. Users can employ this feature with --preload-size and specify the amount of available main memory that they want to use for loading chunks of the database, e.g., --preload-size 8G or --preload-size 500M.

IMPORTANT!  The --preload-size option can only be used with a single input database.

This release also includes an improvement for automatic detection of input format.
The input format (fastq or fasta, bzip2 or gzip compressed) is now detected automatically. No need to use --fasta-input, --fastq-input, --gzip-compressed or --bzip2-compressed switches.

The improvements included in this release are described in the preprint posted on bioRxiv: https://www.biorxiv.org/content/10.1101/2022.06.01.494344v1

## New Release v0.6
This release fixes database preload option. Now --preload option will force loading the database in physical RAM (not swap) if enough physical RAM is available. KrakenUniq (and also Kraken) often ran very slow with really big databases. The problem was that --preload didn't truly force to load the DB in memory, so it spends forever (many days) going back and forth to disk. With the fix included in this release, krakenuniq ran in 16 minutes on a database where before it took >100 hours. 

## Installation
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/krakenuniq/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/krakenuniq/badges/latest_release_date.svg)](https://anaconda.org/bioconda/krakenuniq)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/krakenuniq/badges/platforms.svg)](https://anaconda.org/bioconda/krakenuniq)

KrakenUniq is available in the Anaconda cloud. To install, type:

```
conda install -c bioconda krakenuniq
```
This is the bioconda link for KrakenUniq:  https://anaconda.org/bioconda/krakenuniq

Installation from a release -- recommended.  Check https://github.com/fbreitwieser/krakenuniq/releases for the latest version.  Then run the following commands (\<VERSION\> is the version, e.g. 0.7.3):
```
wget https://github.com/fbreitwieser/krakenuniq/archive/refs/tags/v<VERSION>.tar.gz
tar xzf v<VERSION>.tar.gz
cd krakenuniq-<VERSION>
./install_krakenuniq /PATH/TO/INSTALL_DIR
```

Installation with Docker:
```
docker build -t krakenuniq . 
docker run --rm -it -it krakenuniq krakenuniq --help
```

Installation from source from GitHub (the latest -- may not be stable):
```
git clone https://github.com/fbreitwieser/krakenuniq
cd krakenuniq
./install_krakenuniq /PATH/TO/INSTALL_DIR
```

KrakenUniq requires Jellyfish version 1.x.x to be installed for the database building step (`krakenuniq-build`). Starting with v0.7.3, krakenuniq will download and install jellyfish v1.1.12 automatically by default. To skip installing Jellyfish use the `-s` flag for the `install_krakenuniq.sh` script. Alternatively, you can specify the Jellyfish path to `krakenuniq-build` with `krakenuniq-build --jellyfish-bin /usr/bin/jellyfish1`.

OSX by default links `g++` to `clang` without OpenMP support. When using clang, you may get the error `clang: fatal error: unsupported option '-fopenmp'`. To fix this, install `g++` with HomeBrew and use the `-c` option of `krakenuniq_install.sh` to specify the HomeBrew version of `g++`, which is accessible with `g++-8`: 
``` 
brew install gcc
./install_krakenuniq -c g++-8 /PATH/TO/INSTALL_DIR
```

## Command-line options

This is the output of `krakenuniq --help`:
```
Usage: $PROG --report-file FILENAME [options] <filename(s)>
Options:
  --db NAME               Name for Kraken DB (default: none)
  --threads NUM           Number of threads (default: 1)
  --hll-precision INT     Precision for HyperLogLog k-mer cardinality estimation, between 10 and 18 (default: 12)
  --exact                 Compute exact cardinality instead of estimate (slower, requires memory proportional to cardinality!)
  --quick                 Quick operation (use first hit or hits)
  --min-hits NUM          In quick op., number of hits req'd for classification
                          NOTE: this is ignored if --quick is not specified
  --unclassified-out FILENAME
                          Print unclassified sequences to filename
  --classified-out FILENAME
                          Print classified sequences to filename
  --output FILENAME       Print output to filename (default: stdout); "off" will
                          suppress normal output
  --only-classified-output
                          Print no Kraken output for unclassified sequences
  --preload               Loads the entire DB into memory before classification
  --preload-size SIZE     Loads DB into memory in chunks of SIZE, e.g. 500M or 7G (if RAM is small), overrides --preload flag
  --paired                The two filenames provided are paired-end reads
  --check-names           Ensure each pair of reads have names that agree
                          with each other; ignored if --paired is not specified
  --help                  Print this message
  --version               Print version information
Experimental:
  --uid-mapping           Map using UID database
The file format (fasta/fastq) and compression (gzip/bzip2) do not need to be specified anymore.
The format is detected automatically.
```

## Database building

Note that KrakenUniq natively supports Kraken 1 databases (however not Kraken 2). If you have existing Kraken databases, you may run KrakenUniq directly on them, though for support of taxon nodes for genomes and sequences (see below) you will need to rebuild them with KrakenUniq. For building a custom database, there are three requirements:

1. Sequence files (FASTA format)
2. Mapping files (tab separated format, `sequence header<tab>taxID`
3. NCBI taxonomy files (though a custom taoxnomies may be used, too)

While you may supply this information yourself, `krakenuniq-download` supports a variety of data sources to download the taxonomy, sequence and mapping files. Please find examples below on how to download different sequence sets:

```
## Download the taxonomy and the sequences
krakenuniq-download --db DBDIR <PATTERN>

<PATTERN> can be one of
     'contaminants'     Contaminant sequences from UniVec and EmVec.
     'taxonomy'         NCBI taxonomy mappings from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
     'nucleotide'       Download nucleotide sequences using a query specified using --search or --ac.
     'microbial-nt'     Download microbial sequences from nt database.
     'nt'               Download sequences from nt database, specified via --taxa.
     'viral-neighbors'  Download viral strain sequences from the NCBI Viral Genome Resource.
                        (Search: \"$vir_nbr_search_term\").
     'genbank/DOMAIN'   Download all complete genomes for DOMAIN from GenBank.
     'refseq/DOMAIN'    Download all complete genomes for DOMAIN from RefSeq.
     'refseq/DOMAIN/ASS_LEVEL'
     'refseq/DOMAIN/ASS_LEVEL/COLUMN=value1(/COLUMN=value2)*' 
        Possible values for DOMAIN: @ALL_GENOMES.
        Possible values for ASS_LEVEL: Any, Complete_Genome, Chromosome, Scaffold and Contig. 
        Possible values for COLUMN: Any column in the NCBI assembly_summary.txt, e.g. species_taxid or assembly_accession.
        
EXAMPLES:

## Download all human assemblies
krakenuniq-download --db DBDIR 'refseq/vertebrate_mammalian/Any/species_taxid=9606'

## Download all complete bacterial and archaeal genomes genomes in RefSeq using 10 threads, and masking low-complexity sequences in the genomes
krakenuniq-download --db DBDIR --threads 10 --dust refseq/bacteria refseq/archaea

## Download all contaminant sequences from UniVec and EmVec, plus the human reference genome
krakenuniq-download --db DBDIR refseq/vertebrate_mammalian/Chromosome/species_taxid=9606

## Download all viral genomes from RefSeq plus viral 'neighbors' in NCBI Nucleotide
krakenuniq-download --db DBDIR refseq/viral/Any viral-neighbors

## Download all microbial (including eukaryotes) sequences in the NCBI nt database
krakenuniq-download --db DBDIR --dust microbial-nt
```

To build the database indices on the downloaded files, run `krakenuniq-build --db DBDIR`.  To build a database with a *k*-mer length of 31 (the default), adding virtual taxonomy nodes for genomes and sequences (off by default), run `krakenuniq-build` with the following parameters:
```
krakenuniq-build --db DBDIR --kmer-len 31 --threads 10 --taxids-for-genomes --taxids-for-sequences

```

For more information on taxids for genomes and sequences, look at the [manual](MANUAL.md). The building step may take up to a couple of days on large sequence sets such as nt.

## Classification

To run classification on a pair of FASTQ files, use `krakenuniq`.

```
krakenuniq --db DBDIR --threads 10 --report-file REPORTFILE.tsv > READCLASSIFICATION.tsv
```

It can be advantegeous to preload the database prior to the first run. KrakenUniq uses mmap to map the database files into memory, which reads the file on demand. `krakenuniq --preload` reads the full database into memory, so that subsequent runs can benefit from the mapped pages. You do not need to specify preload before every run, but only after restarting the machine or when using a new database.

```
krakenuniq --db DBDIR --preload --threads 10
krakenuniq --db DBDIR --threads 10 --report-file REPORTFILE.tsv > READCLASSIFICATION.tsv
...
```


## FAQ

### Memory requirements

Stating with version 0.7.1, KrakenUniq can efficiently classify reads with databases exceeding the avalilable RAM using --preload-size switch.  10-16Gb of RAM is enough to classify with ~400Gb database in a reasonable amount of time.

### KrakenUniq vs Kraken vs Kraken 2

KrakenUniq was built on top of Kraken, and supports Kraken 1 databases natively. Kraken 2 is a new development that has a different database format, which is not supported by KrakenUniq.

### Differences to `kraken`
 - Use `krakenuniq --report-file FILENAME ...` to write the kraken report to `FILENAME`.
 - Use `krakenuniq --db DB1 --db DB2 --db DB3 ...` to first attempt, for each k-mer, to assign it based on DB1, then DB2, then DB3. You can use this to prefer identifications based on DB1 (e.g. human and contaminant sequences), then DB2 (e.g. completed bacterial genomes), then DB3, etc. Note that this option is incompatible with `krakenuniq-build --taxids-for-genomes --taxids-for-sequences` since the taxDB between the databases has to be absolutely the same.
 - Add a suffix `.gz` to output files to generate gzipped output files

### Differences to `kraken-build`
 - Use `krakenuniq-build --taxids-for-genomes --taxids-for-sequences ...` to add pseudo-taxonomy IDs for each sequence header and genome assembly (when using `krakenuniq-download`). 
 - `seqid2taxid.map` mapping sequence IDs to taxonomy IDs does NOT parse or require `>gi|`, but rather the sequence ID is the header up to just before the first space
 
### Building a microbial nt database

KrakenUniq supports building databases on subsets of the NCBI nucleotide collection nr/nt, which is most prominently the standard database for BLASTn. On the command line, you can specify to extract all bacterial, viral, archaeal, protozoan, fungal and helminth sequences. The list of protozoan taxa is based on [Kaiju's](https://raw.githubusercontent.com/bioinformatics-centre/kaiju/master/util/taxonlist.tsv).

Example command line:
```
krakenuniq-download --db DB --taxa "archaea,bacteria,viral,fungi,protozoa,helminths" --dust --exclude-environmental-taxa microbial-nt
```


### Custom databases with NCBI taxonomy
To build a custom database with the NCBI taxonomy, first download the taxonomy files with
```
krakenuniq-download --db DBDIR taxonomy
```
Then you can add the desired sequence files to the `DBDIR/library` directory:
```
cp SEQ1.fa SEQ2.fa DBDIR/library
```
KrakenUniq needs a _sequence ID to taxonomy ID mapping_ for each sequence. This mappings can be provided in the `DBDIR/library/*.map` - KrakenUniq pools all `.map` files inside of the `library/` folder prior to database building. Format: three tab-separated fields that are, in order, the sequence ID (i. e. the sequence header without '>' up to the first space), the taxonomy ID and the genome or assembly name:
```
Strain1_Chr1_Seq     <tab> 562 <tab> E. Coli Strain Foo
Strain1_Chr2_Seq     <tab> 562 <tab> E. Coli Strain Foo
Strain1_Plasmid1_Seq <tab> 562 <tab> E. Coli Strain Foo
Strain2_Chr1_Seq     <tab> 621 <tab> S. boydii Strain Bar
Strain2_Plasmid1_Seq <tab> 621 <tab> S. boydii Strain Bar
```
The third column is optional, and used by KrakenUniq only when `--taxids-for-genomes` is specified for `krakenuniq-build` to add new nodes in the taxonomy tree for the genome. If you'd like to have the sequences identifier in the taxonomy report, too, specifiy `--taxids-for-sequences` for `krakenuniq-build` as well.

Finally, run `krakenuniq-build`:
```
krakenuniq-build --db DBDIR --taxids-for-genomes --taxids-for-sequences
```

Note that for custom databases with fewer sequences you might want to choose a smaller k (default: `--kmer-len 31`) and minimizer length (default: `--minimizer-len 15`).

### Custom databases with custom taxonomies

When using custom taxonomies, please provide `DBDIR/taxonomy/nodes.dmp` and `DBDIR/taxonomy/names.dmp` according to the format of NCBI taxonomy dumps.

## License

The code adpated from Kraken 1 is licensed under GPL 3.0. All code added in this project (such as the
HyperLogLog algorithm code) is dual-licensed under MIT and GPL 3.0 (or any later version), unless stated otherwise.
You can choose between one of them if you use that work.
