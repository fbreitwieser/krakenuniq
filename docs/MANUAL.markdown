Introduction
============

[Kraken] is a taxonomic sequence classifier that assigns taxonomic
labels to short DNA reads. It does this by examining the $k$-mers
within a read and querying a database with those $k$-mers. This database
contains a mapping of every $k$-mer in [Kraken]'s genomic library to the
lowest common ancestor (LCA) in a taxonomic tree of all genomes that
contain that $k$-mer. The set of LCA taxa that correspond to the $k$-mers
in a read are then analyzed to create a single taxonomic label for the
read; this label can be any of the nodes in the taxonomic tree.
[Kraken] is designed to be rapid, sensitive, and highly precise. Our
tests on various real and simulated data have shown [Kraken] to have
sensitivity slightly lower than Megablast with precision being slightly
higher. On a set of simulated 100 bp reads, [Kraken] processed over 1.3
million reads per minute on a single core in normal operation, and over
4.1 million reads per minute in quick operation.

The latest released version of Kraken will be available at the
[Kraken website], and the latest updates to the Kraken source code
are available at the [Kraken GitHub repository].

If you use [Kraken] in your research, please cite the [Kraken paper].
Thank you!

[Kraken]:                     http://ccb.jhu.edu/software/kraken/
[Kraken website]:             http://ccb.jhu.edu/software/kraken/
[Kraken paper]:               http://genomebiology.com/2014/15/3/R46
[Kraken GitHub repository]:   https://github.com/DerrickWood/kraken


System Requirements
===================

Note: Users concerned about the disk or memory requirements should
read the paragraph about MiniKraken, below.

* **Disk space**: Construction of Kraken's standard database will require at
    least 160 GB of disk space. Customized databases may require
    more or less space.  Disk space used is linearly proportional
    to the number of distinct $k$-mers; as of Feb. 2015, Kraken's
    default database contains just under 6 billion (6e9) distinct $k$-mers.

    In addition, the disk used to store the database should be
    locally-attached storage. Storing the database on a network
    filesystem (NFS) partition can cause Kraken's operation to be
    very slow, or to be stopped completely. As NFS accesses are
    much slower than local disk accesses, both preloading and database
    building will be slowed by use of NFS.

* **Memory**: To run efficiently, Kraken requires enough free memory to
    hold the database in RAM. While this can be accomplished using a
    ramdisk, Kraken supplies a utility for loading the database into
    RAM via the OS cache. The default database size is 75 GB (as of
    Feb. 2015), and so you will need at least that much RAM if you want
    to build or run with the default database.

* **Dependencies**: Kraken currently makes extensive use of Linux utilities
    such as sed, find, and wget. Many scripts are written using the
    Bash shell, and the main scripts are written using Perl. Core
    programs needed to build the database and run the classifier are
    written in C++, and need to be compiled using g++.  Multithreading
    is handled using OpenMP.  Downloads of NCBI data are performed by
    wget and in some cases, by rsync.  Most Linux systems that have any
    sort of development package installed will have all of the above
    listed programs and libraries available.

    Finally, if you want to build your own database, you will need to
    install the [Jellyfish] $k$-mer counter.  Note that Kraken only
    supports use of Jellyfish version 1.  Jellyfish version 2 is not
    yet compatible with Kraken.

* **Network connectivity**: Kraken's standard database build and download
    commands expect unfettered FTP and rsync access to the NCBI FTP
    server. If you're working behind a proxy, you may need to set
    certain environment variables (such as `ftp_proxy` or `RSYNC_PROXY`)
    in order to get these commands to work properly.

* **MiniKraken**: To allow users with low-memory computing environments to
    use Kraken, we supply a reduced standard database that can be
    downloaded from the Kraken web site. When Kraken is run with a
    reduced database, we call it MiniKraken.

    The database we make available is only 4 GB in size, and should
    run well on computers with as little as 8 GB of RAM. Disk space
    required for this database is also only 4 GB.

[Jellyfish]:  http://www.cbcb.umd.edu/software/jellyfish/


Installation
============

To begin using Kraken, you will first need to install it, and then
either download or create a database.

Kraken consists of two main scripts ("`kraken`" and "`kraken-build`"),
along with several programs and smaller scripts.  As part of the
installation process, all scripts and programs are installed in
the same directory.  After installation, you can move the main
scripts elsewhere, but moving the other scripts and programs
requires editing the scripts and changing the "`$KRAKEN_DIR`" variables.

Once a directory is selected, you need to run the following
command in the directory where you extracted the Kraken
source:

    ./install_kraken.sh $KRAKEN_DIR

(Replace "`$KRAKEN_DIR`" above with the directory where you want to
install Kraken's programs/directories.)

The `install_kraken.sh` script should compile all of Kraken's code
and setup your Kraken data directory.  Installation is successful
if you see the message "`Kraken installation complete.`"

Once installation is complete, you may want to copy the two main
Kraken scripts into a directory found in your `PATH` variable 
(e.g., "`$HOME/bin`"):

    cp $KRAKEN_DIR/bin/kraken $HOME/bin
    cp $KRAKEN_DIR/bin/kraken-build $HOME/bin

After installation, you're ready to either create or download a
database.


Kraken Databases
================

A Kraken database is a directory containing at least 4 files:

* `database.kdb`: Contains the $k$-mer to taxon mappings
* `database.idx`: Contains minimizer offset locations in database.kdb
* `taxonomy/nodes.dmp`: Taxonomy tree structure + ranks
* `taxonomy/names.dmp`: Taxonomy names

Other files may be present as part of the database build process.

In interacting with Kraken, you should not have to directly reference
any of these files, but rather simply provide the name of the directory
in which they are stored.  Kraken allows both the use of a standard
database as well as custom databases; these are described in the sections
[Standard Kraken Database] and [Custom Databases] below, respectively.


Standard Kraken Database
========================

To create the standard Kraken database, you can use the following command:

    kraken-build --standard --db $DBNAME

(Replace "`$DBNAME`" above with your preferred database name/location.
Please note that the database will use approximately 160 GB of
disk space during creation.)

This will download NCBI taxonomic information, as well as the
complete genomes in RefSeq for the bacterial, archaeal, and
viral domains.  After downloading all this data, the build
process begins; this is the most time-consuming step.  If you
have multiple processing cores, you can run this process with
multiple threads, e.g.:

    kraken-build --standard --threads 16 --db $DBNAME

Using 16 threads on a computer with 122 GB of RAM, the build
process took approximately an hour and a half (steps with an asterisk
have some multi-threading enabled) in February 2015:

     7m48s  *Step 1 (create set)
       n/a   Step 2 (reduce database, optional and skipped)
    53m16s  *Step 3 (sort set)
     1m04s   Step 4 (GI number to sequence ID map)
     0m27s   Step 5 (Sequence ID to taxon map)
    29m20s  *Step 6 (set LCA values)
    ------
    91m55s   Total build time

Note that if any step (including the initial downloads) fails,
the build process will abort.  However, `kraken-build` will
produce checkpoints throughout the installation process, and
will restart the build at the last incomplete step if you
attempt to run the same command again on a partially-built
database.

To create a custom database, or to use a database from another
source, see [Custom Databases].

Notes for users with lower amounts of RAM:

1) If you encounter problems with Jellyfish not being able
to allocate enough memory on your system to run the build
process, you can supply a smaller hash size to Jellyfish
using `kraken-build`'s `--jellyfish-hash-size` switch.  Each space
in the hash table uses approximately 6.9 bytes, so using
"`--jellyfish-hash-size 6400M`" will use a hash table size of
6.4 billion spaces and require 44.3 GB of RAM.

2) Kraken's build process will normally attempt to minimize
disk writing by allocating large blocks of RAM and operating
within them until data needs to be written to disk.  However,
this extra RAM usage may exceed your capacity.  In such cases,
you may want to use `kraken-build`'s `--work-on-disk` switch.  This
will minimize the amount of RAM usage and cause Kraken's build
programs to perform most operations off of disk files.  This
switch can also be useful for people building on a ramdisk or
solid state drive.  Please note that working off of disk files
can be quite slow on some computers, causing builds to take
several days if not weeks.


Classification
==============

To classify a set of sequences (reads), use the `kraken` command:

    kraken --db $DBNAME seqs.fa

Output will be sent to standard output by default.  The files
containing the sequences to be classified should be specified
on the command line.  Sequences can also be provided through
standard input using the special filename `/dev/fd/0`.

Note that to obtain optimum speeds, Kraken's database should be
loaded into RAM first.  This can be done through use of a ramdisk,
if you have superuser permissions.  Failing that, you can use
the `--preload` switch to `kraken`, e.g.:

    kraken --preload --db $DBNAME seqs.fa

The database files will be loaded before classification using this
switch.  See [Memory Usage and Efficiency] for more information.

The `kraken` program allows several different options:

* **Multithreading**: Use the `--threads NUM` switch to use multiple
    threads.

* **Quick operation**: Rather than searching all $k$-mers in a sequence,
    stop classification after the first database hit; use `--quick`
    to enable this mode.  Note that `--min-hits` will allow you to
    require multiple hits before declaring a sequence classified,
    which can be especially useful with custom databases when testing
    to see if sequences either do or do not belong to a particular
    genome.

* **Sequence filtering**: Classified or unclassified sequences can be
    sent to a file for later processing, using the `--classified-out`
    and `--unclassified-out` switches, respectively.

* **Output redirection**: Output can be directed using standard shell
    redirection (`|` or `>`), or using the `--output` switch.

* **FASTQ input**: Input is normally expected to be in FASTA format, but
    you can classify FASTQ data using the `--fastq-input` switch.

* **Compressed input**: Kraken can handle gzip and bzip2 compressed
    files as input by specifying the proper switch of `--gzip-compressed`
    or `--bzip2-compressed`.

* **Input format auto-detection**: If regular files are specified on
    the command line as input, Kraken will attempt to determine the
    format of your input prior to classification.  You can disable this
    by explicitly specifying `--fasta-input`, `--fastq-input`,
    `--gzip-compressed`, and/or `--bzip2-compressed` as appropriate.
    Note that use of the character device file `/dev/fd/0` to read
    from standard input (aka `stdin`) will **not** allow auto-detection.

* **Paired reads**: Kraken does not query $k$-mers containing ambiguous
    nucleotides (non-ACGT).  If you have paired reads, you can use this
    fact to your advantage and increase Kraken's accuracy by concatenating
    the pairs together with a single `N` between the sequences.  Using the
    `--paired` option when running `kraken` will automatically do this for
    you; simply specify the two mate pair files on the command line.  We
    have found this to raise sensitivity by about 3 percentage points over
    classifying the sequences as single-end reads.

To get a full list of options, use `kraken --help`.


Output Format
=============

Each sequence classified by Kraken results in a single line of
output.  Output lines contain five tab-delimited fields; from
left to right, they are:

1) "C"/"U": one letter code indicating that the sequence was
   either classified or unclassified.
2) The sequence ID, obtained from the FASTA/FASTQ header.
3) The taxonomy ID Kraken used to label the sequence; this is
   0 if the sequence is unclassified.
4) The length of the sequence in bp.
5) A space-delimited list indicating the LCA mapping of each $k$-mer
   in the sequence.  For example, "562:13 561:4 A:31 0:1 562:3"
   would indicate that:
     - the first 13 $k$-mers mapped to taxonomy ID #562
     - the next 4 $k$-mers mapped to taxonomy ID #561
     - the next 31 $k$-mers contained an ambiguous nucleotide
     - the next $k$-mer was not in the database
     - the last 3 $k$-mers mapped to taxonomy ID #562

For users who want the full taxonomic name associated with each input
sequence, we provide a script named `kraken-translate` that produces two
different output formats for classified sequences.  The script operates
on the output of `kraken`, like so:

    kraken --db $DBNAME sequences.fa > sequences.kraken
    kraken-translate --db $DBNAME sequences.kraken > sequences.labels

(The same database used to run `kraken` should be used to translate the
output; see [Kraken Environment Variables] below for ways to reduce
redundancy on the command line.)

The file `sequences.labels` generated by the above example is a text file
with two tab-delimited columns, and one line for each classified sequence
in `sequences.fa`; unclassified sequences are not reported by
`kraken-translate`.  The first column of `kraken-translate`'s output are the
sequence IDs of the classified sequences, and the second column contains
the taxonomy of the sequence.  For example, an output line from `kraken` of:

    C     SEQ1    562     36      562:6

Would result in a corresponding output line from `kraken-translate` of:

    SEQ1  root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia;Escherichia coli

Alternatively, `kraken-translate` accepts the option `--mpa-format` which
will report only levels of the taxonomy with standard rank assignments
(superkingdom, kingdom, phylum, class, order, family, genus, species),
and uses pipes to delimit the various levels of the taxonomy.  For example,
`kraken-translate --mpa-format --db $DBNAME` with the above example output
from `kraken` would result in the following line of output:

    SEQ1  d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli

Taxonomy assignments above the superkingdom (`d__`) rank are represented as
just "root" when using the `--mpa-report` option with `kraken-translate`.


Custom Databases
================

We realize the standard database may not suit everyone's needs.  Kraken
also allows creation of customized databases.

To build a custom database:

1) Install a taxonomy.  Usually, you will just use the NCBI taxonomy,
    which you can easily download using:
    
        kraken-build --download-taxonomy --db $DBNAME

    This will download the GI number to taxon map, as well as the
    taxonomic name and tree information from NCBI.  These files can
    be found in `$DBNAME/taxonomy/` .  If you need to modify the taxonomy,
    edits can be made to the `names.dmp` and `nodes.dmp` files in this directory;
    the `gi_taxid_nucl.dmp` file will also need to be updated appropriately.

2) Install a genomic library.  Four sets of standard genomes are
    made easily available through `kraken-build`:

    - bacteria: RefSeq complete bacterial/archaeal genomes
    - plasmids: RefSeq plasmid sequences
    - viruses: RefSeq complete viral genomes
    - human: GRCh38 human genome

    To download and install any one of these, use the `--download-library`
    switch, e.g.:

        kraken-build --download-library bacteria --db $DBNAME

    Other genomes can also be added, but such genomes must meet certain
    requirements:
    - Sequences must be in a FASTA file (multi-FASTA is allowed)
    - Each sequence's ID (the string between the `>` and the first
      whitespace character on the header line) must contain either
      a GI number to allow Kraken to lookup the correct taxa, or an
      explicit assignment of the taxonomy ID using `kraken:taxid` (see below).

    Replicons not downloaded from NCBI may need their taxonomy information
    assigned explicitly.  This can be done using the string `kraken:taxid|XXX`
    in the sequence ID, with `XXX` replaced by the desired taxon ID.  For
    example, to put a known adapter sequence in taxon 32630 ("synthetic
    construct"), you could use the following:

        >sequence16|kraken:taxid|32630  Adapter sequence
        CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA

    The `kraken:taxid` string must begin the sequence ID or be immediately
    preceded by a pipe character (`|`).  Explicit assignment of taxonomy IDs
    in this manner will override the GI number mapping provided by NCBI.

    If your genomes meet the requirements above, then you can add each
    replicon to your database's genomic library using the `--add-to-library`
    switch, e.g.:

        kraken-build --add-to-library chr1.fa --db $DBNAME
        kraken-build --add-to-library chr2.fa --db $DBNAME

    Note that if you have a list of files to add, you can do something like
    this in `bash`:

        for file in chr*.fa
        do
            kraken-build --add-to-library $file --db $DBNAME
        done

    Or even add all `*.fa` files found in the directory `genomes`:

        find genomes/ -name '*.fa' -print0 | \
            xargs -0 -I{} -n1 kraken-build --add-to-library {} --db $DBNAME

    (You may also find the `-P` option to `xargs` useful to add many files in
    parallel if you have multiple processors.)

3) Once your library is finalized, you need to build the database.
    Depending on your size requirements, you may want to adjust the
    $k$-mer and/or minimizer lengths from the defaults.  Except for some
    small bookkeeping fields, a Kraken database will use 
    $sD$ + $8(4^{M})$
    bytes, where $s$ is the number of bytes used to store the $k$-mer/taxon 
    pair (usually 12, but lower for smaller $k$-mers), $D$ is the number of
    distinct $k$-mers in your library and
    $M$ is the length (in bp) of the minimizers.  Although $D$ does increase
    as $k$ increases, it is impossible to know exactly how many distinct
    $k$-mers will exist in a library for a given $k$ without actually
    performing the count.  By default, $k$ = 31 and $M$ = 15.

    The minimizers serve to keep $k$-mers that are adjacent in query
    sequences close to each other in the database, which allows
    Kraken to exploit the CPU cache.  Changing the value of $M$ can
    significantly affect the speed of Kraken, and neither increasing
    or decreasing $M$ will guarantee faster or slower speed.

    To build the database, you'll use the `--build` switch:

        kraken-build --build --db $DBNAME

    As noted above, you may want to also use any of `--threads`,
    `--kmer-len`, or `--minimizer-len` to adjust the database build
    time and/or final size.

4) Shrinking the database: The "--shrink" task allows you to take
    an existing Kraken database and create a smaller MiniKraken database
    from it.  The use of this option removes all but a specified number of
    $k$-mer/taxon pairs to create a new, smaller database.  For example:

        kraken-build --shrink 10000 --db $DBNAME --new-db minikraken

    This will create a new database named `minikraken` that contains
    10000 $k$-mers selected from across the original database (`$DBNAME`).

    The `--shrink` task is only meant to be run on a completed database.
    However, if you know before you create a database that you will
    only be able to use a certain amount of memory, you can use the
    `--max-db-size` switch for the `--build` task to provide a maximum
    size (in GB) for the database.  This allows you to create a MiniKraken
    database without having to create a full Kraken database first.

A full list of options for `kraken-build` can be obtained using
`kraken-build --help`.

After building a database, if you want to reduce the disk usage of
the database you can use `kraken-build`'s `--clean` switch to remove
all intermediate files from the database directory.

Memory Usage and Efficiency
===========================

Kraken's execution requires many random accesses to a very large file.
To obtain maximal speed, these accesses need to be made as quickly as
possible.  This means that the database must be in physical memory
during execution.  Although we provide the `--preload` option to Kraken
for users who cannot use a ramdisk, the ramdisk is likely the simplest
option, and is well-suited for installations on computers where Kraken
is to be run a majority of the time.  In addition, using a ramdisk
allows the initial start-up of Kraken to be accomplished much more quickly.
If a ramdisk is used, the `--preload` switch should not be used.

We also note that in some cases, `--preload` may not be needed (or even
advisable).  If you know that your database is already in memory (for
example, if it has been recently read or unzipped, then it should be in
your operating system cache, which resides in physical memory), then there
is no need to perform this step.  We have noticed that in low-memory (~8 GB)
situations, preloading a MiniKraken DB is actually much slower than simply
using `cat minikraken/database.* > /dev/null`.  The selection of the best way
to get the database into memory is dependent on several factors, including
your total amount of RAM, operating system, and current free memory.  For this
reason, you may need to experiment with your own setup to find a good solution
for you.

To create a ramdisk, you will need to have superuser (root) permission.
As root, you can use the following commands to create a ramdisk:

    mkdir /ramdisk
    mount -t ramfs none /ramdisk

Optionally, you may have a trusted user who you want to be able to copy
databases into this directory.  In that case, you'll need to make that
user the owner of the directory via chown.

To put the database on the ramdisk, simply copy the database directory
to the ramdisk directory:

    cp -a $DBNAME /ramdisk

And then you can use it with Kraken by specifying the database copy on
the ramdisk, e.g.:

    kraken --db /ramdisk/$DBNAME seqs.fa

Note that anything copied into a ramdisk will be deleted if the ramdisk
is unmounted or the computer is restarted, so make sure that you have a
copy of the database on a hard disk (or other non-volatile storage).


Note that when using the `--paired` option, Kraken will not (by default)
make any attempt to ensure that the two files you specify are indeed
matching sets of paired-end reads.  To verify that the names of each
read do indeed match, you can use the `--check-names` option in
combination with the `--paired` option.


Sample Reports
==============

To get an idea as to Kraken's results across an entire sample, we provide
the `kraken-report` script.  It is used like this:

    kraken-report --db $DBNAME kraken.output

Note that the database used must be the same as the one used to generate
the output file, or the report script may encounter problems.  Output is
sent to standard output.

The output of `kraken-report` is tab-delimited, with one line per taxon.
The fields of the output, from left-to-right, are as follows:

1) Percentage of reads covered by the clade rooted at this taxon
2) Number of reads covered by the clade rooted at this taxon
3) Number of reads assigned directly to this taxon
4) A rank code, indicating (U)nclassified, (D)omain, (K)ingdom,
    (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
    All other ranks are simply '-'.
5) NCBI taxonomy ID
6) indented scientific name 

The scientific names are indented using spaces, according to the tree
structure specified by the taxonomy.

By default, taxa with no reads assigned to (or under) them will not have
any output produced.  However, if you wish to have all taxa displayed,
you can use the `--show-zeros` switch to do so.  This can be useful if
you are looking to do further downstream analysis of the reports, and
want to compare samples.  Sorting by the taxonomy ID (using `sort -nf5`)
can provide a consistent line ordering between reports.

In addition, we also provide the program `kraken-mpa-report`; this program
provides output in a format similar to MetaPhlAn's tab-delimited output.
For `kraken-mpa-report`, multiple Kraken output files can be specified on
the command line and each will be treated as a separate sample.  For each
taxon at the standard ranks (from domain to species), the count of reads
in each sample assigned to any node in the clade rooted at that taxon is
displayed.  `kraken-mpa-report` is run in the same manner as `kraken-report`,
and its output is also sent to standard output.


Confidence Scoring
==================

At present, we have not yet developed a confidence score with a solid
probabilistic interpretation for Kraken.  However, we have developed a
simple scoring scheme that has yielded good results for us, and we've
made that available in the `kraken-filter` script.  The approach we use
allows a user to specify a threshold score in the [0,1] interval; the
`kraken-filter` script then will adjust labels up the tree until the
label's score (described below) meets or exceeds that threshold.  If
a label at the root of the taxonomic tree would not have a score exceeding
the threshold, the sequence is called unclassified by kraken-filter.

A sequence label's score is a fraction $C$/$Q$, where $C$ is the number of
$k$-mers mapped to LCA values in the clade rooted at the label, and $Q$ is the
number of $k$-mers in the sequence that lack an ambiguous nucleotide (i.e.,
they were queried against the database).  Consider the example of the
LCA mappings in Kraken's output given earlier:

"562:13 561:4 A:31 0:1 562:3" would indicate that:

* the first 13 $k$-mers mapped to taxonomy ID #562
* the next 4 $k$-mers mapped to taxonomy ID #561
* the next 31 $k$-mers contained an ambiguous nucleotide
* the next $k$-mer was not in the database
* the last 3 $k$-mers mapped to taxonomy ID #562

In this case, ID #561 is the parent node of #562.  Here, a label of #562
for this sequence would have a score of $C$/$Q$ = (13+3)/(13+4+1+3) = 16/21.
A label of #561 would have a score of $C$/$Q$ = (13+4+3)/(13+4+1+3) = 20/21.
If a user specified a threshold over 16/21, kraken-filter would adjust the
original label from #562 to #561; if the threshold was greater than 20/21,
the sequence would become unclassified.

`kraken-filter` is used like this:

    kraken-filter --db $DBNAME [--threshold NUM] kraken.output

If not specified, the threshold will be 0.  `kraken-filter`'s output is
similar to `kraken`'s, but a new field between the length and LCA mapping
list is present, indicating the new label's score (or the root label's
score if the sequence has become unclassified).

To give some guidance toward selecting an appropriate threshold, we
show here the results of different thresholds on the MiSeq metagenome
from the [Kraken paper] \(see the paper for more details; note that the
database used here is more recent than that used in the paper).
Precision, sensitivity, and F-score are measured at the genus rank:

<div id="confidence-score-table">

  Thres   Prec     Sens     F-score
  -----   -----    -----    -------
  0       95.43    77.32    85.43
  0.05    97.28    76.31    85.53
  0.10    98.25    75.13    85.15
  0.15    98.81    73.87    84.54
  0.20    99.13    72.82    83.96
  0.25    99.38    71.74    83.33
  0.30    99.55    70.75    82.71
  0.35    99.61    69.53    81.90
  0.40    99.66    68.35    81.09
  0.45    99.70    66.93    80.09
  0.50    99.71    65.49    79.06

</div>

As can be seen, with no threshold (i.e., Kraken's original labels),
Kraken's precision is fairly high, but it does increase with the
threshold.  Diminishing returns apply, however, and there is a loss
in sensitivity that must be taken into account when deciding on the
threshold to use for your own project.


Kraken Environment Variables
============================

The Kraken programs (with the exception of `kraken-build`) support the
use of some environment variables to help in reducing command line
lengths:

* **`KRAKEN_NUM_THREADS`**: this variable is only used by `kraken`; if the
    `--threads` option is not supplied to `kraken`, then the value of this
    variable (if it is set) will be used as the number of threads to run
    `kraken`.

* **`KRAKEN_DB_PATH`**: much like the `PATH` variable is used for executables
    by your shell, `KRAKEN_DB_PATH` is a colon-separated list of directories
    that will be searched for the database you name if the named database
    does not have a slash (`/`) character.  By default, Kraken assumes the
    value of this variable is "`.`" (i.e., the current working directory).
    This variable can be used to create one (or more) central repositories
    of Kraken databases in a multi-user system.  Example usage in bash:

        export KRAKEN_DB_PATH="/home/user/my_kraken_dbs:/data/kraken_dbs:"

    This will cause three directories to be searched, in this order:

    1) `/home/user/my_kraken_dbs`
    2) `/data/kraken_dbs`
    3) the current working directory (caused by the empty string as
         the third colon-separated field in the `KRAKEN_DB_PATH` string)

    The search for a database will stop when a name match is found; if
    two directories in the `KRAKEN_DB_PATH` have databases with the same
    name, the directory of the two that is searched first will have its
    database selected.

    If the above variable and value are used, and the databases
    `/data/kraken_dbs/mainDB` and `./mainDB` are present, then

        kraken --db mainDB sequences.fa

    will classify `sequences.fa` using `/data/kraken_dbs/mainDB`; if instead
    you wanted to use the `mainDB` present in the current directory,
    you would need to specify a directory path to that database in order
    to circumvent searching, e.g.:

        kraken --db ./mainDB sequences.fa

    Note that the `KRAKEN_DB_PATH` directory list can be skipped by the use
    of any absolute (beginning with `/`) or relative pathname (including
    at least one `/`) as the database name.

* **`KRAKEN_DEFAULT_DB`**: if no database is supplied with the `--db` option,
    the database named in this variable will be used instead.  Using this
    variable, you can avoid using `--db` if you only have a single database
    that you usually use, e.g. in bash:

        export KRAKEN_DEFAULT_DB="/home/user/krakendb"
        kraken sequences.fa | kraken-report > sequences.kreport

    This will classify `sequences.fa` using the `/home/user/krakendb` directory.

    Note that the value of `KRAKEN_DEFAULT_DB` will also be interpreted in
    the context of the value of `KRAKEN_DB_PATH` if you don't set
    `KRAKEN_DEFAULT_DB` to an absolute or relative pathname.  Given the earlier
    example in this section, the following:

        export KRAKEN_DEFAULT_DB="mainDB"
        kraken sequences.fa

    will use `/data/kraken_dbs/mainDB` to classify `sequences.fa`.


Upgrading Databases to v0.10+
=============================

The minimizer ordering in Kraken versions prior to v0.10.0-beta was a
simple lexicographical ordering that provided a suboptimal distribution
of k-mers within the bins.  Ideally, the bin sizes would be uniform,
but simple lexicographical ordering creates a bias toward low-complexity
minimizers.  To resolve this, the ordering is now "scrambled" by XORing all
minimizers with a predefined constant to toggle half of each minimizer's
bits before sorting.  The more evenly distributed bins provide better
caching performance, but databases created in this way are not compatible
with earlier versions of Kraken.  Kraken versions from v0.10.0-beta up to
(but not including) v1.0 will support the use of the older databases, but
we nonetheless recommend one of the two following options:

1) Build a new database.  This is the preferred option, as a newly-created
    database will have the latest genomes and NCBI taxonomy information.

2) Re-sort an existing database.  If you have a custom database, you may
    want to simply reformat the database to provide you with Kraken's
    increased speed.  To do so, you'll need to do the following:

        kraken-build --upgrade --db $DBNAME

    (**Note**: the `--threads` switch is both valid and encouraged with this
    operation.)

    This command will **not** delete your existing `$DBNAME/database.*`
    files, but will simply rename them.  If you're satisfied with the new
    database's performance, then you can use `kraken-build`'s `--clean`
    option to remove the old files and save space.

    Sorting the database is step 3 of the build process, so you should
    expect a database upgrade to take about as long as step 3 took when
    building the original database.

Note that the rest of Kraken v0.10.0-beta's speed improvements are available
without upgrading or changing your database.
