KrakenHLL taxonomic sequence classification system with unique k-mer counting
===============================================

[Kraken](https://github.com/DerrickWood/kraken) is a fast taxonomic classifier for metagenomics data. This project, kraken-hll, adds some additional functionality - most notably a unique k-mer count using the HyperLogLog algorithm. Spurious identifications due to sequence contamination in the dataset or database often leads to many reads, however they usually cover only a small portion of the genome. 

KrakenHLL computes the number of unique k-mers observed for each taxon, which allows to filter more false positives.  Here's a small example of a classification against a viral database with k=25. There are three species identified by just one read - Enterobacteria phage BP-4795, Salmonella phage SEN22, Sulfolobus monocaudavirus SMV1. Out of those, the identification of Salmonella phage SEN22 is the strongest, as there read was matched with 116 k-mers that are unique to the sequence, while the match to Sulfolobus monocaudavirus SMV1 is only based on a single 25-mer.

```
99.0958 2192    2192    255510  272869  no rank 0   unclassified
0.904159    20  0   2361    2318    no rank 1   root
0.904159    20  0   2361    2318    superkingdom    10239     Viruses
0.904159    20  0   2361    2318    no rank 35237       dsDNA viruses, no RNA stage
0.768535    17  0   2074    2063    order   548681        Herpesvirales
0.768535    17  0   2074    2063    family  10292           Herpesviridae
0.768535    17  0   2074    2063    subfamily   10374             Gammaherpesvirinae
0.768535    17  0   2074    2063    genus   10375               Lymphocryptovirus
0.768535    17  16  2001    1987    species 10376                 Human gammaherpesvirus 4
0.045208    1   1   4   4   sequence    1000041143                  KC207814.1 Human herpesvirus 4 strain Mutu, complete genome
0.0904159   2   0   254 254 order   28883         Caudovirales
0.045208    1   0   28  28  family  10699           Siphoviridae
0.045208    1   0   28  28  genus   186765            Lambdavirus
0.045208    1   0   28  28  no rank 335795              unclassified Lambda-like viruses
0.045208    1   1   28  28  species 196242                Enterobacteria phage BP-4795
0.045208    1   0   116 116 family  10744           Podoviridae
0.045208    1   0   116 116 no rank 196895            unclassified Podoviridae
0.045208    1   0   116 116 no rank 1758253             Escherichia phage phi191 sensu lato
0.045208    1   1   116 116 species 1647458               Salmonella phage SEN22
0.045208    1   0   1   1   no rank 51368         unclassified dsDNA viruses
0.045208    1   1   1   1   species 1351702         Sulfolobus monocaudavirus SMV1
```

## Usage

For usage, see `krakenhll --help`. Note that you can use the same database as Kraken with one difference - instead of the files `DB_DIR/taxonomy/nodes.dmp` and `DB_DIR/taxonomy/names.dmp` than kraken relies upon, `kraken-hll` needs the file `DB_DIR/taxDB`. This can be generated with the script `build_taxdb`: `KRAKEN_DIR/build_taxdb DB_DIR/taxonomy/names.dmp DB_DIR/taxonomy/nodes.dmp > DB_DIR/taxDB`. The code behind the taxDB is based on [k-SLAM](https://github.com/aindj/k-SLAM).

### Differences to `kraken`
 - Use `krakenhll --report-file FILENAME ...` to write the kraken report to `FILENAME`.
 - Use `krakenhll --db DB1 --db DB2 --db DB3 ...` to first attempt, for each k-mer, to assign it based on DB1, then DB2, then DB3. You can use this to prefer identifications based on DB1 (e.g. human and contaminant sequences), then DB2 (e.g. completed bacterial genomes), then DB3, etc. Note that this option is incompatible with `krakenhll-build --generate-taxonomy-ids-for-sequences` since the taxDB between the databases has to be absolutely the same.
 - Add a suffix `.gz` to output files to generate gzipped output files

### Differences to `kraken-build`
 - Use `krakenhll-build --generate-taxonomy-ids-for-sequences ...` to add pseudo-taxonomy IDs for each sequence header. An example for the result using this is in the ouput above - one read has been assigned specifically to `KC207814.1 Human herpesvirus 4 strain Mutu, complete genome`.
 - `seqid2taxid.map` mapping sequence IDs to taxonomy IDs does NOT parse or require `>gi|`, but rather the sequence ID is the header up to just before the first space
 
 ## FAQ
 
 ### Installing KrakenHLL on MacOS
OSX by default links `g++` to `clang` without OpenMP support. You can install `g++` with HomeBrew and use the `-c` option of `krakenhll_install.sh` to specify the HomeBrew `g++`: 
``` 
brew install g++
./install_krakenhll -c g++-7
```

### Custom databases with NCBI taxonomy
To build a custom database with the NCBI taxonomy, first download the taxonomy files with
```
krakenhll-download --db DBDIR taxonomy
```
Then you can add the desired sequence files to the `DBDIR/library` directory:
```
cp SEQ1.fa SEQ2.fa DBDIR/library
```
KrakenHLL needs a _sequence ID to taxonomy ID mapping_ for each sequence. This mapping can be provided in the `DBDIR/library/seqid2taxid.map`. Format: three tab-separated fields that are, in order, the sequence ID (i. e. the sequence header without '>' up to the first space), the taxonomy ID and the genome or assembly name:
```
Strain1_Chr1_Seq     <tab> 562 <tab> E. Coli Strain Foo
Strain1_Chr2_Seq     <tab> 562 <tab> E. Coli Strain Foo
Strain1_Plasmid1_Seq <tab> 562 <tab> E. Coli Strain Foo
Strain2_Chr1_Seq     <tab> 621 <tab> S. boydii Strain Bar
Strain2_Plasmid1_Seq <tab> 621 <tab> S. boydii Strain Bar
```
The third column is optional, and used by KrakenHLL only when `--taxids-for-genomes` is specified for `krakenhll-build` to add new nodes in the taxonomy tree for the genome. If you'd like to have the sequences identifier in the taxonomy report, too, specifiy `--taxids-for-sequences` for `krakenhll-build` as well.

Finally, run `krakenhll-build`:
```
krakenhll-build --db DBDIR --taxids-for-genomes --taxids-for-sequences
```

Note that for custom databases with fewer sequences you might want to choose a smaller k (default: `--kmer-len 31`) and minimizer length (default: `--minimizer-len 15`).

### Custom databases with custom taxonomies

When using custom taxonomies, please provide `DBDIR/taxonomy/nodes.dmp` and `DBDIR/taxonomy/names.dmp` according to the format of NCBI taxonomy dumps.
