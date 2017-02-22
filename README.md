Kraken taxonomic sequence classification system with Unique K-mer Counting
===============================================

[Kraken](https://github.com/DerrickWood/kraken) is a fast taxonomic classifier for metagenomics data. This project, kraken-hll, adds some additional functionality - most notably a unique k-mer count. Spurious identifications due to sequence contamination in the dataset or database often leads to many reads, however they usually cover only a small portion of the genome. 

kraken-hll adds two additional columns to the Kraken report - total number of k-mers observed for taxon, and the total number of unique k-mers observed for taxon (columns 3 and 4, resp.). 

Here's a small example of a classification against a viral database with k=25. There are three species identified by just one read - Enterobacteria phage BP-4795, Salmonella phage SEN22, Sulfolobus monocaudavirus SMV1. Out of those, the identification of Salmonella phage SEN22 is the strongest, as there read was matched with 116 k-mers that are unique to the sequence, while the match to Sulfolobus monocaudavirus SMV1 is only based on a single 25-mer.

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

For usage, see `kraken_hll --help`. Note that you can use the same database as Kraken with one difference - instead of the files `DB_DIR/taxonomy/nodes.dmp` and `DB_DIR/taxonomy/names.dmp` than kraken relies upon, `kraken-hll` needs the file `DB_DIR/taxDB`. This can be generated with the script `build_taxdb`: `KRAKEN_DIR/build_taxdb DB_DIR/taxonomy/names.dmp DB_DIR/taxonomy/nodes.dmp > DB_DIR/taxDB`. The code behind the taxDB is based on [k-SLAM](https://github.com/aindj/k-SLAM).

### Differences to `kraken`
 - Use `kraken_hll --report-file FILENAME ...` to write the kraken report to `FILENAME`.
 - Use `kraken_hll --db DB1 --db DB2 --db DB3 ...` to first attempt, for each k-mer, to assign it based on DB1, then DB2, then DB3. You can use this to prefer identifications based on DB1 (e.g. human and contaminant sequences), then DB2 (e.g. completed bacterial genomes), then DB3, etc. Note that this option is incompatible with `kraken_hll-build --generate-taxonomy-ids-for-sequences` since the taxDB between the databases has to be absolutely the same.
 - Add a suffix `.gz` to output files to generate gzipped output files

### Differences to `kraken-build`
 - Use `kraken_hll-build --generate-taxonomy-ids-for-sequences ...` to add pseudo-taxonomy IDs for each sequence header. An example for the result using this is in the ouput above - one read has been assigned specifically to `KC207814.1 Human herpesvirus 4 strain Mutu, complete genome`.
 - `seqid2taxid.map` mapping sequence IDs to taxonomy IDs does NOT parse or require `>gi|`, but rather the sequence ID is the header up to just before the first space
