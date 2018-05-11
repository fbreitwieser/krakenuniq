Introduction
============

[KrakenHLL] is a taxonomic sequence classifier that assigns taxonomic labels to short DNA reads. It's based on [Kraken] ([Wood and Salzberg, Genome Biology 2014])
extended with unique k-mer counting using the HyperLogLog algorithm for better precision, as well as a couple other additional features.

KrakenHLL has in principle the same usage and system requirements as Kraken, as defined in the [Kraken manual]. This manual thus descibes the additional features of KrakenHLL (unique k-mer counting, hierarchical classification, taxonomy extended to include nodes for ) and how to switch from Kraken to KrakenHLL. [Pavian] works great for visualizing and analyzing KrakenHLL results! For a more detailed of the method please refer to the [KrakenHLL preprint]. Thank you!

[KrakenHLL]:                  http://ccb.jhu.edu/software/krakenhll/
[KrakenHLL GitHub repository]:   https://github.com/fbreitwieser/krakenhll
[Pavian]:   https://github.com/fbreitwieser/pavian
[KrakenHLL preprint]:            https://www.biorxiv.org/content/early/2018/04/03/262956   
[Kraken]:                     http://ccb.jhu.edu/software/kraken/
[Wood and Salzberg, Genome Biology 2014]:               http://genomebiology.com/2014/15/3/R46
[Kraken manual]:              http://ccb.jhu.edu/software/kraken/MANUAL.html



Switching from Kraken to KrakenHLL
==================================
KrakenHLL can be used as drop-in replacement to Kraken on a Kraken database. The first run will take longer as KrakenHLL builds its own taxonomy index and counts all k-mers in the database. Note that certain features, such as assembly and sequence identifications, require a full database download and build using KrakenHLL, but the unique k-mer counting works out of the box with a standard Kraken database. Note that `--report-file` on the command line is a required option. 
```
krakenhll --db DB --report-file REPORT_FILE --output KRAKEN_FILE
```
The output file of KrakenHLL is identical to Kraken. The report file has a couple of modifications - namely a header and three additional columns:

- kmers: number of unique k-mers
- dup: average number of times each unique k-mer has been seen
- cov: coverage of the k-mers of the clade in the database


Hierarchical read classification with multiple databases
========================================================
KrakenHLL allows using multiple databases hierarchically in order of confidence. In the following example each k-mer is matched first against the HOST, then the PROK, then the EUK_DRAFT database.
```
krakenhll --db HOST --db PROK --db EUK_DRAFT 
```
Note that the KrakenHLL databases need to share the same taxonomy database. If taxonomy nodes are added for genomes or sequences during the database build (see below), the databases have to be built consecutively using the taxDB file from the previous build.


Storing strain genomes with assembly project and sequence accessions
====================================================================
Kraken stores a NCBI taxonomic identifier for each k-mer in its database. This strategy worked well when new taxonomy IDs were assigned to each new microbial strain in GenBank. However, in 2014 the NCBI Taxonomy project stopped assigning new IDs to microbial strains; since then, only novel species get new taxonomy IDs (Federhen, et al., 2014). New microbial genomes, therefore, have the taxonomy ID of the species, or the taxonomy ID of a strain that was added bef¬ore 2014. Microbes that have been intensively surveyed, such as Escherichia coli or Salmonella spp., have hundreds of genomes indexed with the same taxonomy ID, and are thus indistinguishable by Kraken. An alternative way of identifying bacterial strains is to use the Bioproject, Biosample and Assembly accession codes (Breitwieser, et al., 2017). KrakenHLL thus adds new nodes to the taxonomy tree as children of the assigned taxon. A taxonomic node may also be added for each sequence; e.g., specific bacterial chromosomes or plasmids. Those new nodes in the taxonomy tree are given taxonomy IDs starting at 1,000,000,000. Having these extended nodes can help identify specific strains as well as bad database sequences (see Table 2 and Suppl. Table 6).

The additional information can be useful to detect the source of false positive identifications, too. In the reanalysis of the patient samples (Salzberg, et al., 2016) with database ‘std’, Salmonella enterica is detected in every sample with up to 233 reads. This species was not detected in the original analysis, and its ubiquity as well as a very low k-mer count hint that it is a false-positive hit or contaminant. If the only available information was the taxonomy ID, the search for the source of these hits would be difficult: There are 349 complete genomes in RefSeq for Salmonella enterica (taxonomy ID 28901) and still 23 complete genomes for the strain Salmonella enterica subsp. enterica serovar Typhimurium (taxonomy ID 90371). Supplementary Table 3 shows a part of the report KrakenHLL generated for PT8. Most of the reads going to Salmonella enterica hit one specific plasmid in one strain assembly. With standard Kraken output, neither the number of unique k-mers nor the sequence ID would have been known, and additional investigation such as re-alignment of the reads would have been required. 

| Reads |	Taxon Reads |	Kmers |	TaxID |	Rank |	Name |
|-|-|-|-|-|-|
|233	|0|	41|	590|	genus|	Salmonella
|233	|0|	41	|28901|	species|	·Salmonella enterica
|232	|0|	33	|59201|	subspecies|	··Salmonella enterica subsp. Enterica
|204	|0|	19	|90371|	no rank|	···Salmonella enterica subsp. enterica serovar Typhimurium
|203	|0|	8	|1000014850|	assembly|	····GCF_001617585.1 Salmonella enterica subsp. enterica serovar Typhimurium strain=RM9437
|203	|203	|8	|1000014852| 	sequence	| ·····NZ_CP014577.1 Salmonella enterica subsp. enterica serovar Typhimurium strain RM9437 plasmid pRM9437, complete sequence |

Table: Part of KrakenHLL output for PT8 (Salzberg, et al., 2016). Salmonella enterica is likely a false positive identification, and this is indicated by two factors: (1) the unique k-mer count is low. (2) The majority of reads hit a plasmid of one specific strain. 

To enable both features, call `krakenhll-build` with the options `--taxids-for-genomes` and `--taxids-for-sequences`. There is an important drawback to enabling these options: The pseudo-taxonomy IDs - e.g. 1000014850 in the table - are unique to the database build. Special precautions have to be taken when comparing results from different databases, or when using hierarchical mapping.

Integrating viral strain genomes in the database
================================================
The RefSeq project curates viral genomes (Brister, et al., 2015), which are included in the default databases of many metagenomics classifiers. RefSeq includes only one reference genome per viral species, and classifiers that use RefSeq (Kraken and others) therefore only consider those genomes. However, there are thousands of viral strain sequences in GenBank, and the chosen reference genome is often an established but old strain. For example, for HIV-1 the reference is a genome assembly from 1999 (AC GCF_000864765.1), and for JC polyomavirus the reference is the strain Mad1 (AC GCF_000863805.1) assembled in 1993. As many viruses exhibit high strain variabilty, including just the reference genomes in the Kraken database leads to a loss of sensitivity in the detection of strains.

KrakenHLL's database-building script includes the viral strain genomes from the NCBI viral genome resource (Brister, et al., 2015), which maintains a list of ‘neighbors’ to the viral reference genomes. This list has 112,148 sequences from viral strains, as compared to the 7497 viral genomes in RefSeq (as of October 2017). For example, there are over 2500 additional sequences for HIV-1, and over 640 for JC Polyomavirus. In total, these sequences add 100 million (+33%) novel k-mers to the database with k=31. Based on simulated reads from these viral sequences, 21.2% of the reads would not be classified when searching against a database which includes only the RefSeq viral reference genomes.

New taxonomy database format
============================
KrakenHLL has a new taxonomy format based on code from k-SLAM (Ainsworth, et al., 2017). The taxDB file lists the taxa in the following form:
```
Taxonomy ID<tab>Parent Taxonomy ID<tab>Rank<tab>Scientific Name
```
KrakenHLL reports all 27 ranks defined in the NCBI taxonomy, instead of just five abbreviated ranks in Kraken (‘D’ for superkingdom, ‘O’ for order, ‘P’ for phylum, ‘F’ for family, ‘G’ for genus, ‘S’ for species). For example, there are species groups and subgroups, subfamilies and varietas.

References
==========
-Ainsworth, D., et al. k-SLAM: accurate and ultra-fast taxonomic classification and gene identification for large metagenomic data sets. Nucleic Acids Res. 2017;45(4):1649-1656.
-Breitwieser, F.P. and Salzberg, S.L. Pavian: Interactive analysis of metagenomics data for microbiomics and pathogen identification. BioRxiv 2016.
-Brister, J.R., et al. NCBI viral genomes resource. Nucleic Acids Res 2015;43(Database issue):D571-577.
-Ertl, O. New Cardinality Estimation Methods for HyperLogLog Sketches. arXiv:1706.07290 2017.
-Flajolet, P., et al. HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm. In, AofA: Analysis of Algorithms. Juan les Pins, France: Discrete Mathematics and Theoretical Computer Science; 2007. p. 137-156.
-Heule, S., Nunkesser, M. and Hall, A. HyperLogLog in practice. 2013:683.
-Langmead, B. and Salzberg, S.L. Fast gapped-read alignment with Bowtie 2. Nat Methods 2012;9(4):357-359.
-Li, H., et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 2009;25(16):2078-2079.
-McIntyre, A.B.R., et al. Comprehensive benchmarking and ensemble approaches for metagenomic classifiers. Genome biology 2017;18(1).
-Menzel, P., Ng, K.L. and Krogh, A. Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nat. Commun. 2016;7:11257.
-Salzberg, S.L., et al. Next-generation sequencing in neuropathologic diagnosis of infections of the nervous system. Neurology(R) neuroimmunology & neuroinflammation 2016;3(4):e251.
-Whang, K.-Y., Vander-Zanden, B.T. and Taylor, H.M. A linear-time probabilistic counting algorithm for database applications. ACM Trans. Database Syst. 1990;15(2):208-229.


