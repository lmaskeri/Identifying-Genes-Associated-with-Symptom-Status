## Installation

### Bioconda
If you use [Conda](https://conda.io/docs/install/quick.html)
you can use the [Bioconda channel](https://bioconda.github.io/):
```
conda install -c conda-forge -c bioconda -c defaults prokka
```

### Brew
If you are using the [MacOS Brew](http://brew.sh/) 
or [LinuxBrew](http://brew.sh/linuxbrew/) packaging system:
```
brew install brewsci/bio/prokka
```

### Docker
Maintained by https://hub.docker.com/u/staphb
```

docker pull staphb/prokka:latest
docker run staphb/prokka:latest prokka -h
```

### Singularity
```
singularity build prokka.sif docker://staphb/prokka:latest
singularity exec prokka.sif prokka -h
```

### Ubuntu/Debian/Mint
```
sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
sudo cpan Bio::Perl
git clone https://github.com/tseemann/prokka.git $HOME/prokka
$HOME/prokka/bin/prokka --setupdb
```

### Centos/Fedora/RHEL
```
sudo yum install git perl-Time-Piece perl-XML-Simple perl-Digest-MD5 perl-App-cpanminus git java perl-CPAN perl-Module-Build
sudo cpanm Bio::Perl
git clone https://github.com/tseemann/prokka.git $HOME/prokka
$HOME/prokka/bin/prokka --setupdb
```

### MacOS
```
sudo cpan Time::Piece XML::Simple Digest::MD5 Bio::Perl
git clone https://github.com/tseemann/prokka.git $HOME/prokka
$HOME/prokka/bin/prokka --setupdb
```

## Test

* Type `prokka` and it should output its help screen.
* Type `prokka --version` and you should see an output like `prokka 1.x`
* Type `prokka --listdb` and it will show you what databases it has installed to use.

## How to use
Prokka is very simple to use, just prepare the fasta format file of the nucleic acid sequence that needs to be annotated. Here I will use the reference genome egd of *Listeria monocytogenes* as an example.
```
# Download genome data
$ esearch -db nuccore -query "1639[txid] AND egd" | efetch -format fasta > egd.fasta

# Annotated genome sequence
$ prokka egd.fasta
```
`--listdb` parameter can view prokka database information.

```

# View available databases
$ prokka --listdb

Looking for databases in: ...
* Kingdoms: Archaea Bacteria Mitochondria Viruses
* Genera: Enterococcus Escherichia Staphylococcus
* HMMs: HAMAP
* CMs: Bacteria Viruses
```
If you need a more complete sequence annotation, and the results can be more convenient for submission, you can add corresponding parameters to standardize annotations and call different tools to scan ncRNA, etc.

```
# Comment Listeria monocytogenes standard strain egd
$ prokka --outdir egd --prefix egd --addgenes --addmrna --compliant --centre CDC --genus Listeria --species "Listeria monocytogenes" --strain egd
--kingdom Bacteria --usegenus --cpus 4 --rfam --rnammer --force egd.fasta
```

## Output Files

| Extension | Description |
| --------- | ----------- |
| .gff | This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV. |
| .gbk | This is a standard Genbank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-Genbank, with one record for each sequence. |
| .fna | Nucleotide FASTA file of the input contig sequences. |
| .faa | Protein FASTA file of the translated CDS sequences. |
| .ffn | Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA) |
| .sqn | An ASN1 format "Sequin" file for submission to Genbank. It needs to be edited to set the correct taxonomy, authors, related publication etc. |
| .fsa | Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines. |
| .tbl | Feature Table file, used by "tbl2asn" to create the .sqn file. |
| .err | Unacceptable annotations - the NCBI discrepancy report. |
| .log | Contains all the output that Prokka produced during its run. This is a record of what settings you used, even if the --quiet option was enabled. |
| .txt | Statistics relating to the annotated features found. |
| .tsv | Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product |

## .tsv files
The tsv file sorts the results of the comments in the order of locus_tag.

```
#Display tsv content (only in this case)
$ head PROKKA_04172021.tsv

locus_tag ftype length_bp gene EC_number COG product
DJECODEN_00001 CDS 1356 dnaA COG0593 Chromosomal replication initiator protein DnaA
DJECODEN_00001 gene 1356 dnaA
DJECODEN_00001 mRNA 1356 dnaA
DJECODEN_00002 CDS 1146 dnaN COG0592 Beta sliding clamp
DJECODEN_00002 gene 1146 dnaN
DJECODEN_00002 mRNA 1146 dnaN
DJECODEN_00003 CDS 1344 yeeO_1 COG0534 putative FMN/FAD exporter YeeO
DJECODEN_00003 gene 1344 yeeO_1
DJECODEN_00003 mRNA 1344 yeeO_1

```
locus_tag: the locus name of the annotated gene
ftype: Type, the default is CDS, if you turn on the --addgenes and --addmrna parameters, it will distinguish between RNA and coding genes
length_bp: sequence length
gene: According to the gene name corresponding to the database annotation, if it is multiple copies, use_1_2, etc. to distinguish
EC_number: EC value corresponding to the gene
COG: COG corresponding to the gene
product: The protein product encoded by the gene
