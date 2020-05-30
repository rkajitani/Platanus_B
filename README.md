# Platanus_B README.md

## Description
Platanus_B is a de novo assembler for isolated bacterial genomes. The features of this tool are as follows:
(1) It requires at least one Illuimina paired-end library. This library is useful for large-scale and/or high-resolution analysis. 
(2) It also can accept Oxford-Nanopore and/or PacBio long reads ("iterate" command).
(3) Implementing the iteration of sequence-extensions, gap-closing and error-removals, it can archive assemblies with contiguity and accuracy for many cases.
(4) As an utility function, it can combine multiple assemblies through the "combine" command.


## Version
v1.3.2

## Web site
<http://platanus.bio.titech.ac.jp/>

## Author
Rei Kajitani at Tokyo Institute of Technology wrote key source codes.  
Address for this tool: <platanus@bio.titech.ac.jp>


## Requirements
* OpenMP 
    - To compile the source code.

* Minimap2
    - <https://github.com/lh3/minimap2>
    - Included in this package.
    - Only required to use Oxford-Nanopore/PacBio long reads.

* Perl
    - To execute the scripts in this package, which  .
   
   
## Installation
```sh
make
cp platanus_b <installation_path>
```
Note that the absolute path of "sub_bin", which consists of Perl-scripts and minimap2 for "iterate" and "combine" commands, are written in the platanus_b executable.
Please re-complile this if you move sub_bin to another directory.
For macOS, please install OpenMP if a complilation fails. This can be installed using Homebrew with the command of "brew install libomp" or "brew install llvm". 


## Synopsis
### Inputs
* Illumina paired-end: PE_1.fq PE_2.fq (mandatory)
* Oxford Nanopor long-reads: ONT.fq (optional)

### Commands
```
platanus_b assemble -f PE_1.fq PE_2.fq 2>assemble.log
platanus_b iterate -c out_contig.fa -IP1 PE_1.fq PE_2.fq -ont ONT.fq 2>iterate.log
```

### Final output
    out_iterativeAssembly.fa



---
## Contig assembly usage
### Command
```sh
platanus_allee assemble [OPTIONS] 2>log
```
### Options
    -o STR               : prefix of output files (default out, length <= 200)
    -f FILE1 [FILE2 ...] : reads file (fasta or fastq, number <= 100)
    -k INT               : initial k-mer size (default 32)
    -K FLOAT             : maximum-k-mer factor (maximum-k = FLOAT*read-length, default  0.5)
    -s INT               : step size of k-mer extension (>= 1, default 10)
    -n INT               : initial k-mer coverage cutoff (default 0, 0 means auto)
    -c INT               : minimun k-mer coverage (default 1)
    -a FLOAT             : k-mer extension safety level (default 10.0)
    -u FLOAT             : maximum difference for bubble crush (identity, default 0)
    -d FLOAT             : maximum difference for branch cutting (coverage ratio, default 0.5)
    -e FLOAT             : k-mer coverage depth (k = initial k-mer size specified by -k) of homozygous region (default auto)
    -t INT               : number of threads (<= 100, default 1)
    -m INT               : memory limit for making kmer distribution (GB, >=1, default 16)
    -tmp DIR             : directory for temporary files (default .)
    -kmer_occ_only       : only output k-mer occurrence table (out_kmer_occ.bin; default off)
    -repeat              : mode to assemble repetitive sequences (e.g. 16s rRNA))


### Input format:
    Uncompressed and compressed (gzip or bzip2) files are accepted for -f option.

### Outputs:
    PREFIX_contig.fa
    PREFIX_kmerFrq.tsv

PREFIX is specified by -o
  
  
## Iteration of sequence-extension, gap-closeing and error-removal.
### Command
```sh
platanus_b iterate [OPTIONS] 2>log
```
### Options
    -o STR                             : prefix of output file and directory (do not use "/", default out, length <= 200)
    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format)
    -i INT                             : number of iterations (default 6)
    -l INT                             : -l value of "scaffold" step
    -u FLOAT                           : maximum difference for bubble crush (identity, default 0)
    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)
    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)
    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)
    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)
    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)
    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)
    -gc FILE1 [FILE2 ...]              : Guiding contig file; i.e. other assemblies, synthetic long-reads or corrected reads (fasta or fastq)
    -t INT                             : number of threads (default 1)
    -m INT                             : memory limit for making kmer distribution (GB, >=1, default 1)
    -tmp DIR                           : directory for temporary files (default .)
    -sub_bin DIR                       : directory for binary files which platanus_b use internally (e.g. minimap2) (default compilation-dir/sub_bin)
    -keep_file                         : keep intermediate files (default, off)
    -trim_overlap                      : trim overlapping edges of scaffolds (default, off)


### Input format:
   Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -p, -x and -X option.

### Outputs:
    PREFIX_allPhasedScaffold.fa (including sequences below)
    PREFIX_primaryBubble.fa
    PREFIX_secondaryBubble.fa
    PREFIX_nonBubbleHetero.fa
    PREFIX_nonBubbleOther.fa
    PREFIX_consensusInput.fa (for "consensus" command (-c))


PREFIX is specified by -o
  
  
## Consensus scaffold construction usage
### Command
```sh
platanus_allee consensus [OPTIONS] 2>log
```

### Options
    -o STR                             : prefix of output file (default out, length <= 200)
    -c FILE1 [FILE2 ...]               : input_scaffolds (fasta format)
    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)
    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)
    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)
    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)
    -p PAIR1 [PAIR2 ...]               : long-read file (PacBio, Nanopore) (fasta or fastq)
    -x PAIR1 [PAIR2 ...]               : linked-reads files (paired-ends, 10x Genomics) (interleaved, fasta or fastq)
    -X FWD1 REV1 [FWD2 REV2 ...]       : linked-reads files (paired-ends, 10x Genomics) (separate forward and reverse files, fasta or fastq)
    -t INT                             : number of threads (default 1)
    -tmp DIR                           : directory for temporary files (default .)
    -e FLOAT                           : coverage depth of homozygous region (default auto)
    -L INT                             : maximum fragment length of tag (10x Genomics) (default 200000)
    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default 32 64 96)
    -l INT                             : minimum number of links to scaffold (default 3)
    -mapper FILE                       : path of mapper executable file (default, minimap, only effective with -p option)
    -minimap2_sensitive                : sensitive mode for minimap2 (default, off; only effective with -p option)
    -no_partial                        : not close gaps partially, i.e. only close ones completely (default, off)

### Input format:
    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -ont, -p and -gc options.

### Outputs
    PREFIX_iterativeAssembly.fa (final)

PREFIX is specified by -o
  
  
## Dividing erroneous sequences usage
### Command
```sh
platanus_allee divide [OPTIONS] 2>log
```

### Options
    -o STR                             : prefix of output file and directory (do not use "/", default out, length <= 200)
    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format)
    -gc FILE1 [FILE2 ...]              : Guiding contig file; i.e. other assemblies, synthetic long-reads or corrected reads (fasta or fastq)
    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)
    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)
    -t INT                             : number of threads (default 1)
    -tmp DIR                           : directory for temporary files (default .)
    -sub_bin DIR                       : directory for binary files which platanus_b use internally (e.g. minimap2) (default compilation-dir/sub_bin)
    -no_gap_close                      : not close gaps by guiding contigs (default, false)
    -keep_file                         : keep intermediate files (default, off)

### Input format:
    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -gc, -ont and -p options.

### Outputs
    PREFIX_combined.fa (final)

PREFIX is specified by -o


---
## Notes
* Options related to run time
Although -t (number of threads) of all commands and -m (memory amount) of the "assemble" command are 
not mandatory to run, it is recommended to set the values adjusting your machine-environment.
These options may severely effect the run time.  
e.g.,  
Available number of threads and memory amount are 4 and 16GB, respectively.  
->  -t 4 -m 16

* Compressed input files
Both uncompressed and compressed (gzip or bzip2) FASTA/FASTQ files are accepted.
Formats are auto-detected. Internally, "file -bL", "gzip -cd" and "bzip2 -cd" commands, which can be
used in most of the UNIX OSs, are utilized.

* Paired-end (mate-pair) input  
The "phase" and "consensus" accept paired-end and/or mate-pair libraries. Paired libraries are 
classified into "inward-pair" and "outward-pair" according to the sequence direction. 
For file formats, separate and interleaved files can be input through -IP (-OP) and -ip (-op) 
options, respectively.

Inward-pair (usually called "paired-end", accepted in options "-IP" or "-ip"):

    FWD --->
        5' -------------------- 3'
        3' -------------------- 5'
                        <--- REV 

Outward-pair (usually called "mate-pair", accepted in options "-OP" or "-op"):

                        ---> REV 
        5' -------------------- 3'
        3' -------------------- 5'
    FWD <---

Example inputs:

    Inward-pair (separate, insert=300)   : PE300_1.fq PE300_2.fq
    Inward-pair (interleaved, insert=500): PE500_pair.fq
    Outward-pair (separate, insert=2k)   : MP2k_1.fa MP2k_2.fq

Corresponding options:

    -IP1 PE300_1_pair.fq PE300_2.fq \
    -ip2 PE500_pair.fq \
    -OP3 MP2k_1.fq MP2k_2.fq
