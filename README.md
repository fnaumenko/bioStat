# bioStat
Statistical package for NGS data.<br>
It includes:<br>
 * [**cc**](#biocc)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; advanced correlation calculator for basic bioinformatics file formats<br>
 * [**calldist**](#calldist)&nbsp;&nbsp;&nbsp; calls fragment/read length distribution<br>
 * [**valign**](#valign)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; aligned reads verifier<br>
 * [**fqstatn**](#fqstatn)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; FastQ 'N' statistics calculator

## Usage
`biostat <command> [options] [<file>…]`<br>
or<br>
`<Command> [options] [<file>…]`<br><br>

## Installation
### Executable file

[Download Linux version](https://github.com/fnaumenko/biostat/releases/download/v1.0/biostat-Linux-x64.tar.gz) 
and extract **biostat** folder typing `tar -xf biostat-Linux-x64.tar.gz`.<br>
[Download Windows version](https://github.com/fnaumenko/biostat/releases/download/v1.0/biostat-Windows-x64.zip) 
and extract **biostat** folder using [WinRAR](https://www.win-rar.com/download.html?&L=0) or any other ZIP-archiver.

Alternative download in Linux:<br>
`wget https://github.com/fnaumenko/bioStat/releases/download/v1.0/biostat-Linux-x64.tar.gz`<br>

Add **biostat** to the PATH environment variable.

### Compiling in Linux
Required libraries:<br>
g++<br>
zlib

To compile from Git, type:
```
git clone https://github.com/fnaumenko/bioStat
cd bioStat
make
```
Alternative:
```
wget -O biostat.tar.gz https://github.com/fnaumenko/bioStat/archive/v1.0.tar.gz
tar -xf biostat.tar.gz
cd bioStat-1.0
make
```

---
## bioCC
fast advanced **C**orrelation **C**alculator for basic **bio**informatics data types.<br>
It computes Pearson’s and signal’s correlation coefficients for coverage, features and read densities.<br>
Program allows to obtain correlation coefficients for the whole genome, for each chromosome separately, 
and for predefined regions within the chromosomes. 
In the last case it optionally prints region coefficients frequency histogram. 
This, for example, makes it possible to correlate the densities precisely within peaks.<br>
**bioCC** uses a single-pass range-based correlation algorithm. It also is designed to treat a bunch of files at once.

### Usage
```
  biostat cc [options] file0 file1 ...
  biostat cc [options] -l|--list <file>
```
or
```
  bioCC [options] file0 file1 ...
  bioCC [options] -l|--list <file>
```

### Help
```
Input:
  -a|--align            input bed files are alignments. Ignored for bam and wig
  -g|--gen <name>       chromosome sizes file
  -l|--list <name>      list of multiple input files.
                        First (primary) file in list is comparing with others (secondary)
Processing:
  -c|--chr <name>       treat specified chromosome only
  -r|--cc <P,S>         correlation coefficient, in any order: P - Pearson, S - signal [P]
Region processing:
  -f|--fbed <name>      'template' ordinary bed file which features define compared regions.
                        Ignored for the ordinary beds
  -e|--ext-len <int>    length by which the features in primary file (for ordinary beds) or in
                        'template' (for alignments and wigs) will be extended in both directions
                        before treatment [0]
  -s|--ext-step <int>   step of extending features in primary bed file;
                        if 0 then no step calculation. For the ordinary beds only [0]
Output:
  -R|--pr-cc <LOC,TOT>  print coefficient, in any order:
                        LOC - for each chromosome, TOT - total [LOC]
  -B|--bin-width <float>
                        print frequency histogram with given bin width [0]
  -F|--fcc-sort [<RGN|CC>]
                        print region coefficients, sorted by: RGN - regions, CC - coefficients [CC]
  -V|--verbose  <LAC|NM|CNT|STAT>
                        set verbose level:
                        LAC  - laconic
                        NM   - file names
                        CNT  - file names and number of items
                        STAT - file names and statistics [NM]
  -o|--out              duplicate standard output to bioCC_out.txt file
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```

### Details

#### Input data
Dense continuous data (coverage) are compared using wiggle data in [WIG](https://genome.ucsc.edu/goldenpath/help/wiggle.html) format.
Features are compared using ordinary [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.
Read densities are compared using aligned DNA sequences (alignments) in [BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) or BED format.
The program recognizes the file format automatically by their extension (case-insensitive).

**Coverage**<br>
All type of coverage representation – [BedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html), 
[wiggle](https://genome.ucsc.edu/goldenpath/help/wiggle.html) variable step, fixed step – can be used in any combination.

**Features**<br>
Formally *alignment* (collection of reads) and *ordinary* bed (collection of features) have the same extention bed, 
but they are interpreted and handled differently. 
Since their automatic recognition is generally impossible, a special option is provided to indicate an *alignment*. 
See `-a|--align` option for more details.

**Read density**<br>
Each read is counted not by its length (this will be the read coverage, the corresponding WIG file can be obtained 
with standard tool such as [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html), 
[deepTools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html), 
[peakranger wigpe](http://ranger.sourceforge.net/manual1.18.html)), but by its 5’ position. 
Thus, two internal views of the coverage with span=1 are formed, which are then compared with each other.

BAM/BED do not have to be sorted by chromosomes, 
but all items (reads and features) belonging to a chromosome should be clustered together.  
The simplest way to meet this requirement is to use pre-sorted data. WIG items are sorted by definition.

#### Input data order
Comparable files can be represented both through program parameters, and by means of a file containing their names. 
In both cases the first file in a list – *primary* – is compared to the others – *secondary* – by turns. 
Only common chromosomes for the compared pair are considered (unless it is limited by the option `-c|--chr`). 

Be careful by using standard naming conventions like *abc?.bed*, *a\*.bed*. 
You must be sure that first of the lexically recognized files is really primary, and that all other files really need to be compared.

For more details about list of file names see `-l|--list` option.

#### Output
The program displays (and optionally duplicates to the file) the coefficients for each chromosome and the total, if it is specified.<br>
An example of the extended output for one chromosome, showing regional coefficients and their histograms histograms 
(here`$SM` stores the chromosome sizes filename): 
```
$ bioCC -tF -c 1 -f 2455_th09.bed -e 300 -B 0.1 -g $SM m36.wig 2455_M1.wig
template 2455_th09.bed: 65 features per chrom 1 (00:00)
Pearson CC between
m36.wig: 5909230 intervals per chrom 1 (00:04)
 and
2455_M1.wig: 2331020 intervals per chrom 1 (00:03)

#RGN    CC
43      -0.126805
63      0.52336
28      0.668467
...
54      0.990053
4       0.993
45      0.994721
BIN UP  COUNT
1       51
0.9     9
0.8     2
0.7     1
0.6     1
0.5     0
0.4     0
0.3     0
0.2     1
chr1    0.551579
00:07
```
This example compares the real (2455_M1.wig) and model (m36.wig) coverage at the regions specified in the optional template file 2455_th09.bed. 
In this case, the features of the template represent transcription factor Oct4 motifs, identified with a probability of 0.9. 
Before correlation, motives expand by 300 bp in both directions, thus the features indicate coverage peaks. 
The coefficients for each of 65 features (i.e. peaks) are sorted by their value. 
A histogram of their frequency distribution with a step of 0.1 is also displayed. 
In this case, we see that 81% of features have a correlation coefficient greater than 0.9, which indicates good model data.<br>
Data from different experiments can also be compared in this way.

### Options description
Note: enumerable option values are case-insensitive.<br>

`-a|--align`<br>
indicates that input bed files are *alignments*, so read density correlation would be performed.<br>
Each *alignment* line includes a strand character at the 6th tab position (separated by tabs), however, 
it is not forbidden to also include this symbol at the same position in a *regular bed*. 
Thus, automatic recognition of *alignment* and *ordinary* bed is generally not possible. 
To avoid possible ambiguity, the *alignment* is designated explicitly.<br>
If *alignments* are treated without this option, in most cases **bioCC** will print a cancel message and complete. 
But if input files have strand character at the 6th line positions and are not *alignments*, they will be compared wrongly.<br>
If bed files have strand character at the 6th line  and are treated without this option, a warning message will be printed.

`-g|--gen <name>`<br>
specifies chromosome sizes file. Required for BED and WIG files.

`-l|--list <file>`<br>
specifies a list of compared files. 
The list is a plain text file, in which each file name is located on a separate line.<br>
Lines beginning with ‘#’ are considered as comments and will be skipped.<br>
Empty lines are allowed.<br>
This option abolishes input files as parameters.

`-c|--chr <name>`<br>
treats specified chromosome only. 
The value `name` is the chromosome identifier; it is a number or character, for example, `10`, `X`.<br>
The indication of one chromosome reduces run time on 1.5-20 times depending on how far this chromosome is placed in an input data.<br>
For *ordinary* beds it has no time-improvement effect: any result appears quickly.<br>

`-r|--cc <P,S>`<br>
specifies correlation coefficient, in any order: `P` – Pearson, `S` – signal.<br>
See [Pearson and signal correlation](#pearson-and-signal-correlation).<br>
Default: Pearson

`-f|--fbed <file>`<br>
specifies 'template' *ordinary* bed file with features that defines compared regions.<br>
Correlation coefficients are calculated only within these areas (including their boundaries). 
Thus, the total coefficient for the chromosome will differ from the coefficient calculated without this option.<br>
An example of using this option is given in the [Output](#output) section.<br>
While template is specified, the output is limited only by the chromosomes presented in it.<br>
See also `-e|--ext-len` and `-b|--bin-width` options.<br>
This option is ignored for the *ordinary* bed files.

`-e|--ext-len <int>`<br>
specifies the value by which all features in a 'template' bed file or in a *primary ordinary* bed file should be stretched in both directions before comparison.<br>
If stretched features become intersected, they are joined.<br>
This option is mainly constructed for enriched region comparison while initial binding sites are represented by ‘template’. 
In case of *ordinary* bed, the *primary* file acts as a template. 
An example of using this option is given in the [Output](#output) section.<br>
Range: 0-1000<br>
Default: 0

`-s|--ext-step <int>`<br>
If set, activates the mode of consecutive calculation of the coefficients for stretching features in *primary ordinary* bed file with the stated step. 
The maximum value of the extension is limited by `--e|--ext-len` option.<br>
This option is topical for *ordinary* bed files only.<br>
Range: 0-500<br>
Default: 0 (no step calculation)

`-C|--pr-cc <LOC,TOT>`<br>
print coefficients, in any order: `LOC` - for each chromosome individually, `TOT` - total.<br>
If only one chromosome is specified, the total coefficient is not printed as an identical.<br>
Default: `LOC`.

`-B|--bin-width <float>`<br>
If set, forces to consolidate coefficients into bins, and print histogram values.<br>
Histogram shows count of coefficients (frequency) within some value ranges (bins). 
It is printed as a list of pairs *\<bin upper bound\>\<count in bin\>*. 
Negative coefficients have been turning to absolute during consolidation.<br>
This option defines the width of bin as a part of 1.<br>
Empty bins at the edges are not printed.<br>
An example of using this option is given in the [Output](#output) section.<br>
This option is topical only with option `-f|--fbed`.<br>
Range: 0-1<br>
Default: 0 (no consolidation)

`-F|--fcc-sort [<RGN|CC>]`<br>
If set, forces to print coefficients calculated for each region as a list of pairs *\<number-of-region\>\<coefficient\>*. 
`RGN` value prescribes list to be sorted by region’s number, `CC` – by coefficient.<br>
If both of the coefficients are declared, list is sorted by Pearson.<br>
First region number is 1.<br>
This option is topical only with option `-f|--fbed`.<br>

`-i|--info <LAC|NM|CNT|STAT>`<br>
outputs information about items (features/reads/intervals).<br>
`LAC`:&nbsp;&nbsp; laconic output. This value minimizes the output as possible to remain clear, e.g. for use in batch file.<br>
`NM`:&nbsp;&nbsp;&nbsp; brief output. Prints results and input file names.<br>
`CNT`:&nbsp;&nbsp; prints file names and number of all accepted items.<br>
`STAT`: prints item ambiguities statistics, if they exist. 
In the current version, statistics are displayed only for *ordinary* bed files, including template.<br>
Default: `NM`

`-o|--out`<br>
duplicates standard output to **bioCC_out.txt** file (except alarm messages).<br>
It is analogue of **tee** Linux command and is rather useful by calling **bioCC** under Windows.


### Pearson and signal correlation
For the data under consideration, the Pearson comparison is correct, since the data are linear in nature.<br>
While we consider coverage/density as the data distributed along regular dimension (scaled chromosome’s length), signal's method is appropriate as well.<br>
The only difference between them is the subtraction of the mean value in covariance on Pearson’s coefficient.

What consequences does it entail, and which coefficient is better to choose?

A. Coverage/density distributions.<br>
To be more clear there are 3 illustrations of pair of signals:<br>
![signal-Pearson](https://github.com/fnaumenko/bioStat/blob/master/pict/Signal-Pearson_50.png)<br>
Both of coefficients demonstrate value’s normalization independence (fig 1-3).<br>
Signal method is a bit more sensible, but it is sensitive to the mean amplitude (*DC offset* in terms of signal function), as we can see on (fig 2).
This means that the greater the background level of the compared distributions, the less relevant is the Signal method.<br>
For this reason, Pearson’s method is recommended for the distribution comparison.<br>
But if background’s level is considered part of the measure of similarity (for example, by comparing two replicas for noise level),
 in this case signal method would be preferable.

B.  Bed features.
Bed features can be accounted as discrete function accepted one of two values: 0 or 1. In that case signal method becomes inappropriate due to obvious reason: it counts intersections only. 
For example, by comparison two mutually complementary beds (while function1 have zero value every time when function2 have non-zero and vice versa), signal coefficient would be 0. 
Although the correct answer is -1.<br>
Thus, for the features only the Pearson method is correct.

---
## calldist

*Originally called callDist.* Calculates paired-end fragment size or read variable length distribution parameters.<br>
Examples of frequency profiles and recovered distributions of experimental datasets from NCBI database are shown 
in the ![Frag distributions figure](https://github.com/fnaumenko/bioStat/tree/master/pict/FragPE_distrs.png).<br>
Examples of frequency profiles and recovered read length distributions of experimental datasets from NCBI database are shown 
in the ![Read distributions figure](https://github.com/fnaumenko/bioStat/tree/master/pict/Read_distrs.png).

### Usage

`biostat calldist [options] <in-file>`<br>
or<br>
`calldist [options] <in-file>`

### Options
```
  -i|--inp <FRAG|READ>  input data to call distribution: FRAG - fragments, READ - reads [FRAG]
  -D|--dist <N,LN,G>    called distribution, in any order:
                        N – normal, LN – lognormal, G - Gamma [LN]
  -p|--pr-dist          print frequency distribution
  -o|--out [<name>]     duplicate standard output to specified file
                        or to <in-file>.dist if file is not specified
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```
### Details

#### Input
Fragment size distribution is called based on aligned DNA paired-end sequence in 
[BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm)/
[BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.<br>
Read size distribution is called based on original DNA sequence in [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format 
or aligned DNA sequence in BAM/BED format.<br>
Note that the number of mapped reads can be significantly less than the initial one, which can lead to distortion 
of the read distribution parameters relative to the original one (see, for example, cases 7-8 and 10-11 
in ![Read distributions figure](https://github.com/fnaumenko/bioStat/tree/master/pict/Read_distrs.png)).<br>
The program can also accept a file containing the finished distribution, in order to call its parameters. 
This is a plain text file with *.dist* extension, each line of which corresponds to one distribution point, 
i.e. a <frequency>-<size> pair. 
A similar file is produced when the `-p|--pr-dist` option is activated, and it can also be used as an input.

The program recognizes the file format automatically by their extension (case-insensitive).

#### Output
Called distribution parameters and Pearson correlation coefficient for the original and called distributions, 
calculated on the basis of the \<start of the sequence\> – \<the first frequency value less than 0.1% of the maximum\>.<br>
If the original distribution is assumed to be lognormal, 
the program also outputs the parameters and Pearson's coefficient for the normal distribution if it looks similar.<br>
An example of the extended output:
```
$ callDist -D ln,g -p 5278099.bam
5278099.bam: 4557867 fragments

	 PCC	relPCC	p1\*	p2\*\*	mode	exp.val
Lognorm	0.9811		5.775	0.4631	260	358.6
Gamma	0.95554	-2.6%	4.856	67.43	260	327.4

  \*p1 - mean, or alpha for Gamma
 \*\*p2 - sigma, or beta for Gamma

Original distribution:
length	frequency
60	1
70	1
72	1 
...
```
#### Options description

`-i|--inp <FRAG|READ>`<br>
sets the subject of distribution parameter calling: `FRAG` - fragments, `READ` - reads<br>
WARNING: in the Windows version, while trying to call fragment distribution with single-end BAM file, 
the program will crash silently instead of printing the corresponding message. 
This is due to an incorrectness in the external BamTools library being used. 
With single-end BED alignment **calldist** exits correctly, as well as with both formats under Linux.<br>
Default: `FRAG` for BAM/BED, `READ` for FASTQ

`-D|--dist <N,LN,G>`<br>
specifies the desired distribution type to call: `N` – normal, `LN` – lognormal, `G` – gamma. 
Types can be specified independently of each other and in any order.<br>
For each specified type, the called distribution parameters are displayed, as well as the Pearson correlation coefficient 
(PCC) with the original sequence. The coefficient is calculated on the basis from the beginning of the distribution 
to the first point with an ordinate that is less than 0.1% of the maximum.<br>
The types of distributions are sorted by PCC in descending order, and the ratio of the PCC to the maximum, in percent, 
is indicated as well.<br>
While a lognormal type is specified (default or explicit), the actual sequence is also automatically checked for normal distribution. 
Its parameters are displayed if its PCC exceeds the threshold of PCC_lognorm-2%. 
This is done because the lognormal distribution for certain parameters may differ slightly from the normal one. 
The final judgment is up to the user.<br>
Default: `LN` for BAM/BED, `N` for FASTQ

`-p|--pr-dist`
prints actual fragment/read length frequency distribution as a set of \<frequency\>-\<size\> pairs.
This allows to visualize the distribution using some suitable tool such as Excel.

`-o|--out [<name>]`<br>
duplicates standard output to specified file (except alarm messages).<br>
If file is not specified, duplicates output to file with name, 
constructed as input file short name (without path and extension) with addition of the extension *.freq*.<br>
If, in addition, the input file already has the *.freq* extension, then the "_out" suffix is added to the name.<br>
It is an analogue of the **tee** Linux command and is constructed rather for the execution under Windows.

---
## vAlign
**V**erify **Align**ment is a fast verifier of reads forming the aligned DNA sequence, 
which is recalled from an artificial [FastQ](https://en.wikipedia.org/wiki/FASTQ_format) sequence. 
It compares the original and mapped coordinates of each read and prints statistics of right and wrong mappings.

To do this each read in an initial artificial sequence should keep its location as a part of its name. 
Read’s name format should be \<some_text\>:chr\<ID\>:\<original_start_pos\>.\<uniq_number\> for sigle-end reads and 
\<some_text\>:chr\<ID\>:\<original_frag_start_pos\>-\<original_frag_end_pos\>.\<uniq_number\>/\<mate\> for paired-end reads.<br>
Such template sequence can be generated by [**isChIP**](https://github.com/fnaumenko/isChIP) software. 

### Usage
`biostat valign [options] -g|--gen <name> <in-file>`<br>
or<br>
`vAlign [options] -g|--gen <name> <in-file>`

### Options:
```
Treatment:
  -g|--gen <name>       reference genome library or single nucleotide sequence. Required
  -c|--chr <name>       treat specified chromosome only. For reference genome only
  --min-scr <int>       score threshold for treated reads
  --char-case <OFF|ON>  recognize uppercase and lowercase characters in template and test
                        as different [OFF]
  -o|--out [<name>]     duplicate standard output to specified file
                        or to <in-file>_valign.txt if file is not specified
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help WARNING     print usage information and exit
```
#### Input
Aligned DNA single- or paired end sequence (with fixed read length) 
in [BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) 
or [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.<br>
The program recognizes the file format automatically by their extension (case-insensitive).<br>
The sequence does not have to be sorted by chromosomes, but all reads belonging to a chromosome should be clustered together. 
The simplest way to meet this requirement is to use pre-sorted data.

#### Output
An example of the output(here $G stores the reference genome directory):
```
$ vAlign -to -g $G testVAlign.bam
<in-file> testVAlign.bam: 24534298 reads, from which
	37738 (0.1538%) duplicated reads; accepted
	total accepted: 24534298 (100%) reads (00:18)
chr1
mism	readCnt	quality
precise	872019	1
0	9301	1
1	3	1
2	2	1
...
47	418	1
48	147	1
49	81	1
50	102	1
total reads per chrom 1:	1754353	92.8572%	1889303
reads per different chroms:	134950	7.14285%
chr2
...
```
`mism` – mismatches – means the number of erroneous nucleotides in a read (limited by the length of the read);<br>
`readCnt` – number of reads with given mismatches number;<br>
`quality` – average quality value in relative units for given reads;<br>
`precise` means the number of reads mapped to true coordinates without mismatches. 
Zero mismatches denotes reads without mismatches but mapped to "false" position (different from the original).

#### Options description
`-g|--gen <name>`<br>
specifies reference genome library or single nucleotide sequence.<br>
Genome library is a directory containing nucleotide sequences for each chromosome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.<br>
If `name` is a .fa[.gz] file, **vALign** accepts the corresponding chromosome as the only treated.<br>
Otherwise first the program searches for .fa files in the directory `name`. If there are no such files in this directory, 
**vALign** searches for .fa.gz files.<br>
If chromosome is specified by option `–c|--chr`, the program searches for the corresponding .fa[.gz] file.

One can obtain a genome library in UCSC ftp://hgdownload.soe.ucsc.edu/goldenPath/ or in Ensemble ftp://ftp.ensembl.org/pub/release-73/fasta storage. 
In the second case please copy genomic sequences with the same masked type only, f.e. unmasked (‘dna'), since program does not recognise mask’s types.<br>
This option is required.

`-c|--chr <name>`<br>
treats specified chromosome only. The value `name` is the chromosome identifier; 
it is a number or character, for example, 10 or X.
The indication of one chromosome reduces run time on 1.5-20 times depending on how far this chromosome is placed in an alignment. 
This simultaneously reduces statistical completeness.

`--min-scr <int>`<br>
specifies score threshold for treated reads. Reads with the score equal or less then specified will be ignored.<br>
Default: all reads are accepted.

`--char-case <OFF|ON>`<br>
turns off/on recognition of uppercase and lowercase characters in template and test as different.<br>
Default: `OFF`.

`-o|--out`<br>
duplicates standard output to specified file (except alarm messages). If file is not specified, duplicates output to file with name, 
constructed as input alignment short name (without path and extension) with addition of the suffix “_valign.txt”.<br>
It is an analogue of the **tee** Linux command and is constructed rather for the execution under Windows.

---
## fqStatN

Calculates the statistics of occurrence of ambiguous code N, and patterns of reads including N 
in the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file.<br>
These statistics helps to better evaluate the quality of the sequencer output.<br>
Only fixed-length reads are accepted.

### Usage
`biostat fqstatn [options] <sequence>`<br>
or<br>
`fqStatN [options] <sequence>`

### Options:
```
  -o|--out [<name>]     duplicate standard output to specified file
                        or to <sequence>_statn.txt if file is not specified
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit 
```

### Details

#### Input
DNA sequence in [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format.

#### Output
The program displays the frequency of occurrence of 'N' in the position in the read, as well as the frequency of the template of the reads containing code N.
Example  of output:
```
'N' POSITION STATISTICS
pos     count   % of total 'N'
------------------------------
 0      29588   71.1%
 3      4314    10.4%
20      69      0.166%
21      49      0.118%
...
47      256     0.615%
48      33      0.0793%
49      1904    4.58%

READ TEMPLATE STATISTICS
position  10        20        30        40        50    patterns count
01234567890123456789012345678901234567890123456789
----------------------------------------------------------------------
N.................................................    29566     0.172%
...N..............................................     4313     0.0251%
..............................NN..................       83     <0.001%
.............................NNNN.NNNN.........N.N       10     <0.001%
...
...............................N..NN..............        3     <0.001%
...............................N..N...............        2     <0.001%
.........................N.....N..NN.N...........N        1     <0.001%

'N' relative to the total number of nucleotides: 0.00483%
Reads that include 'N' relative to the total number of reads: 0.211%
```

---
If you face to bugs, incorrect English, or have commentary/suggestions, please do not hesitate to write me on fedor.naumenko@gmail.com
