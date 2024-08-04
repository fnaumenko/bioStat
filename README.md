# bioStat
Cross-platform statistical package for NGS data.<br>
It includes:<br>
 * [**cc**](#biocc)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;advanced correlation calculator for basic bioinformatics file formats<br>
 * [**fgstest**](#fgstest)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;features Gold Sandard test<br>
 * [**calldist**](#calldist)&nbsp;&nbsp;&nbsp;&nbsp;calls fragment/read length distribution<br>
 * [**valign**](#valign)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;aligned reads verifier<br>
 * [**fqstatn**](#fqstatn)&nbsp;&nbsp;&nbsp;&nbsp;fastQ 'N' statistics calculator

version 3.0 in progress

## Usage
`biostat <command> [options] [<file>…]`<br>
or<br>
`<Command> [options] [<file>…]`<br><br>
Yes, the utilities can be invoked separately. They all are packed in one archive bioStat.

*Notes for all utilities:*<br>
Enumerable option values are case-insensitive.<br>
Single letter options with missing values can be merged.<br>
Compressed input files in gzip format are acceptable.

## Installation
### Executable file

[Download Linux version](https://github.com/fnaumenko/biostat/releases/download/v2.1/biostat-Linux-x64.tar.gz) 
and extract **biostat** folder typing `tar -xf biostat-Linux-x64.tar.gz`.<br>
[Download Windows version](https://github.com/fnaumenko/biostat/releases/download/v2.1/biostat-Windows-x64.zip) 
and extract **biostat** folder using [WinRAR](https://www.win-rar.com/download.html?&L=0) or any other ZIP-archiver.

Alternative download in Linux:<br>
`wget https://github.com/fnaumenko/bioStat/releases/download/v2.1/biostat-Linux-x64.tar.gz`<br>

Add **biostat** to the PATH environment variable.

### Compiling
Required components:<br>
g++ (Linux)<br>
cmake<br>
zlib (optionally)

To compile from Git, type:
```
git clone https://github.com/fnaumenko/bioStat
cd bioStat
git submodule init
git submodule update
cmake .
cmake --build . --target ALL_BUILD --config Release
```
Alternative:
```
wget -O biostat.tar.gz https://github.com/fnaumenko/bioStat/archive/v2.0.tar.gz
tar -xf biostat.tar.gz
cd bioStat-2.1
git submodule init
git submodule update
cmake src
cmake --build . --target ALL_BUILD --config Release
```
Quick check the build:
```
build/release/biocc -v
```

---
## bioCC
fast advanced **C**orrelation **C**alculator for basic **bio**informatics data types.<br>
It computes Pearson’s correlation coefficients for coverage, features and read densities.<br>
Program allows to obtain correlation coefficients for the whole genome, for each chromosome separately, 
and for predefined regions within the chromosomes.<br>
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
  -c|--chr <name>       treat specified chromosome only
  -o|--overl <OFF|ON>   allow (and merge) overlapping features. For the ordinary beds only [OFF]
  -d|--dup <OFF|ON>     allow duplicate reads. For the alignments only [ON]
  -l|--list <name>      list of multiple input files.
                        First (primary) file in list is comparing with others (secondary)
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
  -V|--verbose  <LAC|NM|CNT|STAT>s
                        set verbose level:
                        LAC  - laconic
                        NM   - file names
                        ITEM - file names and number of items
                        STAT - file names and items statistics [NM]
  -O|--out [<name>]     duplicate standard output to specified file
                        or to bioCC.output.txt if <name> is not specified
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```

### Details

#### Input data
**Dense continuous data** (*coverage*) are compared using wiggle data in WIG format.<br>
Any WIG type – [BedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html), 
[wiggle](https://genome.ucsc.edu/goldenpath/help/wiggle.html) variable step, fixed step – can be used in any combination.<br>
**Features** are compared using *ordinary* [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.<br>
**Read densities** are compared using aligned DNA sequences (*alignments*) in [BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) or BED format.<br>
Each read is counted by its 5’ position. Thus, internal views of the coverage with span = 1 are formed and compared.<br>
To correlate the actual read coverages , translate alignments into WIG format using standard tools such 
[bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html), 
[deepTools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html), 
[peakranger wigpe](http://ranger.sourceforge.net/manual1.18.html)).

The program recognizes the file format automatically by their extention (case-insensitive).<br>
To distinguish an *alignment* in BED format from *ordinary* BED apply `-a|--align` option.<br>
Datasets must be sorted.

#### Input data order
Comparable files can be represented both through program parameters, and by means of a file containing their names. 
In both cases the first file in a list – *primary* – is compared to the others – *secondary* – by turns. 
Be careful by using standard naming conventions like *abc?.bed*, *a\*.bed*. 
You must be sure that first of the lexically recognized files is really primary, and that all other files really need to be compared.<br>
For more details about list of file names see `-l|--list` option.

#### Output
The program displays (and optionally duplicates to the file) the coefficients for each chromosome and/or the total one.<br>
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
This example compares the real and model coverage at the regions specified in the optional template file. 
In this case, the features of the template represent proved transcription factor motifs. 
Before correlation, motives expand by 300 bp in both directions in order to spot peaks area. 
The coefficients for each of 65 features (i.e. peaks) are sorted by their value. 
A histogram of their frequency distribution with a step of 0.1 is also displayed. 
In this case, we see that 81% of features have a correlation coefficient greater than 0.9.

### Options description

`-a|--align`<br>
indicates that input bed files are *alignments*, so the read density correlation would be performed.

`-g|--gen <name>`<br>
specifies chromosome sizes file. Required for BED and WIG files.

`-c|--chr <name>`<br>
treats specified chromosome only.<br>
`name` identifies chromosome by number or character, e.g. `10` or `X`. Character is case-insensitive.<br>
Specifying one chromosome reduces processing time of multi-chromosomal data by 2-20  times.<br>
*Ordinary* beds are treated quickly in any case.

`-o|--overl <OFF|ON>`<br>
rejects or accept overlapping features for processing.<br>
In the first case, 'overlapping chains' break. 
This means, for example, that if feature #2 overlaps feature #1 and feature #3 overlaps feature #2 but not #1, 
then feature #3 will be accepted.<br>
In the second case, overlapping features are joined.<br>
Adjacent features are also considered to be overlapping.<br>
Makes sense for the *ordinary* beds only, including template.<br>
Default: `ON`

`-d|--dupl <OFF|ON>`<br>
rejects or accept duplicated reads for processing.<br>
Makes sense for the *alignments only*.<br>
Default: `ON`

`-l|--list <file>`<br>
specifies a list of compared files. 
The list is a plain text file, with one file name per line.<br>
Lines starting  with ‘#’ are treated as comments and are ignored, as well as empty lines.<br>
This option abolishes input files as parameters.

`-f|--fbed <file>`<br>
specifies 'template' *ordinary* bed file with features that defines compared regions within chromosomes.<br>
Correlation coefficients are calculated only within these areas (including their boundaries). 
It is the same as we cut these data areas, joined them and compare.<br>
Data for chromosomes not presented in 'template' are ignored.<br>
An example of using this option is given in the [Output](#output) section.<br>
See also `-e|--ext-len` and `-b|--bin-width` options.<br>
This option is ignored for the *ordinary* bed files.

`-e|--ext-len <int>`<br>
specifies the extending value by which all features in a 'template' bed file or in a *primary ordinary* bed file should be stretched in both directions before comparison.<br>
If stretched features become intersected, the extending value is limited to keep the minimum possible gap between nearest features.<br>
An example of using this option is given in the [Output](#output) section.<br>
Range: 0-2000<br>
Default: 0

`-s|--ext-step <int>`<br>
If set, activates the mode of consecutive calculation of the coefficients for stretching features in *primary ordinary* bed file with the stated step. 
The maximum value of the extension is limited by `--e|--ext-len` option.<br>
This option is topical for *ordinary* bed files only.<br>
Range: 0-500<br>

`-C|--pr-cc <LOC,TOT>`<br>
print coefficients, in any order: `LOC` - for each chromosome individually, `TOT` - total.<br>
Default: `LOC`.

`-B|--bin-width <float>`<br>
If set, forces to consolidate coefficients into bins, and print histogram values.<br>
Histogram shows count of coefficients (frequency) within some value ranges (bins). 
It is printed as a list of pairs *\<bin-upper-bound\>\<count-in-bin\>*. 
Negative coefficients have been turning to absolute during consolidation.<br>
This option defines the width of bin as a part of 1.<br>
Empty high and low bins are not printed.<br>
An example of using this option is given in the [Output](#output) section.<br>
This option is topical only with option `-f|--fbed`.<br>
Range: 0-1<br>

`-F|--fcc-sort [<RGN|CC>]`<br>
If set, forces to print coefficients calculated for each region as a list of pairs *\<number-of-region\>\<coefficient\>*. 
`RGN` value prescribes list to be sorted by region’s number, `CC` – by coefficient.<br>
Regions are numbered starting from 1.<br>
This option is topical only with option `-f|--fbed`.<br>

`V|--verbose <LAC|NM|CNT|STAT>`<br>
sets verbose level:<br>
`LAC`:&nbsp;&nbsp; laconic output. This value minimizes the output as possible to remain clear, e.g. for use in batch file.<br>
`NM`:&nbsp;&nbsp;&nbsp; besides results prints input file names.<br>
`ITEM`: besides results prints input file names and number of presented/accepted items.<br>
`STAT`: besides results prints input file names and item ambiguities statistics, if exist.<br>
Default: `NM`

`-O|--out [<name>]`<br>
duplicates standard output to specified file (except alarm messages).<br>
If <name> is not specified, duplicates output to **bioCC.output.txt** file.<br>
If the <name> denotes an existing folder, the output file is created in it according to the rule described above.<br>
It is an analogue of the **tee** Linux command and is constructed rather for the execution under Windows.

---
## FGStest
**F**eatures **G**old **S**andard statistical **Test**

### Usage

`biostat fgstest [options] -S|--sample <name> <in-file>`<br>
or<br>
`fgstest [options] -S|--sample <name> <in-file>`

### Options
```
Input:
  -c|--chr <name>       treat specified chromosome only
  -S|--sample <name>    sample file. Required
  -C|--min-cdev <int>   threshold centre deviation for writing a test feature to an issues file [10]
  -W|--min-wdev <float> threshold width deviation for writing a test feature to an issues file [0]
  -s|--min-scr <float>  threshold score for taking sample features into accounts [0]
  -e|--expand <int>     expand sample features [0]
Output:
  -I|--issues [<name>]  output locused issues to <name>.bed file
                        or to <in-file>.issues.bed file if <name> is not specified
  -O|--out [<name>]     duplicate standard output to <name> file
                        or to <in-file>.output.txt file if <name> is not specified
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```
  
#### Input
[BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file containing test sites of interest<br>
Typically this is the result of peak detectors.

#### Output
In progress.

#### Options description

`-S|--sample <name>`<br>
name of [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file containing the real sites of interest ('gold standard')<br>
Required.

`-c|--chr <name>`<br>
treats specified chromosome only.<br>
`name` identifies chromosome by number or character, e.g. `10` or `X`. Character is case-insensitive.

`-C|--min-cdev <int>`<br>
in progress.<br>
Default: 10.

`-W|--min-wdev <float>`<br>
in progress.<br>
Default: 0.

`-s|--min-scr <float>`<br>
specifies template features score threshold for testing.<br>
Default: 0.

`-e|--expand <int>`<br>
in progress.<br>
Default: 0.

`-I|--issues [<name>]`<br>
specifies output file containing issued features in BED format.<br>
The name of output file is constructed as *`name`.bed*, possible extension in <name> is truncated.<br>
If `name` is not specified, the output file name is constructed as *`in-file`.issues.bed*,
 where `in-file` is a program input file. Its extention is also truncated.<br>
If `name` denotes an existing folder, the output file is created in it using `in-file` template, according to the rule described above.<br>
File contents:<br>
The first five fields are standard.<br>
The feature name is the designation of the issue:<br>
*FN* - False Negative case<br>
*FP* - False Positive case<br>
*cD* - critical centre deviation case<br>
*wD* - critical width deviation case<br>
The sixth and seventh fields are non-standard.<br>
The 6th field contains critical deviation values: integer for *cD* and floating point for *wD*. 
In the latter case, it is the ratio of the test feature width to the reference one<br>
The 7th field contains the issue locus (coordinates for display in the genome browser).<br>
The issues file is auxiliary and is intended mainly for quick viewing of problematic cases in the genome browser.

`-O|--out [<name>]`<br>
duplicates standard output to specified file (except alarm messages).<br>
If <name> is not specified, the output file name is constructed as *<in-file>.ioutput.txt*,
 where <in-file> is a program input file. Its extention is also truncated.<br>
If <name> denotes an existing folder, the output file is created in it using <in-file> template, according to the rule described above.<br>
It is an analogue of the **tee** Linux command and is constructed rather for the execution under Windows.

---
## callDist

Calculates paired-end fragment size or read variable length distribution parameters.<br>
Can check it against *normal*, *lognormal* and *gamma* distributions.<br>
Examples of frequency profiles and recovered distributions of experimental datasets from NCBI database are shown 
in the ![Frag distributions figure](https://github.com/fnaumenko/bioStat/tree/master/pict/FragPE_distrs.png).<br>
Examples of frequency profiles and recovered read length distributions of experimental datasets from NCBI database are shown 
in the ![Read distributions figure](https://github.com/fnaumenko/bioStat/tree/master/pict/Read_distrs.png).

### Usage

`biostat calldist [options] <in-file>`<br>
or<br>
`callDist [options] <in-file>`

### Options
```
Input:
  -i|--inp <FRAG|READ>  input data to call distribution: FRAG - fragments, READ - reads [FRAG]
  -c|--chr <name>       treat specified chromosome only
  -D|--dist <N,LN,G>    called distribution (can be combined in any order):
                        N - normal, LN - lognormal, G - Gamma [LN]
  -d|--dup <OFF|ON>     allow duplicates [ON]
Processing:
  -p|--pr-dist          print obtained frequency distribution to file
  -s|--stats            print input item issues statistics
  -O|--out [<name>]     duplicate standard output to specified file
                        or to <in-file>.dist if <name> is not specified
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```
### Details

#### Input
*Fragment* size distribution is called based on aligned sorted DNA paired-end sequence in 
[BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm)/
[BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.<br>
*Read* size distribution is called based on original DNA sequence in [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format 
or aligned sorted DNA sequence in BAM/BED format.<br>
The number of mapped reads can be significantly less than the initial one, which can lead to distortion 
of the read distribution parameters relative to the original one (see, for example, cases 7-8 and 10-11 
in ![Read distributions figure](https://github.com/fnaumenko/bioStat/tree/master/pict/Read_distrs.png)).<br>
The program can also accept a file containing the finished distribution, in order to call its parameters.<br>
This is a plain text file with *.dist* extention, each line of which corresponds to one distribution point, 
i.e. a pair \<size\>&#x2011;\<frequency\>. Both values should be integers.<br>
The first lines of the file that do not contain such a pair are ignored.<br>
A similar file is produced when the `-p|--pr-dist` option is activated, and it can also be used as an input.<br>
Input file with *.dist* extention ignores `-p|--pr-dist` option (but not `-O|--out` one).

The program recognizes the file format automatically by their extention (case-insensitive).

#### Output
Called distribution parameters and Pearson correlation coefficient (PCC) for the original and called distributions, 
calculated on the basis of the \<start of the sequence\>–\<the first frequency value less than 0.1% of the maximum\>.<br>
If the original distribution is assumed to be lognormal, 
the program also outputs the parameters and PCC for the normal distribution if it looks similar.<br>
An example of the output:
```
$ callDist -D ln,g -p 5278099.bam
5278099.bam: 4557867 fragments

	 PCC	relPCC	p1*	p2**	mode	exp.val
Lognorm	0.9811		5.775	0.4631	260	358.6
Gamma	0.95554	-2.6%	4.856	67.43	260	327.4

  *p1 - mean, or alpha for Gamma
 **p2 - sigma, or beta for Gamma

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
This option is topical for BAM/BED files only.<br>
Default: `FRAG` for BAM/BED, `READ` for FASTQ

`-c|--chr <name>`<br>
treats specified chromosome only.<br>
`name` identifies chromosome by number or character, e.g. `10` or `X`. Character is case-insensitive.<br>
Specifying one first chromosome gives a difference from the distribution parameters of the whole sequence 
of less than 2% (for a reliable number of fragments, exceeding thousand), 
but significantly speeds up processing (e.g. about 8 times for the mouse genome).

`-D|--dist <N,LN,G>`<br>
specifies the desired distribution type to call: `N` – normal, `LN` – lognormal, `G` – gamma.<br>
Сan be assigned independently of each other in any order. Characters are case-insensitive.<br>
For each specified type, the called distribution parameters are displayed, as well as the Pearson correlation coefficient 
(PCC) with the original sequence. The coefficient is calculated on the basis from the beginning of the distribution 
to the first point with an ordinate that is less than 0.1% of the maximum.<br>
The types of distributions are sorted by PCC in descending order, and the ratio of the PCC to the maximum, in percent, 
is indicated as well.<br>
While a lognormal type is specified (default or explicit), 
the original sequence is also automatically checked for normal distribution. 
If normal PCC is no less than lognormal PCC minus 2%, its parameters are also printed. 
This is done because for certain parameters these distribution may differ very slightly.<br>
Default: `LN` for BAM/BED (assuming fragments), `N` for FASTQ (assuming reads)

`-d|--dup <OFF|ON>`<br>
rejects/allows duplicate fragments/reads.<br>
This option is topical for BAM/BED files only.<br>
Default: `ON`

`-p|--pr-dist`<br>
prints original (actual) fragment/read length frequency distribution as a set of \<size\>-\<frequency\> pairs.<br>
This allows to visualize the distribution using some suitable tool such as Excel, etc.<br>
Printing is performed only to a file that duplicates the standard output (see `-O|--out` option).<br>
If duplicating output is not set, it is activated automatically.<br>
Input file with *.dist* extention ignores this option, but not the explicit `-O|--out` option.

`-s|--stats`<br>
prints input item issues statistics

`-O|--out [<name>]`<br>
duplicates standard output to specified file (except alarm messages).<br>
If <name> is not specified, duplicates output to file with name, 
constructed as input file short name with extention *.dist*.<br>
If the <name> denotes an existing folder, the output file is created in it according to the rule described above.<br>
If the input file itself has a *.dist* extention, then the "_out" suffix is added to the name.<br>
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
Output:
  -O|--out [<name>]     duplicate standard output to specified file
                        or to <in-file>.output.txt if <name> is not specified
  -T|--sep              use 1000 separator in output
  -V|--verbose <TOT|LAC|DET>
                        set output verbose level:
                        TOT - only total detailed,
                        LAC - laconic for each chromosome and total detailed,
                        DET - detailed for each chromosome [LAC]
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help WARNING     print usage information and exit
```
#### Input
Aligned DNA single- or paired end sequence (with fixed read length) 
in [BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) 
or [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.<br>
The program recognizes the file format automatically by their extention (case-insensitive).<br>
The sequence must be sorted.

#### Output
An example of the output(here $G stores the reference genome directory):
```
$ vAlign -to -g $G mInp-rqLow.B1.bam
mInp-rqLow.B1.bam
chrom 1
mismCnt	readCnt
-----------------
precise	872019
0	9301
1	3
2	2
...
49	81
50	102
-----------------
reads total per chrom 1:              1889303  (including 2834 (0.150%) duplicates)
reads mapped  to the correct chrom 1: 1754353  92.86%
from wich:
      mapped without mismathes:        881320  46.65%
      mapped with mismathes:           873033  46.21%
reads mapped to different chroms:      134950  7.143%
chrom 2
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
treats specified chromosome only.<br>
`name` identifies chromosome by number or character, e.g. `10` or `X`. Character is case-insensitive.<br>
The indication of one chromosome reduces run time on 1.5-20 times but reduces statistical completeness.

`--min-scr <int>`<br>
specifies score threshold for treated reads. Reads with the score equal or less then specified will be ignored.<br>
Default: all reads are accepted.

`--char-case <OFF|ON>`<br>
turns off/on recognition of uppercase and lowercase characters in template and test as different.<br>
Default: `OFF`.

`-O|--out [<name>]`<br>
duplicates standard output to specified file (except alarm messages).<br>
If <name> is not specified, duplicates output to file with name, 
constructed as input file short name with addition of the suffix “.output.txt”.<br>
If the <name> denotes an existing folder, the output file is created in it according to the rule described above.<br>
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
  -O|--out [<name>]     duplicate standard output to specified file
                        or to <sequence>.output.txt if <name> is not specified
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit 
```

### Details

#### Input
DNA sequence in [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format.
-t -O outputSE\x3_2 -d 0 -V dbg -g $(SM) -R 40 -f 160 outputSE\G-3_SE_frag.wig
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
<br>
If you face to bugs, or have questions/suggestions, please do not hesitate to write on fedor.naumenko@gmail.com
