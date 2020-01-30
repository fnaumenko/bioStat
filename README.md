# bioStat
Statistical package for NGS data.<br>
It includes next utilities:<br>
 * [**cc**](#biocc)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; advanced correlation calculator for basic bioinformatics file formats<br>
 * [**fragdist**](#fragdist)&nbsp;&nbsp;&nbsp; calls fragment size distribution<br>
 * [**readdens**](#readdens)&nbsp; read density profile calculator<br>
 * [**valign**](#valign)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; alignment verifier<br>
 * [**fqstatn**](#fqstatn)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; fastq 'N' statistics calculator

## Usage
`biostat <command> [options] [<file>…]`<br>
or<br>
`<Command> [options] [<file>…]`<br>
In the second case, the output to the terminal is carried out immediately, which can be informative when processing large input files.

## Installation
### Executable file

[Download Linux version](https://github.com/fnaumenko/biostat/releases/download/1.0/biostat-Linux-x64.tar.gz) 
and extract **biostat** folder typing `tar -xf biostat-Linux-x64.tar.gz`.<br>
[Download Windows version](https://github.com/fnaumenko/biostat/releases/download/1.0/biostat-Windows-x64.zip) 
and extract **biostat** folder using [WinRAR](https://www.win-rar.com/download.html?&L=0) or another ZIP-archiver.

Alternative download in Linux:<br>
`wget https://github.com/fnaumenko/bioStat/releases/download/1.0/biostat-Linux-x64.tar.gz`<br>

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
fast advanced **C**orrelation **C**alculator for basic **bio**informatics file formats.<br>
It computes Pearson’s and signal’s correlation coefficients for densities, coverage and features.<br>
Program allows to obtain correlation coefficients for the whole genome, for each chromosome separately 
and for predefined regions within the chromosomes. 
It can print coefficients for each predefined region and coefficients frequency histogram.<br>
**bioCC** is designed to treat a bunch of files at once.

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
  -a|--align            input bed files are alignments
  -g|--gen <name>       chromosome sizes file
  -l|--list <name>      list of multiple input files.
                        First (primary) file in list is comparing with others (secondary)
Processing:
  -c|--chr <name>       treat specified chromosome only
  -r|--cc <P,S>         correlation coefficient, in any combination: P - Pearson, S - signal [P]
  -s|--space <int>      resolution: span in bp by which reads will be counted to define a density.
                        For the alignments only [100]
Region processing:
  -f|--fbed <name>      'template' ordinary bed file which features define compared regions.
                        Ignored for the ordinary beds
  -e|--ext-len <int>    length by which the features in primary file (for ordinary beds) or in 'template'
                        (for alignments and wigs) will be extended in both directions before treatment [0]
  --ext-step <int>      step of extending features in primary bed file; if 0 then no step calculation.
                        For the ordinary beds only [0]
  --norm <OFF|ON>       normalize regions before calculation. Ignored for the ordinary beds [ON]
Output:
  -С|--pr-cc <IND,TOT>  print coefficient, in any combination:
                        IND - for each chromosome individually, TOT - total [IND]
  -B|--bin-width <float>
                        print frequency histogram with given bin width [0]
  -S|--sort <RGN|CC>    print region coefficients, sorted by: RGN - regions, CC - coefficients
  -i|--info <LAC|NM|CNT|STAT>
                        print information about file:
                        LAC - laconic, NM - name only, CNT - number of items, STAT - statistics [NM]
  -o|--out              duplicate standard output to bioCC_out.txt file
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```

### Details

#### Input data types
Alignment read densities should be presented by [BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) 
or [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.<br>
Sequencing coverage should be presented by [WIG](https://genome.ucsc.edu/goldenpath/help/wiggle.html) format.<br>
Features are also presented by [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.

**BED**<br>
Formally bed file as a result of peak calling process, and bed file represents alignment have the common required fields, but they differ in the interpretation. 
For greater certainty, the file of the first type will be called *ordinary* bed, the file of the second type will be called *alignment* bed.<br>
For more information about read’s density see option `-s|--space`.

**WIG**<br>
*bedGraph* type is not supported in this version, and the only *wiggle* in *variableStep* format is acceptable.<br>
Note, that one of the fastest *wiggle* generator [PeakRanger](http://ranger.sourceforge.net/manual1.18.html) 
(it also supports the separate strand generation) produces data with some peculiarity, 
namely each its interval has initial span = 1. 
This does not allow to compare coverage.<br>
To level this feature, use the [wigReg](https://github.com/fnaumenko/wigReg) software.<br>
**wigReg** also allows to shorten *wig* files generated by [MACS](http://liulab.dfci.harvard.edu/MACS/00README.html) 2-10 times.

All input files should have chromosome's items (features/reads/intervals) been clustering together.<br>
For faster processing, items belonging to the same chromosome also should be sorted in position ascending order.<br>
*Wig* files usually meet these requirements. 
Opposite, *bed* files often have messy initialization. 
The simplest way to prepare *bed* files is to sort them, for example by using **sortBed** utility from [bedtools](http://bedtools.readthedocs.io/en/latest/) package.<br>
It concerns 'template' *bed* file as well (see `-f|--fbed` option).

**bioCC** recognizes the file formats automatically by their extension, so the extensions should be BAM, BED or WIG (case-insensitive). 
To distinguish between *ordinary beds* and *alignments*, a special option `-a|--align` is provided. 

#### Input data order
Comparable files can be represented both through program parameters, and by means of a file containing their names. 
In both cases the first file in a list – *primary* – is compared to the others – *secondary* – by turns.
The *primary* file specifies a set of compared chromosomes (unless it is limited by the option `-c|--chr`). 
If one or more chromosomes are represented in the *primary* file, only they will be compared.<br>
Also, in the case of *wiggle*, the *primary* defines a resolution (see option ```-s|--space``` for more information).

Be careful by using standard naming conventions like *abc?.bed*, *a\*.bed*. 
You must be sure that first of the lexically recognized files is really primary, and that all other files really need to be compared.<br>
For more details about list of file names as a file see option ```-l|--list```.

### Options description
Note:<br>
Enumerable option values are case-insensitive.<br>
Compressed input files in gzip format are acceptable.

`-a|--align`<br>
indicates that input bed files are *alignments*, so density correlation would be performed.<br>
Since *alignment* may be considered as an 'instance' of *ordinary bed*, in theory both cases should give the same result. 
It is true while the *alignments* as well as the *ordinary beds* have no *ambiguous* reads/features (see ```–i|--info``` option for definition). 
But in practice this condition is almost never achieved, and these ambiguities are resolved by different way. 
Consequently the results can be dramatically different.<br>
In addition, *ordinary beds* are compared using a separate, ultra-fast algorithm.<br>
If alignments are treated without this option, bioCC will print a warning message.<br>
If *ordinary bed* files have features of different sizes and are treated with this option, 
**bioCC** will print a cancel message and complete.

`-g|--gen <name>`<br>
specifies chromosome sizes file. Required for BED and WIG files.<br>

`-l|--list <file>`<br>
specifies an external list of compared files. 
The list is a plain text file, in which each file name is located on a separate line.<br>
Lines beginning with ‘#’ are considered as comments and will be skipped.<br>
Empty lines are allowed.<br>
This option abolishes input files as parameters.

`-c|--chr <name>`<br>
treats specified chromosome only. 
The value `name` is the chromosome identifier; it is a number or character, for example, `10`, `X`.<br>
The indication of one chromosome reduces run time on 1.5-20 times depending on how far this chromosome is placed in an input *alignment* or *wig*.<br>
For *ordinary beds* it has no time-improvement effect: any result appears quickly.<br>

`-r|--cc <P,S>`<br>
specifies correlation coefficient, in any combination: `P` – Pearson, `S` – signal.<br>
For more details see [Pearson and signal correlation](#pearson-and-signal-correlation).<br>
Default: Pearson

`-s|--space <int>`<br>
specifies an alignment resolution: the length of windows (span) in bp by which reads will be counted to define a density.<br>
Read is considered belonging to span if it`s centre is placed within the span. 
Thus each read is counted once. 
Then, when using the reference genome from the input sequences, the reference undefined regions (gaps) are preliminarily removed from the compared sequences: 
the regions separated by a gap are merged.<br>
As a result, the program compares the actual read density distributions.<br>
This option is topical for the *alignments* only.<br>
*Wiggles* are already presented according to their resolution.
**bioCC** takes the resolution of *primary wiggle* as the basic for treatment. 
If *secondary wiggle* resolution is differ, it will be reduced to the base. 
For example, if the primary resolution is 10 and the secondary resolution is 1, then every 10 partitions in the secondary data will be merged so finally the secondary resolution will also be 10. 
Contrariwise, if secondary resolution is 100, then each partition of the secondary data will be divided into 10 with the same density.<br>
The best way is to use *wiggles* with the same resolutions.<br>
Range: 2-1000<br>
Default: 100

`-f|--fbed <file>`<br>
specifies 'template' *ordinary* bed file with features that defines compared regions.<br>
In practice the most typical case is comparing detected peaks across the distributions.<br>
If options `-B|--bin-width` and `--f|--fbed` are not defined, only total coefficients for the all regions are printed.<br>
This option is ignored for the *ordinary bed* files.

`-e|--ext-len <int>`<br>
specifies the value by which all features in a 'template' *bed* file or in a *primary ordinary bed* file should be stretched in both directions before comparison.<br>
If set, all the features from 'template' or *primary bed* file will be stretched before processing: 
*start* positions are decreased for this value, *end* positions are increased.<br>
If stretched features become intersected, they are joined.<br>
This option is mainly constructed for enriched region comparison 
while initial binding sites are represented by ‘template’ or *primary* features. 
It is only relevant in addition to the option `-f|--fbed`.<br>
Range: 0-1000<br>
Default: 0

`--ext-step <int>`<br>
If set, activates the mode of consecutive calculation of the coefficients for stretching features in *primary ordinary* bed file with the stated step. 
The maximum value of the extension is limited by `--e|--ext-len` option.<br>
This option is topical for *ordinary bed* files only.<br>
Range: 0-500<br>
Default: 0 (no step calculation)

`--norm <OFF|ON>`<br>
turns off/on regions normalization before calculation. <br>
Normalization means levelling distribution’s patterns by their maximal values across regions.<br>
In spite the fact that of each pair of regions will be always compared correctly, the common result for all the regions may been unfairly falling down. It occurs when patterns levels are appreciably differ from region to region.<br>
This option is assigned to eliminate such effect.<br>
This option is topical only with option `-f|--fbed`.<br>
Default: `ON`

`-C|--pr-cc <IND,TOT>`<br>
print coefficients, in any combination: `IND` - for each chromosome individually, `TOT` - total.<br>
If only one chromosome is specified, the total coefficient is not output as an identical.<br>
Default: `IND`.

`-B|--bin-width <float>`<br>
If set, forces to consolidate coefficients into bins, and print histogram values.<br>
Histogram shows count of coefficients (frequency) within some value ranges (bins). 
It is printed as a list of pairs *\<bin upper bound\>\<count in bin\>*. 
Negative coefficients have been turning to absolute during consolidation.<br>
This option defines the width of bin as a part of 1.<br>
Empty bins at the edges are not printed.<br>
For example, if –b value is 0.1, and all coefficients are placing in the range 0.65 to 0.85, only three bins [0.6-] 0.7, [0.7-] 0.8, [0.8-] 0.9 would be printed.<br>
This option is topical only with option `-f|--fbed`.<br>
Range: 0-1<br>
Default: 0 (no consolidation)

```--sort <RGN|CC>```<br>
If set, forces to output coefficients calculated for each region as a list of pairs *\<number-of-region\>\<coefficient\>*. 
`RGN` value prescribes list to be sorted by region’s number, `CC` – by coefficient.<br>
If both of the coefficients are declared, list is sorted by Pearson.<br>
First region number is 1.<br>
This option is topical only with option `-f|--fbed`.<br>

`-i|--info <LAC|NM|CNT|STAT>`<br>
outputs information about items (features/reads/intervals).<br>
`LAC`:  laconic output. This value minimizes the output as possible to remain clear. It is constructed mainly for using in batch file.<br>
`NM`:   brief output. Prints file names without any additional information.<br>
`CNT`:  prints file names and number of all and accepted items, if they are different.<br>
`STAT`: prints item ambiguities statistics, if they exist.<br>
When we are speaking about *bed* files, there are some issues which can be normal or may be not – that depends of file’s destination. 
For instance, duplicated and crossed features are normal for the *alignments*, but rather unusual for the *ordinary beds*. 
On the other hand, features with different length are typical for the *ordinary beds*, but they are a rare case for the *alignments*. 
All these situations we call a conditional term *ambiguities*.<br>
It concerns to 'template' *bed* as well.<br>
**bioCC** treats ambiguities the most natural way. 
For *ordinary bed* it merges crossed or adjacent features, and ignores submerging and duplicated features. 
For *alignment* it accepts duplicated reads by default (in fact, all the ambiguities are permitted for the alignment, except different length of reads).<br>
Thus, not all records present in the file can be accepted.<br>
In some circumstances you need to be aware of these issues. 
The `STAT` value provides the summary method. 
It forces to display number of all recognized certain type ambiguities, and appropriate treatment.<br>
Note, that in contrast to *bed* with their features/reads, the number of intervals does not match the number of lines in the *wiggle* file. 
All contiguous bases with the same data value are considered to be a single interval, although it may correspond to several file lines.<br>
Default: `NM`

`-o|--out`<br>
duplicates standard output to **bioCC_out.txt** file.<br>
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
Signal method is a bit more sensible, but it is sensitive to the mean amplitude (*DC offset* in terms of signal function), as it it we can see on (fig 2).
This means that the greater the background level of the compared distributions, the less relevant is the Signal method.<br>
For this reason, Pearson’s method is recommended for the distributions comparison.<br>
But if background’s level is considered part of the measure of similarity (for example, by comparing two replicas for noise level),
 in this case signal method would be preferable.

B.  Bed features.
Bed features can be accounted as discrete function accepted one of two values: 0 or 1. In that case signal method becomes inappropriate due to obvious reason: it counts intersections only. 
For example, by comparison two mutually complementary beds (while function1 have zero value every time when function2 have non-zero and vice versa), signal coefficient would be 0. 
Although the correct answer is -1.<br>
Thus, for the features only the Pearson method is correct.

---
## fragDist

Calculates paired-end fragment size lognormal distribution parameters and frequency profile.<br>
Examples of frequency profiles and recovered distributions of experimental datasets from NCBI database are shown 
in the ![figure](https://github.com/fnaumenko/bioStat/tree/master/pict/PEdistribs_medium.png).

### Usage

`biostat fragdist [options] <in-file>`<br>
or<br>
`fragDist [options] <in-file>`

### Options
```
  -c|--chr <name>       treat specified chromosome only
  -D|--dist             print fragment size frequency distribution
  -o|--out [<name>]     duplicate standard output to specified file
                        or to <in-file>.dist if file is not specified
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```
### Details

#### Input
See Input section [here](#input_).<br>
The program can also accept a file with a ready-made lognormal size frequency distribution, in order to determine its parameters. 
Distribution format is described in the comments to the option `-D|--dist`. This file should have *.dist* extension.

#### Output
Called (calculated) lognormal mean, sigma, Mode and Mean (expected size).<br>
Fragment size frequency distribution is printed optionally.

#### Options description

`-c|--chr <name>`<br>
treats specified chromosome only. The value name is the chromosome identifier; 
it is a number or character, for example, 10 or  X.<br>
The indication of one chromosome reduces run time on 1.5-20 times depending on how far this chromosome is placed in an alignment. 
However, this simultaneously reduces the result accuracy. 
Practically the best choice is the first chromosome, which provides the minimum time and highest accuracy relative to other chromosomes.

`-D|--dist`<br>
prints fragment size frequency distribution as a set of \<frequency\>-\<size\> pairs.

`-o|--out [<name>]`<br>
duplicates standard output to specified file. If file is not specified, duplicates output to file with name, 
constructed as input alignment short name (without path and extension) with addition of the extension .freq.<br>
If the option is last and its value is not set, the next (and the last) token is treated as a program parameter. 
It is an analogue of the **tee** Linux command and is constructed rather for the execution under Windows.

---
## readDens

Calculates alignment density profile and precise mean density into inside and outside given regions.<br>
Density profile is NOT a coverage. 
It means a set of frequencies of the observed equal parts of the sequence with the same density. 
The program splits each given region into non-overlapping equal parts (windows), 
and then counts the number of windows containing the same number of reads.<br>
'Precise' means that all the undefined regions in the reference genome are excluded from consideration. 
If the input regions are not defined, then only mean density is calculated for each chromosome.<br>
Example of density profiles of 8 experimental datasets from NCBI database are shown in the ![figure](https://github.com/fnaumenko/bioStat/tree/master/pict/readDensProfile_Oct4-Sox2.png).

### Usage
`biostat readdens [options] <in-file>`<br>
or<br>
`readDens [options] <in-file>`

### Options:
```
Treatment:
  -g|--gen <name>       reference genome library or chromosome sizes file
  -c|--chr <name>       treat specified chromosome only
  -f|--fbed <name>      'template' bed file which features define treated regions
  -e|--ext-len <int>    length by which the features in the 'template' bed file
                        will be extended in both directions before treatment [200]
  --gap-len <int>       minimal length of undefined nucleotides region in genome
                        which is declared as a gap. Ignored for the genome size file [1000]
  --min-scr <int>       score threshold for treated reads
  -s|--space <int>      resolution: span in bp by which reads will be counted
                        to define a density [200]
  --serv <name>         folder to store service files
Output:
  -W|--win-freq         print windows frequency distribution
  -o|--out [<name>]     duplicate standard output to specified file
                        or to <in-file>_dens.txt if file is not specified
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help WARNING     print usage information and exit
```

### Details

#### Input_
Aligned DNA sequence in [BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) 
or [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. 
Format is automatically recognized by file extension, so it should be BAM or BED (case-insensitive).<br>
Reads does not have to be sorted, but must be grouped by chromosomes. 
This means that reads belonging to the same chromosome must be arranged sequentially, in a single group. 
The simplest way to ensure this is to pre-sort the alignment.<br>
BAM files are read 2-4 times slower than even zipped BED files.

#### Output
Mean density is measured in read per kilobase.<br>
Density profile is printed as a set of a pairs \<number of read in window\>\<count of window\>.<br>
The results are calculated for each chromosome separately, though the total mean density is also printed.

#### Options description
`-g|--gen <file>`<br>
Chromosome sizes file or reference genome library (directory containing nucleotide chromosome sequences in FASTA format).<br>
If chromosome sizes file is specified, the density will be considered over the entire chromosome length.<br>
If reference library is specified, all the undefined regions in the reference genome will be excluded from processing, 
which means counting the actual density. 
Undefined regions (gaps) are regions entirely composed of the ambiguous reference code N. 
The minimal length of accounting gaps is managed by `--gap-len` option.<br>
For example, chromosome 1 from mm9 library contains 14 regions, separated by gaps with length more than 400 bps, 
and 10 regions, separated by gaps with length more than 1000.<br>
Gaps are determined by scanning the reference genome once and are saved in service files with ‘region’ extension 
in the same directory by default. To assign another service directory, use `--serv` option. 
If this option is not set, and the reference directory is write-protected, scanning will occur every time.

The reference library may contain .fa and/or .fa.gz files. First program searches for unzipped files in it. 
If there are no such files, or the file corresponded to chromosome specified by option `–c|--chr` is absent, 
the program searches for zipped files.
One can obtain a genome library in UCSC: ftp://hgdownload.soe.ucsc.edu/goldenPath/ or in Ensemble: ftp://ftp.ensembl.org/pub/release-73/fasta storage. 
In the second case please copy genomic sequences with the same masked type only, e.g. unmasked (‘dna'), since program does not recognise mask’s types.

This option or the `--serv` option replacing it is mandatory for the BED file. 
For an BAM file, it can be omitted, which is equivalent to specifying chromosome sizes file (since the BAM/SAM file contains the chromosome sizes itself).

`-c|--chr <name>`<br>
treats specified chromosome only. The value `name` is the chromosome identifier; 
it is a number or character, for example, 10 or X.
The indication of one chromosome reduces run time on 1.5-20 times depending on how far this chromosome is placed in an alignment. 
However, this simultaneously reduces the result accuracy.

`-f|--fbed <name>`<br>
specifies 'template' ordinary BED file with features that defines treated regions. 
The density profile will be constructed on the region resulting from the merger of all regions specified by the 'template' regions.<br>
This option abolishes the merge of defined regions specified by `-g|--gen` option.

`-e|--ext-len <int>`<br>
specifies value by which all features in 'template' should be stretched in both directions before the density count.<br>
If set, all the features from 'template' bed file will be stretched before processing: 
*start* positions are decreased for this value, *end* positions are increased.<br>
If stretched features become intersected, they are joined.<br>
This option is mainly constructed for enriched region comparison while initial binding sites are represented by ‘template’ features. 
It is only relevant in addition to the option `-f|--fbed`.<br>
Range: 0-1000<br>
Default: 0

`--gap-len <int>`<br>
Minimal length of undefined nucleotides region which is taken as a gap.<br>
Ignored for chromosome sizes file, specified by `-g|--gen` option or embedded in BAM file.<br>
Range: 10-1e5<br>
Default: 1000

`--min-scr <int>`<br>
specifies score threshold for treated reads. Reads with the score equal or less then stated will be ignored.<br>
Range: 0-1000<br>
Default: all reads are accepted

`-s|--space <int>`<br>
Resolution: the length of windows (span) in bp by which reads will be counted to define a density.<br>
Read is considered belonging to span if it`s centre is placed within the span. 
Thus each read is counted once. Then, while using the reference genome from the input sequences, 
and if 'template' is not specified, the reference undefined regions (gaps) are preliminarily removed from the compared sequences. 
This means that the regions separated by a gap are merged.<br>
As a result, the program compares the actual read density distributions.<br>
Range: 2-1e4<br>
Default: 200

`--serv <name>`<br>
specifies the service directory – a place for keeping service files *chr\<x\>.region* and chromosome sizes file. 
These files are created once by scanning reference genome, and after that are using as the gap data source. 
By default they are keeping in the reference genome folder, but if this folder is closed for writing, 
or you want to store these files separately for your own reasons, that is the place.<br>
Once the service files have been created, you can use this option instead of `-g|--gen`.<br>
Note that in order to avoid ambiguity, each reference genome library should have its own service directory.<br>
If this option is not set, and the reference genome folder is write-protected, scanning will occur every time.

`-o|--out`<br>
duplicates standard output to specified file. If file is not specified, duplicates output to file with name, 
constructed as input alignment short name (without path and extension) with addition of the suffix “_dens.txt”.<br>
If the option is last and its value is not set, the next (and the last) token is treated as a program parameter.<br>
It is an analogue of the **tee** Linux command and is constructed rather for the execution under Windows.

---
## vAlign
**V**erify **Align**ment is a fast verifier of reads forming the aligned DNA sequence, 
which is recalled from an artificial [FastQ](https://en.wikipedia.org/wiki/FASTQ_format) sequence. 
It compares the original and actual coordinates of each read and prints statistics of right and wrong mappings.

To do this each read in an initial artificial sequence should keep its location as a part of its name. 
Such template sequence can be generated by [**isChIP**](https://github.com/fnaumenko/isChIP) software. 
Both single end and paired end reads are possible. 
If initial reads do not answer this requirement, the treatment will be rejected.

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
See Input section [here](#input_).

#### Output
**vAlign** outputs number of exactly matched reads, and number of wrong placed reads with 0, 1, 2, … N mismatches, where N is length of read.

#### Options description
`-g|--gen <name>`<br>
specifies reference genome library or single nucleotide sequence.<br>
Genome library is a directory containing nucleotide sequences for each chromosome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.<br>
If `name` is a .fa[.gz] file, **vALign** accepts the corresponding chromosome as the only treated.<br>
Otherwise first the program searches for .fa files in the directory `name`. If there are no such files in this directory, **vALign** searches for .fa.gz files.<br>
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
duplicates standard output to specified file. If file is not specified, duplicates output to file with name, 
constructed as input alignment short name (without path and extension) with addition of the suffix “_valign.txt”.<br>
If the option is last and its value is not set, the next (and the last) token is treated as a program parameter.<br>
It is an analogue of the **tee** Linux command and is constructed rather for the execution under Windows.

---
## fqStatN

Calculates the statistics of occurrence of ambiguous code N, and patterns of reads including N 
in the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file.<br>
These statistics helps to better evaluate the quality of the sequencer output.

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
. . .
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
. . .
...............................N..NN..............        3     <0.001%
...............................N..N...............        2     <0.001%
.........................N.....N..NN.N...........N        1     <0.001%

'N' relative to the total number of nucleotides: 0.00483%
Reads that include 'N' relative to the total number of reads: 0.211%
```

---
If you face to bugs, incorrect English, or have commentary/suggestions, please do not hesitate to write me on fedor.naumenko@gmail.com
