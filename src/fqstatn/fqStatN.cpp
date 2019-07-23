 /*
fqStatN outputs statistic of ambiguous reference characters ‘N’ in fq-file.

Copyright (C) 2018 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 20.06.2019
-------------------------

This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the	GNU General Public License for more details.

Sample output:

'N' POSITION STATISTICS
pos     count   % of total 'N'
-----------------------------
	1      29588   71.1%
	4      4314    10.4%
21      69      0.166%
...
READ TEMPLATE STATISTICS
position  10        20        30        40          patterns count
12345678901234567890123456789012345678901234567890
----------------------------------------------------------------------
N.................................................    29566     0.172%
...N..............................................     4313     0.0251%
...............................N..NN.............N      213     0.00124%
...
'N' relative to the total number of nucleotides: 0.00483%
Reads that include 'N' relative to the total number of reads: 0.211%
 */

#include "fqStatN.h"
#include <fstream>
#include <algorithm>    // std::sort

using namespace std;

const string Product::Title = "fqStatN";
const string Product::Version = "1.0";
const string Product::Descr = "fastq 'N' statistics calculator";

const char*	ProgParam = "<sequence>";	// program parameter tip
const string OutFileSuff = "_statn.txt";
const string HelpOutFile = "duplicate standard output to specified file\nor to "
	+ string(ProgParam) + OutFileSuff + " if file is not specified";

enum eOptGroup	{ oOPTION };	// oOTHER should be the last
const char* Options::OptGroups [] = { NULL };
const BYTE Options::GroupCount = sizeof(Options::OptGroups)/sizeof(char*);

// { char, str, Signs (8: hidden), type, group, 
//	defVal (if NO_DEF then no default value printed),
//	minVal (if NO_VAL then value is prohibited), maxVal, strVal, descr, addDescr }
Options::Option Options::List [] = {
	{ 'o', "out",	2,	tNAME,	oOPTION, NO_DEF,	0,	0, NULL, HelpOutFile.c_str(), NULL },
	{ 't', "time",	0,	tENUM,	oOPTION, FALSE,	vUNDEF,	2, NULL, "print run time", NULL },
	{ 'v', Version,	0,	tVERS,	oOPTION, NO_DEF, NO_VAL,0, NULL, "print program's version", NULL },
	{ HPH, "summ",	8,	tSUMM,	oOPTION, vUNDEF, vUNDEF,0, NULL, "print program's summary", NULL },
	{ 'h', "help",	0,	tHELP,	oOPTION, vUNDEF, vUNDEF,0, NULL, "print usage information", NULL }
};
const BYTE	Options::OptCount = sizeof(Options::List)/sizeof(Options::Option);

const Options::Usage Options::Usages[] = {
	{ vUNDEF, ProgParam, true, "nucleotide sequence in fastq format" }
};
const BYTE Options::UsageCount = sizeof(Options::Usages)/sizeof(Options::Usage);

dostream dout;	// stream's duplicator

int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv, ProgParam);
	if( fileInd < 0 )	return 1;		// wrong option
	int ret = 0;						// main() return code

	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* sName = FS::CheckedFileName(argv[fileInd]);
		dout.OpenFile(Options::GetSVal(oOUTFILE, sName, OutFileSuff));	
		dout << sName << SepCl;	fflush(stdout);
		if(!FS::HasExt(sName, FT::Ext(FT::FQ))) Err("wrong format").Throw();
		FqFile fq(sName);

		StatN::Scan(fq);
	}
	catch(const Err &e)			{ ret = 1;	cerr << e.what() << EOL; }
	catch(const exception &e)	{ ret = 1;	cerr << e.what() << EOL; }
	timer.Stop("wall-clock: ", false, true);
	return ret;
}

/************************  class StatN ************************/

void StatN::Scan	(FqFile & fq)
{
	UINT i, n;	// n used as N counter in Read in 1st part (counting), and as counter in 2th (output)
	ULONG cntTotalN = 0, cntTotalReads = 0;
	readlen k;
	bool insert;
	const char *read;
	vector<TemplN> templs;

	templs.reserve(20);
	fq.GetSequence();
	readlen rLen = fq.ReadLength();
	Array<char>	buf(rLen+1);
	Array<ULONG> distr(rLen);	// array of 'N' frequencies
	//vector<char>	buf(rLen+1);
	//vector<ULONG> distr(rLen);	// array of 'N' frequencies

	ULONG cnt = 0;
	// GET OCCURENCES
	do {
		read = fq.GetCurrRead();
		buf.Clear();
		for(n=0, i=0; i<rLen; i++)
			if( read[i] == cN ) {
				buf[n++] = i;
				distr[i]++;
			}
		// insert new line in statistic
		if(n > 0) {
			insert = false;
			for(i=0; i<templs.size(); i++)
				if( templs[i].Count == n
				&& !strcmp(buf.Data(), templs[i].Pos.Data()) ) {
					templs[i].CountRead++;
					cntTotalReads++;
					insert = true;
					break;
				}
			if( !insert )
				templs.push_back(TemplN(n, buf));
			cntTotalN += n;
		}
	}
	while( fq.GetSequence() );

	// OUTPUT RESULT
	dout << "total " << fq.Count() << " reads length of " << rLen << EOL;;
	if( templs.size() ) {
		dout<< "'N' POSITION STATISTICS\n"
			<< "pos\tcount\t% of total 'N'\n"
			<< "------------------------------\n";
		for(k=0; k<rLen; k++)
			if( distr[k] )
				dout << setw(2) << int(k) << TAB << distr[k] << TAB
						<< sPercent(Percent(distr[k], cntTotalN), 3, 0, false) << EOL;

		dout << "\nREAD PATTERN STATISTICS\n";
		//sort(templs.begin(), templs.end(), Compare);
		// OUTPUT HEADER
		dout << "position";
		// output tens
		for(k=1, n=7; k<=rLen/10; n=0, k++) {
			for(; n<8; n++) 	dout << ' ';
			dout << k*10;
		}
		dout << setfill(BLANK) << "\tcount\n";
		// output units
		for(k=0; k<rLen/10; k++) {
			for(n=1; n<10; n++)		dout << n;	// !!!
			dout << 0;
		}
		// output rest of units
		k = rLen%10;
		for(n=0; n<k; n++)		dout << n;
		dout << EOL;
		for(k=0; k<rLen+2*8+4; k++)	dout << HPH;
		dout << EOL;

		// OUTPUT ENTRIES
		for(i=0; i<templs.size(); i++) {
			for(n=k=0; k<rLen; k++)
				if(templs[i].Pos[n]-1 == k) {	dout << cN; n++; }
				else							dout << DOT;
			dout << setw(9) << templs[i].CountRead;
			dout << TAB << sPercent(Percent(templs[i].CountRead, fq.Count()), 3, 0, false) << EOL;
		}
	
		dout << "\n'N' relative to the total number of nucleotides: " 
			<< sPercent(Percent(cntTotalN, fq.Count()*rLen), 3, 0, false)  << EOL;
		dout << "Reads that include 'N' relative to the total number of reads: " 
			<< sPercent(Percent(cntTotalReads, fq.Count()), 3, 0, false)  << EOL;
	}
	else
		dout << "No reads included 'N'\n";
}

/************************  end of class StatN ************************/

