/************************************************************************************
fragDist - fragment's distribution - is designed to calculate the frequency distribution of paired-end fragments.
 
The frequency distribution is displayedin the form of lines, 
each of which contains one pair <fragment length><TAB><frequency>.
Then this output is transferred to the Excel, which allows to plot it.
 
Copyright (C) 2018 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 20.06.2019
-------------------------
 
This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See GNU General Public License for more details.
************************************************************************************/

#include <iostream>
#include "fragDist.h"
#include <fstream>	// for dout

using namespace std;

const string Product::Title = "fragDist";
const string Product::Version = "1.0";
const string Product::Descr = "paired-end fragment size distribution parameters and frequency profile calculator";

const char*	ProgParam = "<in-file>";	// program parameter tip
const string OutFileExt = FT::Ext(FT::DIST);
const string HelpOutFile = "duplicate standard output to specified file\nor to "
	+ string(ProgParam) + OutFileExt + " if file is not specified";

enum eOptGroup	{ oOTHER };	// oOTHER should be the last 
const char* Options::OptGroups [] = { NULL };
const BYTE Options::GroupCount = sizeof(Options::OptGroups)/sizeof(char*);

// { char, str, Signs (8: hidden), type, group, 
//	defVal (if NO_DEF then no default value printed),
//	minVal (if NO_VAL then value is prohibited), maxVal, strVal, descr, addDescr }
Options::Option Options::List [] = {
	{ 'c',Chrom::Abbr,0,tNAME, oOTHER,	NO_DEF, 0, 0, NULL, "treat specified chromosome only", NULL },
	{ 'd', "no-dup",8,	tENUM,	oOTHER,	FALSE,	NO_VAL, 0, NULL, "reject duplicate reads", NULL },
	{ 'i', "info",	8,	tENUM, oOTHER,	FALSE,	NO_VAL,	0, NULL,
	"print statistical info about alignment", NULL },
	{ 'w', "warn",	8,	tENUM,	oOTHER, FALSE,	NO_VAL, 0, NULL,
	"print each read ambiguity, if they exist" },
	{ 'D', "dist",	0,	tENUM,	oOTHER,	FALSE,	NO_VAL, 0, NULL,
	"print fragment size frequency distribution", NULL },
	{ 'o', "out",	2,	tNAME,	oOTHER,	NO_DEF,	0,	0, NULL, HelpOutFile.c_str(), NULL },
	{ 't', "time",	0,	tENUM,	oOTHER,	FALSE,	NO_VAL, 0, NULL, "print run time", NULL },
	{ 'v', Version,	0,	tVERS,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, "print program's version", NULL },
	{ HPH, "summ",	8,	tSUMM,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, "print program's summary", NULL },
	{ 'h', "help",	0,	tHELP,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, "print usage information", NULL }
};
const BYTE	Options::OptCount = sizeof(Options::List)/sizeof(Options::Option);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, ProgParam, true, "paired-end alignment in bam/bed format" }
};
const BYTE Options::UsageCount = sizeof(Options::Usages)/sizeof(Options::Usage);

dostream dout;	// stream's duplicator

/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv, ProgParam);
	if( fileInd < 0 )	return 1;		// wrong option or tip output
	int ret = 0;						// main() return code
	
	Chrom::SetCustomOption(oCHROM, true);
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* iName = FS::CheckedFileName(argv[fileInd]);		// input
		bool isDist = FT::GetType(iName)==FT::DIST;
		if( dout.OpenFile(Options::GetSVal(oOUTFILE, iName, OutFileExt))
		&& isDist && !Options::GetSVal(oOUTFILE))
			Err("if input file is distribution, the output name should be different").Throw();

		if(isDist)		FragFreq(iName).Print(dout, Options::GetBVal(oPR_DIST));
		else {
			ChromSizes cSizes;
			Reads reads(ProgParam, iName, cSizes, 
				Obj::eInfo (Options::GetBVal(oINFO) ? Obj::iSTAT : Obj::iEXT), 
				true, false, Options::GetBVal(oALARM),
				true, !Options::GetBVal(oNO_DUPL), -1);

			FragDist(reads).PrintDist( Options::GetBVal(oPR_DIST) );
		}
	}
	catch(Err &e)				{ ret = 1; cerr << e.what() << endl; }
	catch(const exception &e)	{ ret = 1; cerr << e.what() << EOL; }
	timer.Stop("wall-clock: ", false, true);
	return ret;
}


FragDist::FragDist(Reads& test)
{
	Reads::cItemsIter rit, end;
	Reads::cIter cit=test.cBegin();

	//_freq.reserve(1000);
	_waits.rehash( test.Count(CID(cit)) / 2 );

	Read::Len = test.ReadLen();
	for(; cit!=test.cEnd(); cit++) {
		//if(!Chrom::StatedAll && Chrom::StatedID() != CID(cit))	continue;
		//cout << Chrom::AbbrName(CID(cit)) << "...";
		//fflush(stdout);
		end = test.ReadsEnd(cit);
		for(rit=test.ReadsBegin(cit); rit!=end; rit++)
			Add(rit);
		Clear();
		//cout << Done << EOL;
	}
}
