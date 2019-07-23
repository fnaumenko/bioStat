/************************************************************************************
readDens calculates density profile and precise mean density of aligned DNA sequences
inside and outside of a given regions. 
	
Copyright (C) 2018 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 15.06.2019
-------------------------
This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the	GNU General Public License for more details.
************************************************************************************/

#include "readDens.h"
#include <fstream>	// for dout

using namespace std;

const string Product::Title = "readDens";
const string Product::Version = "1.0";
const string Product::Descr = "read density profile calculator";

const char*	ProgParam = "<in-file>";	// program parameter tip
const string OutFileSuff = "_dens.txt";
const string HelpOutFile = "duplicate standard output to specified file\nor to "
	+ string(ProgParam) + OutFileSuff + " if file is not specified";

enum { oTREAT, oOUTPUT, oOTHER };	// oOTHER should be the last 

const char* Options::OptGroups [] = { "Treatment", "Output", "Other" };
const BYTE Options::GroupCount = sizeof(Options::OptGroups)/sizeof(char*);

// --info option: types of info notations
const char* infos [] = { "NM", "CNT", "STAT" };	// corresponds to eInfo; iNONE and iLAC are hidden

// { char, str, Signs (8: hidden), type, group, 
//	defVal (if NO_DEF then no default value printed),
//	minVal (if NO_VAL then value is prohibited), maxVal, strVal, descr, addDescr }
Options::Option Options::List [] = {
	{ 'g', "gen",	0,	tNAME,	oTREAT, NO_DEF, 0, 0, NULL,
	"reference genome library or chromosome sizes file", NULL },
	{ 'c', Chrom::Abbr,	0,	tNAME,	oTREAT, NO_DEF, 0, 0, NULL, "treat specified chromosome only", NULL },
	{ 'f', "fbed",	0,	tNAME,	oTREAT, NO_DEF,	0,	0, NULL,
	"'template' bed file which features define treated regions", NULL },
	{ 'e',"ext-len",0,	tINT,	oTREAT, 200, 0, 1e3, NULL,
	"length by which the features in the 'template' bed file\nwill be extended in both directions before treatment", NULL },
	{ HPH,"gap-len",0,	tINT,	oTREAT, 1000, 10, 1e5, NULL,
	"minimal length of undefined nucleotides region in genome\nwhich is declared as a gap. Ignored for the genome size file", NULL },
	{ 'd', "dupl",	8,	tENUM,	oTREAT, TRUE,	0, 2, (char*)Options::Booleans,
	"accept duplicate reads", NULL },
	//{ HPH, "diff-sz",	0,	tENUM,	oTREAT, FALSE,	0, 2, (char*)Options::Booleans,
	//"allow to ignore reads with different size", NULL },
	{ HPH, "min-scr",	0,	tINT,	oTREAT, NO_DEF, 0, 1000, NULL,
	"score threshold for treated reads", NULL },
	{ 's',"space",	0,	tINT,	oTREAT, 200, 2, 1e4, NULL,
	"resolution: span in bp by which reads will be counted\nto define a density", NULL },
	{ HPH, "serv",	0,	tNAME,	oTREAT, NO_DEF, 0, 0, NULL,	"folder to store service files", NULL },
	{ 'i', "info",	8,	tENUM,	oOUTPUT,Obj::iEXT, Obj::iNM, Obj::iSTAT, (char*)infos,
	"print information about file:\n? - name only, ? - number of reads, ? - statistics", NULL },
	{ 'w', "warn",	8,	tENUM,	oOUTPUT,FALSE,	NO_VAL, 2, NULL,
	"print read ambiguities, if they exist", NULL },
	{ 'W', "win-freq",0,tENUM,	oOUTPUT,FALSE,	NO_VAL, 2, NULL,
	"print windows frequency distribution", NULL },
	{ 'o', "out",	2,	tNAME,	oOUTPUT,NO_DEF,	0,	0, NULL, HelpOutFile.c_str(), NULL },
	{ 't', "time",	0,	tENUM,	oOTHER,	FALSE,	NO_VAL, 2, NULL, "print run time", NULL },
	{ 'v', Version,	0,	tVERS,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, "print program's version", NULL },
	{ HPH, "summ",	8,	tSUMM,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, "print program's summary", NULL },
	{ 'h', "help",	0,	tHELP,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, "print usage information", NULL }
};
const BYTE	Options::OptCount = sizeof(Options::List)/sizeof(Options::Option);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, ProgParam, true, "alignment in bam/bed format" }
};
const BYTE Options::UsageCount = sizeof(Options::Usages)/sizeof(Options::Usage);

dostream dout;		// stream's duplicator

/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv, ProgParam);
	if( fileInd < 0 )	return 1;		// wrong option or tip output
	int ret = 0;						// main() return code
	Features* templ = NULL;
	
	Chrom::SetCustomOption(oCHROM);
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		// check file names first of all
		const char* aName = FS::CheckedFileName(argv[fileInd]);	// input alignment
		const char* tName = FS::CheckedFileName(oFBED);			// template
		const char* sName = FS::CheckedDirName(oSERV);			// service dir

		dout.OpenFile(Options::GetSVal(oOUTFILE, aName, OutFileSuff));	
		if(!sName && FS::HasExt(aName, FT::Ext(FT::ABED)))
			Err("for BED file " + Options::OptionToStr(oGEN) + 
			" or " + Options::OptionToStr(oSERV) + " is required", aName).Throw();

		Obj::eInfo info	= Obj::eInfo(Options::GetIVal(oINFO));	// if print ambiguities info
		ChromReadDistrib::WinLen = Options::GetIVal(oSPACE);
		ChromReadDistrib::PrintFreq = Options::GetBVal(oFREQ);

		ChromSizes cSizes(Options::GetSVal(oGEN), Options::GetSVal(oSERV), true);
		Reads align(ProgParam, aName, cSizes, info, true, true, 
			Options::GetBVal(oALARM), Options::GetBVal(oDUPL), Options::GetIVal(oMINSCR));
		Read::Len = align.ReadLen();
		if( tName ) {
			templ = new Features(sTemplate, tName, cSizes, info, true, true, Options::GetBVal(oALARM));
			templ->Extend(Options::GetIVal(oEXTLEN), cSizes, info);
			templ->CheckFeaturesLength(Options::GetIVal(oSPACE), "space", "template");
			align.SetCommonChroms(*templ, false, true);
		}
		dout << EOL;

		DefRegions gRgns(cSizes, Options::GetIVal(oGAPLEN));

		GenReadDistrib(align, gRgns, templ).PrintDensity();
	}
	catch(Err &e)		{ ret = 1;	cerr << e.what() << EOL; }
	catch(exception &e)	{ ret = 1;	cerr << " exception " << e.what() << EOL; }

	if(templ)	delete templ;
	timer.Stop(true);
	return ret;
}

chrlen	ChromReadDistrib::WinLen;
bool	ChromReadDistrib::IsFG;
bool	ChromReadDistrib::PrintFreq;	// true if windows frequency should be printed
