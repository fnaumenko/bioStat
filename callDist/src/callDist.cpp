/************************************************************************************
callDist is designed to calculate the size distribution parameters of
paired-end fragments OR reads

The frequency distribution is displayedin the form of lines,
each of which contains one pair <fragment length><TAB><frequency>.

Copyright (C) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 07/28/2024
-------------------------
************************************************************************************/

#include "callDist.h"

const string Product::Title = "callDist";
const string Product::Version = "2.0";
const string Product::Descr = "PE-fragment-size/read-length distribution parameters caller";

const char* ProgParam = "<in-file>";	// program parameter tip
const string OutFileExt = FT::Ext(FT::eType::DIST);

const char* inputs[] = { "FRAG","READ" };	// input option; corresponds to Inp
const char* dTypes[] = { "N","LN","G" };	// input distrib option; corresponds to InType	

// *** Options definition

enum eOptGroup { gTREAT, gOUTPUT };
const char* Options::OptGroups[] = { "Processing", "Output" };
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

const BYTE Options::Option::IndentInTabs = 3;

// { char, str, Signs (8: hidden), type, group, 
//	defVal (if NO_DEF then no default value printed),
//	minVal (if NO_VAL then value is prohibited), maxVal, strVal, descr, addDescr }
Options::Option Options::List[] = {
	{ 'i', "inp",	tOpt::NONE,	 tENUM,	gTREAT, float(InpType::FRAG), float(InpType::FRAG), ArrCnt(inputs), (char*)inputs,
	"input data to call distribution: ? - fragments, ? - reads", NULL },
	{ 'c', sChrom,	tOpt::NONE,tNAME,	gTREAT, NO_DEF, 0, 0, NULL,	sHelpChrom, NULL },
	{ 'D',"dist",	tOpt::NONE,	 tCOMB,	gTREAT, Distrib::LNORM, Distrib::NORM, ArrCnt(dTypes), (char*)dTypes,
	"called distribution (can be combined in any order):\n? - normal, ? - lognormal, ? - Gamma", NULL },
	{ 'd', "dup",	tOpt::NONE,	 tENUM,	gTREAT,	TRUE,	0, 2, (char*)Booleans, "allow duplicates", NULL },
	{ 'p', "pr-dist",tOpt::NONE, tENUM,	gOUTPUT,FALSE,	NO_VAL, 0, NULL, "print obtained frequency distribution", NULL },
	{ 's', "stats",	tOpt::NONE,	 tENUM,	gOUTPUT,FALSE,	NO_VAL, 0, NULL, "print input item issues statistics", NULL },
	{ 'O', sOutput,	tOpt::FACULT,tNAME,	gOUTPUT,NO_DEF,	0,	0, NULL, DoutHelp(ProgParam, OutFileExt), NULL },
	{ 't',	sTime,	tOpt::NONE,	 tENUM,	gOUTPUT,FALSE,	NO_VAL, 0, NULL, sHelpTime, NULL },
	{ HPH,	sSumm,	tOpt::HIDDEN,tSUMM,	gOUTPUT,NO_DEF, NO_VAL, 0, NULL, sHelpSummary, NULL },
	{ 'v',	sVers,	tOpt::NONE,	 tVERS,	gOUTPUT,NO_DEF, NO_VAL, 0, NULL, sHelpVersion, NULL },
	{ 'h',	sHelp,	tOpt::NONE,	 tHELP,	gOUTPUT,NO_DEF, NO_VAL, 0, NULL, sHelpUsage, NULL }
};
const BYTE Options::OptCount = ArrCnt(Options::List);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, ProgParam, true, "paired-end alignment in bam/bed format OR reads in fq/bam/bed format" }
};
const BYTE Options::UsageCount = ArrCnt(Options::Usages);

dostream dout;	// stream's duplicator

// Returns explicitly defined by user combo type, otherwise default combo type
Distrib::eCType GetType(Distrib::eCType defType)
{
	return Options::Assigned(oDTYPE) ? Distrib::eCType(Options::GetIVal(oDTYPE)) : defType;
}

/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv, ProgParam);
	if (fileInd < 0)	return 1;		// wrong option or tip output
	int ret = 0;						// main() return code

	Chrom::SetUserChrom(Options::GetSVal(oCHROM));
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* iName = FS::CheckedFileName(argv[fileInd]);		// input
		const bool prDist = Options::GetBVal(oPR_DIST);
		const FT::eType ftype = FT::GetType(iName);

		Options::SetDoutFile(oDOUT_FILE, iName, prDist && ftype != FT::eType::DIST, OutFileExt);

		// take distribution
		const bool fragType = InpType(Options::GetIVal(oINPUT)) == InpType::FRAG;
		const string sWFormat = "wrong format";
		const string sRExt = "should be BED or BAM";

		switch (ftype) {
		case FT::eType::BED:
		case FT::eType::BAM:
			if (fragType)
				FragDist(iName, Options::GetBVal(oPR_STATS)).Print(GetType(Distrib::LNORM), prDist);
			else
				ReadDist(iName, Options::GetBVal(oPR_STATS)).Print(GetType(Distrib::NORM), prDist);
			break;
		case FT::eType::FQ:
			dout << iName;	cout.flush();
			if (fragType && Options::Assigned(oINPUT))
				Err(sWFormat + " for fragment distribution; " + sRExt).Throw();
			FqReadDist(iName).Print(GetType(Distrib::NORM), prDist);
			break;
		case FT::eType::DIST:
			dout << iName;
			Distrib(iName, dout).Print(dout, Distrib::eCType(Options::GetIVal(oDTYPE)), true, false);
			break;
		default:
			Err(sWFormat + SepSCl + sRExt + " or FQ", iName).Throw();
			break;
		}
	}
	catch (const Err& e) { ret = 1; cerr << e.what() << endl; }
	catch (const exception& e) { ret = 1; cerr << e.what() << endl; }
	timer.Stop();
	return ret;
}
