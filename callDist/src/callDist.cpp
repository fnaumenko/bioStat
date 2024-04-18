/************************************************************************************
callDist is designed to calculate the size distribution parameters of
paired-end fragments OR reads

The frequency distribution is displayedin the form of lines,
each of which contains one pair <fragment length><TAB><frequency>.
Then this output is transferred to the Excel, which allows to plot it.

Copyright (C) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 11/10/2023
-------------------------

This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See GNU General Public License for more details.
************************************************************************************/

#include "callDist.h"

using namespace std;

const string Product::Title = "callDist";
const string Product::Version = "2.0";
const string Product::Descr = "PE-fragment-size/read-length distribution parameters caller";

const char* ProgParam = "<in-file>";	// program parameter tip
const string OutFileExt = FT::Ext(FT::eType::DIST);
const string HelpOutFile = sFileDuplBegin + string(ProgParam) + OutFileExt + sFileDuplEnd;

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
	{ 'i', "inp",	tOpt::NONE,	tENUM,	gTREAT, float(InpType::FRAG), float(InpType::FRAG), ArrCnt(inputs), (char*)inputs,
	"input data to call distribution: ? - fragments, ? - reads", NULL },
	{ 'D',"dist",	tOpt::NONE,	tCOMB,	gTREAT, Distrib::LNORM, Distrib::NORM, ArrCnt(dTypes), (char*)dTypes,
	"called distribution (in any order): ? - normal, ? - lognormal, ? - Gamma", NULL },
	{ 'd', "dup",	tOpt::NONE,	tENUM,	gTREAT,	FALSE,	0, 2, (char*)Options::Booleans, "allow duplicates", NULL },
	{ 'p', "pr-dist",tOpt::NONE,tENUM,	gOUTPUT,FALSE,	NO_VAL, 0, NULL, "print obtained frequency distribution", NULL },
	{ 's', "stats",	tOpt::NONE,	tENUM,	gOUTPUT,FALSE,	NO_VAL, 0, NULL, "print input item issues statistics", NULL },
	{ 'o', sOutput,	tOpt::FACULT,tNAME,	gOUTPUT,NO_DEF,	0,	0, NULL, HelpOutFile.c_str(), NULL },
	{ 't',	sTime,	tOpt::NONE,	tENUM,	gOUTPUT,FALSE,	NO_VAL, 0, NULL, sPrTime, NULL },
	{ HPH,	sSumm,	tOpt::HIDDEN,tSUMM,	gOUTPUT,NO_DEF, NO_VAL, 0, NULL, sPrSummary, NULL },
	{ 'v',	sVers,	tOpt::NONE,	tVERS,	gOUTPUT,NO_DEF, NO_VAL, 0, NULL, sPrVersion, NULL },
	{ 'h',	sHelp,	tOpt::NONE,	tHELP,	gOUTPUT,NO_DEF, NO_VAL, 0, NULL, sPrUsage, NULL }
};
const BYTE Options::OptCount = ArrCnt(Options::List);
//template <typename T, int N>
//static constexpr int array_size(T(&a)[N]) { return N; }

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, ProgParam, true, "paired-end alignment in bam/bed format OR reads in fq/bam/bed format" }
};
const BYTE Options::UsageCount = ArrCnt(Options::Usages);

dostream dout;	// stream's duplicator

// Returns explicitly defined by user combo type, otherwise default combo type
Distrib::eCType GetType(Distrib::eCType defType) {
	return Options::Assigned(oDTYPE) ? Distrib::eCType(Options::GetIVal(oDTYPE)) : defType;
}

void pr(float v) {
	if (v-int(v)!=0)	cout << fixed << setprecision(2);
	cout << v << LF;
}
/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv, ProgParam);
	if (fileInd < 0)	return 1;		// wrong option or tip output
	int ret = 0;						// main() return code

	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* iName = FS::CheckedFileName(argv[fileInd]);		// input
		const char* oName = Options::GetSVal(oOUTFILE);				// output
		const FT::eType ftype = FT::GetType(iName);

		// while input file is distribution, add "_out" to output file name -
		// only if isn't specified or the same as input
		if (Options::Assigned(oOUTFILE)
			&& !dout.OpenFile(Options::GetFileName(oOUTFILE, iName,
				(ftype == FT::eType::DIST &&
					(!oName || !strcmp(FS::FileNameWithoutExt(iName).c_str(), oName)) ?
					"_out" : strEmpty)
				+ OutFileExt))
			)
			Err(Err::FailOpenOFile).Throw();

		// take distribution
		const bool fragType = InpType(Options::GetIVal(oINPUT)) == InpType::FRAG;
		const bool prDist = Options::GetBVal(oPR_DIST);
		const bool prStat = Options::GetBVal(oPR_STATS);
		const string sWFormat = "wrong format";
		const string sRExt = "should be BED or BAM";

		dout << iName;	cout.flush();
		switch (ftype) {
		case FT::eType::BED:
		case FT::eType::BAM:
			if (fragType)
				FragDist(iName, prStat).Print(GetType(Distrib::eCType::LNORM), prDist);
			else
				ReadDist(iName, prStat).Print(GetType(Distrib::eCType::NORM), prDist);
			break;
		case FT::eType::FQ:
			if (fragType && Options::Assigned(oINPUT))
				Err(sWFormat + " for fragment distribution; " + sRExt).Throw();
			FqReadDist(iName).Print(GetType(Distrib::eCType::NORM), prDist);
			break;
		case FT::eType::DIST:
			dout << LF;
			Distrib(iName).Print(dout, Distrib::eCType(Options::GetIVal(oDTYPE)), prDist);
			break;
		default:
			Err(sWFormat + SepSCl + sRExt + " or FQ").Throw();
			break;
		}
	}
	catch (const Err& e) { ret = 1; cerr << e.what() << endl; }
	catch (const exception& e) { ret = 1; cerr << e.what() << endl; }
	//timer.Stop("wall-clock: ", false, true);
	timer.Stop();
	return ret;
}

// *********************** FragDist *********************************

// treats current read
bool FragDist::operator()()
{
	if (_uncheckedPE) {
		if (!File().IsPairedRead())
			Err("only paired-end reads are acceptable to call fragments distribution",
				File().CondFileName()).Throw();
		_uncheckedPE = false;
	}
	const Read read(File());
	const auto itMate = _waits.find(read.Numb);	// look for the read with given Numb

	if (itMate == _waits.end())					// is read not on the waiting list?
		_waits.emplace(read.Numb, read);		// add read the waiting list
	else {										// a mate
		const Read& mate = itMate->second;
		if (mate.Start != _pos[mate.Strand] || read.Start != _pos[read.Strand])	// not a duplicate
			AddFrag(mate, read);				// add uniq fragment into statistics
		else {
			if (_dupl)	AddFrag(mate, read);	// add dupl fragment into statistics
			_issues[0].Cnt++;
		}
		_pos[mate.Strand] = mate.Start;
		_pos[read.Strand] = read.Start;
		_waits.erase(itMate);					// remove read from the waiting list
		_cnt++;
	}
	return true;
}

// *********************** end of FragDist *********************************
