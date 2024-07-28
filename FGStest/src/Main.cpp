/************************************************************************************
FGStest - Features Gold Standard test
-------------------------
Last modified: 07/27/2024
-------------------------
This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See GNU General Public License for more details.
************************************************************************************/

#include "Main.h"
#include "Options.h"
#include "Features.h"

using namespace std;

const string Product::Title = "FGStest";
const string Product::Version = "1.0";
const string Product::Descr = "Features Gold Standard test";

const char* ProgParam = "<in-file>";	// program parameter tip

// *** Options definition

enum eOptGroup { gOTHER };						// the only member - no gropus in help
const char* Options::OptGroups[] = { NULL };	// no gropus in help
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

const BYTE Options::Option::IndentInTabs = 3;

// { char, str, Signs (8: hidden), type, group, 
//	defVal (if NO_DEF then no default value printed),
//	minVal (if NO_VAL then value is prohibited), maxVal, strVal, descr, addDescr }
Options::Option Options::List[] = {
	//{ 'g', sGen,	tOpt::NONE,	tNAME,	gOTHER, NO_DEF, 0, 0, NULL, "chromosome sizes file" },
	{ 'c', sChrom,	tOpt::NONE,	tNAME,	gOTHER,	NO_DEF, 0, 0, NULL, sHelpChrom },
	{ 'S',"sample",	tOpt::OBLIG,tNAME,	gOTHER, NO_DEF, 0, 0, NULL, "sample file." },
	{ 'd',"min-dev",tOpt::NONE,	tINT,	gOTHER, 10, 0, 1000, NULL, "threshold deviation for writing a test feature to a dump file" },
	{ 's',"min-scr",tOpt::NONE,	tFLOAT,	gOTHER, 0, 0, 1, NULL, "threshold score for taking sample features into accounts" },
	{ 'w', "warn",	tOpt::HIDDEN,tENUM,	gOTHER, FALSE,	NO_VAL, 0, NULL, "print each feature ambiguity, if they exist" },
	{ 'D', "dump",	tOpt::NONE,	tNAME,	gOTHER,	NO_DEF,	0,	0, NULL, "output dump file name" },
	{ 'O', sOutput,	tOpt::FACULT,tNAME,	gOTHER,	NO_DEF,	0,	0, NULL, DoutHelp(ProgParam) },
	{ 't',	sTime,	tOpt::NONE,	tENUM,	gOTHER,	FALSE,	NO_VAL, 0, NULL, sHelpTime },
	{ 'v',	sVers,	tOpt::NONE,	tVERS,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sHelpVersion },
	{ HPH,	sSumm,	tOpt::HIDDEN,tSUMM,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sHelpSummary },
	{ 'h',	sHelp,	tOpt::NONE,	tHELP,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sHelpUsage },
};
const BYTE Options::OptCount = ArrCnt(Options::List);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ NO_DEF, ProgParam, true, ".bed[.gz] file containing peaks" }
};
const BYTE Options::UsageCount = ArrCnt(Options::Usages);

dostream dout;	// stream's duplicator

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
		const char* iName = FS::CheckedFileName(argv[fileInd]);	// input name
		//const char* gName = Options::GetSVal(oGEN);				// chrom sizes

		Options::SetDoutFile(oDOUT_FILE, iName);

		const Features smpl(FS::CheckedFileName(Options::GetSVal(oTEMPL)),
			nullptr, false, eOInfo::STD, true);
		if (!smpl.ChromCount())	return 0;
		const Features test(iName, nullptr, false, eOInfo::STD, true);
		if (!test.ChromCount())	return 0;

		FeaturesStatTuple fst(
			smpl,
			test,
			Options::GetFVal(oMIN_SCORE),
			short(Options::GetIVal(oMIN_DEV)),
			Options::GetSVal(oDUMP_FILE)
		);
		chrid	chrCount = 0;	// count of threated chromosomes

		FeaturesStatTuple::PrintHeader();
		for (auto it0 = smpl.cBegin(); it0 != smpl.cEnd(); it0++) {
			auto it1 = test.GetIter(CID(it0));
			if (it1 == test.cEnd())		continue;	// chrom not found
			fst.GetChromStat(it0, it1);
			chrCount++;
		}
		if (!chrCount)
			dout << Chrom::NoChromMsg() << " common to sample and test\n";
		else if (chrCount > 1)
			fst.PrintTotalStat();
	}
	catch (const Err& e) { ret = 1; cerr << e.what() << endl; }
	catch (const exception& e) { ret = 1; cerr << e.what() << endl; }
	timer.Stop();
	return ret;
}


const char* FeaturesStatTuple::BC::titles[2]{ "FP", "FN" };