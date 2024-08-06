/************************************************************************************
FGStest - Features Gold Standard test
-------------------------
Last modified: 08/06/2024
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
const string IssFileSuffix = ".issues";

// *** Options definition

enum eOptGroup { gINPUT, gOUTPUT, gOTHER };	// gOTHER should be the last 
const char* Options::OptGroups[] = { "Input", "Output", "Other" };
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

const BYTE Options::Option::IndentInTabs = 3;

// { char, str, Signs (8: hidden), type, group, 
//	defVal (if NO_DEF then no default value printed),
//	minVal (if NO_VAL then value is prohibited), maxVal, strVal, descr, addDescr }
Options::Option Options::List[] = {
	//{ 'g', sGen,	tOpt::NONE,	tNAME,	gINPUT, NO_DEF, 0, 0, NULL, "chromosome sizes file" },
	{ 'c', sChrom,	tOpt::NONE,	tNAME,	gINPUT,	NO_DEF, 0, 0, NULL, sHelpChrom },
	{ 'S',"sample",	tOpt::OBLIG,tNAME,	gINPUT, NO_DEF, 0, 0, NULL, "sample file." },
	{ 'C',"min-cdev",tOpt::NONE,tINT,	gINPUT, 10, 0, 1000, NULL, "threshold centre deviation for writing a test feature to an issues file" },
	{ 'W',"min-wdev",tOpt::NONE,tFLOAT,	gINPUT, 0, 1, 100, NULL, "threshold width deviation for writing a test feature to an issues file" },
	{ 's',"min-scr",tOpt::NONE,	tFLOAT,	gINPUT, 0, 0, 1, NULL, "threshold score for taking sample features into accounts" },
	{ 'e', "expand",tOpt::NONE,	tINT,	gINPUT, 0, 0, 500, NULL, "virtually expand sample features while searching for sample/test intersecting" },
	{ 'w', "warn",	tOpt::HIDDEN,tENUM,	gOUTPUT,FALSE,	NO_VAL, 0, NULL, "print each feature ambiguity, if they exist" },
	{ 'I', "issues",tOpt::FACULT,tNAME,	gOUTPUT,NO_DEF,	0,	0, NULL,
	OptFileNameHelp("output locused issues", ProgParam, IssFileSuffix, FT::Ext(FT::BED))},
	{ 'O', sOutput,	tOpt::FACULT,tNAME,	gOUTPUT,NO_DEF,	0,	0, NULL, DoutHelp(ProgParam) },
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


// Returns C-string containing the output issue file name
const char* GetIssFileName(string& oName, const char* defName, const string& suffix)
{
	if (Options::Assigned(oISSUE_FILE)) {
		oName = FS::ComposeFileName(Options::GetSVal(oISSUE_FILE), defName, suffix, FT::Ext(FT::BED));
		if (oName.length()) {
			oName += FT::Ext(FT::BED);
			return oName.c_str();
		}
	}
	return NULL;
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
		const char* iName = FS::CheckedFileName(argv[fileInd]);	// input name
		auto expVal = USHORT(Options::GetFVal(oEXPAND));

		Options::SetDoutFile(oDOUT_FILE, iName);

		dout << "sample: ";
		Features smpl(FS::CheckedFileName(Options::GetSVal(oTEMPL)),
			nullptr, false, eOInfo::STD);
		if (!smpl.ChromCount())	return 0;
		// check expansion value
		if (expVal) {
			auto minDistance = smpl.GetMinDistance();
			if (expVal >= minDistance / 2) {
				cerr << Options::OptionToStr(oEXPAND, true)
					<< " leads to the intersection of the sample features (the minimum distance between them is "
					<< minDistance << ")\n";
				return 1;
			}
		}
		dout << "test:\t";
		const Features test(iName, nullptr, false, eOInfo::STD);
		if (!test.ChromCount())	return 0;

		string oName;			// issues output file name
		FeaturesStatTuple fst(
			smpl,
			test,
			expVal,
			Options::GetFVal(oMIN_SCORE),
			short(Options::GetFVal(oMIN_CDEV)),
			Options::GetFVal(oMIN_WDEV),
			GetIssFileName(oName, iName, IssFileSuffix)
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