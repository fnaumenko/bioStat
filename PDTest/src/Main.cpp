/************************************************************************************
PDTest - Peak Detectors test
-------------------------
Last modified: 07/24/2024
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

const string Product::Title = "PDTest";
const string Product::Version = "1.0";
const string Product::Descr = "Peak detectors Statistical Test";

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
	//{ 'g', sGen,	tOpt::NONE,	tNAME,	gOTHER, NO_DEF, 0, 0, NULL, "chromosome sizes file", NULL },
	//{ 'c',Chrom::Abbr,tOpt::NONE,tNAME,	gOTHER,	NO_DEF, 0, 0, NULL, "treat specified chromosome only", NULL },
	{ 'T',"templ",	tOpt::OBLIG,tNAME,	gOTHER, NO_DEF, 0, 0, NULL, "template BS file.", NULL },
	{ 's',"min-scr",tOpt::NONE,	tFLOAT,	gOTHER, 0, 0, 1, NULL, "threshold score for taking template features into accounts", NULL },
	{ 'w', "warn",	tOpt::HIDDEN,tENUM,	gOTHER, FALSE,	NO_VAL, 0, NULL, "print each read ambiguity, if they exist" },
	{ 'O', sOutput,	tOpt::NONE,	tNAME,	gOTHER,	NO_DEF,	0,	0, NULL, "output file name", NULL },
	{ 't',	sTime,	tOpt::NONE,	tENUM,	gOTHER,	FALSE,	NO_VAL, 0, NULL, sPrTime, NULL },
	{ 'v',	sVers,	tOpt::NONE,	tVERS,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrVersion, NULL },
	{ HPH,	sSumm,	tOpt::HIDDEN,tSUMM,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrSummary, NULL },
	{ 'h',	sHelp,	tOpt::NONE,	tHELP,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrUsage, NULL },
};
const BYTE Options::OptCount = ArrCnt(Options::List);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ NO_DEF, ProgParam, true, ".bed[.gz] file containing peaks" }
};
const BYTE Options::UsageCount = ArrCnt(Options::Usages);


/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv, ProgParam);
	if (fileInd < 0)	return 1;		// wrong option or tip output
	int ret = 0;						// main() return code

	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* iName = FS::CheckedFileName(argv[fileInd]);	// input name
		const char* oName = Options::GetSVal(oOUTFILE);			// output name
		//const char* gName = Options::GetSVal(oGEN);				// chrom sizes

		//auto name = FS::ComposeFileName(oName, iName);
		//cout << name << LF;
		//return 0;

		const Features tmpl(FS::CheckedFileName(Options::GetSVal(oTEMPL)),
			nullptr, false, eOInfo::STD, true);
		const Features test(iName, nullptr, false, eOInfo::STD, true);

		FeaturesStatTuple fst(tmpl, test, Options::GetFVal(oMIN_SCORE), oName);
		chrid	chrCount = 0;	// count of threated chromosomes

		FeaturesStatTuple::PrintHeader();
		for (auto it0 = tmpl.cBegin(); it0 != tmpl.cEnd(); it0++) {
			auto it1 = test.GetIter(CID(it0));
			if (it1 == test.cEnd())		continue;	// chrom not found

			fst.SetChrom(it0, it1);
			fst.Treat();
			fst.PrintStat(CHR);
			fst.ResetChrom();
			chrCount++;
		}
		//if (chrCount > 1) 
			fst.PrintStat(ALL);
	}
	catch (const Err& e) { ret = 1; cerr << e.what() << endl; }
	catch (const exception& e) { ret = 1; cerr << e.what() << endl; }
	timer.Stop();
	return ret;
}


const char* FeaturesStatTuple::BC::titles[2]{ "FP", "FN" };
