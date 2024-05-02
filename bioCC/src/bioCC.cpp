/************************************************************************************
bioCC is a fast advanced correlation calculator for basics bioinformatics file formats.
It computes Signal and Pearson correlation coefficients for densities, coverages and features.
Program allows to know correlation coefficients for the whole genome, for each chromosome
separately and for predefined regions inside chromosomes:
again for the all regions and for each region separately as well.

bioCC is designed to treat a bunch of files at once.

Copyright (C) 2017 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 05/02/2024
-------------------------

This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the	GNU General Public License for more details.
************************************************************************************/

#include "Calc.h"
#include <string>
#include <algorithm>

const string Product::Title = "bioCC";
const string Product::Version = "2.0";
const string Product::Descr = "advanced correlation calculator";

const string OutFile = Product::Title + "_out.txt";
const string HelpOutFile = sFileDuplBegin + OutFile + sFileDuplEnd;
const string InFiles = "input files";

// --ext option
//const Options::PairVals ext(0, 0, 0, 0, 1000, 10000);	// defStep, defLen, minStep, minLen, maxStep, maxLen
// --total option: total coefficients notations
const char* prCCs[] = { "LOC", "TOT" };		// corresponds to eTotal; totalOFF is hidden
// --sort option: sorting type notations
const char* fsort[] = { "RGN", "CC" };		// corresponds to eRS; rsOFF is hidden
// --oinfo option: types of oinfo notations
const char* infos[] = { "LAC", "NM", "ITEM", "STAT" };	// corresponds to eOInfo; iNONE is hidden

const char* ForAligns = "For the alignments only";
const char* IgnoreBed = "Ignored for the ordinary beds";

// *** Options definition

enum eOptGroup { gINPUT, gTREAT, gOUTPUT, gOTHER };	// gOTHER should be the last 
const char* Options::OptGroups[] = {
	"Input", "Region processing", "Output", "Other"
};
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

const BYTE Options::Option::IndentInTabs = 3;

//	{ char,	str,	Signs,	type,	group,	defVal,	minVal,	maxVal,	strVal,	descr, addDescr }
// field 7: vUNDEF if value is prohibited
// field 6: vUNDEF if no default value should be printed
Options::Option Options::List[] = {
	{ 'a', "align",	tOpt::NONE,	tENUM,	gINPUT,	FALSE,	vUNDEF, 2, NULL,
	"input bed files are alignments. Ignored for bam and wig", NULL },
	{ 'g', sGen,	tOpt::NONE,	tNAME,	gINPUT, vUNDEF, 0, 0, NULL,
	"chromosome sizes file", NULL },
	{ 'c', Chrom::Abbr,	tOpt::NONE,	tNAME,	gINPUT, vUNDEF, 0, 0, NULL,	"treat specified chromosome only", NULL },
	{ HPH,"gap-len",tOpt::HIDDEN,tINT,	gINPUT,	1000, 50, 1e5, NULL,
	"minimal length of undefined nucleotide region in genome\nwhich is declared as a gap.\nIgnored for the chromosome sizes file and for the ordinary beds", NULL },
	{ 'd', "dup",	tOpt::NONE,	tENUM,	gINPUT, TRUE,	0, 2, (char*)Booleans, "allow duplicate reads.", ForAligns },
	{ 'o', "overl",	tOpt::NONE,	tENUM,	gINPUT, FALSE,	0, 2, (char*)Booleans,
	"allow (and merge) overlapping features. For the ordinary beds only", NULL },
	{ 'l', "list",	tOpt::NONE,	tNAME,	gINPUT, vUNDEF, 0, 0, NULL,
	"list of multiple input files.\nFirst (primary) file in list is comparing with others (secondary)", NULL },
	{ 'f', "fbed",	tOpt::NONE,	tNAME,	gTREAT, vUNDEF,	0, 0, NULL,
	"'template' ordinary bed file which features define compared regions.\n", IgnoreBed},
	{ 'e', "ext-len",	tOpt::NONE,	tINT,	gTREAT,0, 0, 2e4, NULL,
	"length by which the features in primary file (for ordinary beds) or in\n'template' (for alignments and wigs) will be extended in both directions\nbefore treatment", NULL },
	{ 's', "ext-step",	tOpt::NONE,	tINT,	gTREAT,0, 0, 500, NULL,
	"step of extending features in primary bed file;\nif 0 then no step calculation. For the ordinary beds only", NULL },
	{ 'R',	"pr-cc",	tOpt::NONE,	tCOMB,	gOUTPUT, PrintMngr::LOC, PrintMngr::LOC, PrintMngr::TOT, (char*)prCCs,
	"print coefficient, in any order:\n? - for each chromosome, ? - total", NULL },
	{ 'B', "bin-width",	tOpt::NONE,	tFLOAT,	gOUTPUT,0, 0, 1.0F, NULL,
	"print frequency histogram with given bin width", NULL },
	//{ 'F', "fcc-sort ",	fOptnal,	tENUM,	gOUTPUT, rsOFF,	rsR, rsC, (char*)fsort,
	{ 'F', "fcc-sort",	tOpt::NONE,	tENUM,	gOUTPUT, eRS::rsOFF,	eRS::rsR, eRS::rsC, (char*)fsort,
	"print region coefficients, sorted by: ? - regions, ? - coefficients", NULL },
	{ 'V', "verbose ",	tOpt::NONE,	tENUM, gOUTPUT,	float(eOInfo::NM), float(eOInfo::LAC), float(eOInfo::STAT), (char*)infos,
	"set verbose level:\n?  - laconic\n?   - file names\n? - file names and number of items\n? - file names and items statistics", NULL },
	{ 'w', "write",	tOpt::HIDDEN,tENUM,	gOUTPUT,FALSE,	vUNDEF, 2, NULL,
	"write each inner representation to file with '_out' suffix", NULL },
	{ 'O', sOutput,	tOpt::FACULT,tNAME,	gOUTPUT,NO_DEF,	0,	0, NULL, HelpOutFile.c_str(), NULL },
	{ 't', sTime,	tOpt::NONE,	tENUM,	gOTHER,	FALSE,	vUNDEF, 2, NULL, sPrTime, NULL },
	{ HPH, sSumm,	tOpt::HIDDEN,tSUMM,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, sPrSummary, NULL },
	{ 'v', sVers,	tOpt::NONE,	tVERS,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, sPrVersion, NULL },
	{ 'h', sHelp,	tOpt::NONE,	tHELP,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, sPrUsage, NULL }
};
const BYTE	Options::OptCount = ArrCnt(Options::List);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, "<file1> <file2> ...", true, NULL },
	{ oFILE_LIST, NULL, true, NULL }
};
const BYTE Options::UsageCount = ArrCnt(Options::Usages);

dostream dout;	// stream's duplicator

/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv);
	if (fileInd < 0)	return 1;		// wrong option or tip output
	int ret = 0;						// main() return code

	Chrom::SetUserChrom(Options::GetSVal(oCHROM));
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		char** inFiles;
		short inFilesCnt;
		unique_ptr<FileList> fList;

		// check file-list is set & exist
		const char* fListName = Options::GetSVal(oFILE_LIST);
		if (fListName) {
			fList.reset(new FileList(FS::CheckedFileName(fListName)));
			inFiles = fList->Files();
			inFilesCnt = fList->Count();
			if (!inFilesCnt)
				Err(Err::MISSED, NULL, InFiles + " (no significant line in " + fListName + ')').
				Throw();
		}
		else
			inFiles = argv + fileInd,
			inFilesCnt = argc - fileInd;
		if (inFilesCnt < 2)			// check input files
			Err(Err::MISSED, NULL, inFilesCnt ? "secondary " + InFiles : InFiles).Throw();

		// set genom
		const char* gName = Options::GetSVal(oGENOM);
		if (!gName && !FS::HasExt(inFiles[0], FT::Ext(FT::eType::BAM)))
			Err(Options::OptionToStr(oGENOM) + " is required for all the file types except BAM",
				inFiles[0]).Throw();

		// set output file
		if (Options::Assigned(oOUTFILE))
			if (!dout.OpenFile(FS::ComposeFileName(Options::GetSVal(oOUTFILE), OutFile.c_str())))
				Err(Err::FailOpenOFile).Throw();

		ChromSizes cSizes(gName, true);
		DefRegions gRgn(cSizes, Options::GetIVal(oGAPLEN));
		CorrPair cPair(inFiles[0], gRgn, Options::GetSVal(oFBED), inFilesCnt > 2);
		for (short i = 1; i < inFilesCnt; i++)
			cPair.CalcCC(inFiles[i]);
	}
	catch (const Err & e)		{ ret = 1; cerr << e.what() << LF; }
	catch (const exception & e) { ret = 1; cerr << SPACE << e.what() << LF; }
	timer.Stop();
	return ret;
}

/************************ class FileList ************************/
//#ifdef OS_Windows
//// Returns true if 'name' is file's pattern
//bool	IsFilePattern(const char* name)
//{
//	return strchr(name, '*') != NULL || strchr(name, '?') != NULL;
//}
//
//string GetPath(const LPCTSTR name)
//{
//	const char* pch = strrchr(name, '/');
//	return pch ? string(name, pch - name + 1) : strEmpty;
//}
//#endif	// OS_Windows

FileList::FileList(const char* fName)
{
	TabReader file(fName);
	vector<char*> tmpFiles;	// we don't know the proof capacity

	//== fill tmpFiles
	tmpFiles.reserve(file.EstLineCount());
	while (file.GetNextLine()) {
		const char* src = file.StrField(0);
		const size_t len = 1 + strlen(src);
		char* dst = new char[len];		// will be free in destructor
		memcpy(dst, src, len);
		tmpFiles.push_back(dst);
		_count++;
	}
	//== copy to _files
	_files = new char* [_count];		// will be free in destructor
	move(tmpFiles.data(), tmpFiles.data() + _count, _files);
}

FileList::~FileList()
{
	if (_files) {
		for (int i = 0; i < _count; delete[] _files[i++]);
		delete[] _files;
	}
}

#ifdef MY_DEBUG
void FileList::Print() const
{
	if (_files)
		for (short i = 0; i < _count; i++)
			cout << _files[i] << LF;
	else
		cout << "Empty\n";
}
#endif

/************************ end of class FileList ************************/
