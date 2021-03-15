/************************************************************************************
bioCC is a fast advanced correlation calculator for basics bioinformatics file formats.
It computes Signal and Pearson correlation coefficients for densities, coverages and features.
Program allows to know correlation coefficients for the whole genome, for each chromosome
separately and for predefined regions inside chromosomes:
again for the all regions and for each region separately as well.

bioCC is designed to treat a bunch of files at once.

Copyright (C) 2017 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 3.12.2020
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
const string Product::Version = "1.0";
//const string Product::Descr = "advanced correlation calculator for bam/bed/wig files";
const string Product::Descr = "advanced correlation calculator";

const string OutFile = Product::Title + "_out.txt";
const string HelpOutFile = "duplicate standard output to " + OutFile + " file";
const string InFiles = "input files";

// --cc option: correlation coefficient notations
const char* CCs[] = { "P", "S" };			// corresponds to CCkey::eCC
// --total option: total coefficients notations
const char* prCCs[] = { "LOC", "TOT" };	// corresponds to eTotal; totalOFF is hidden
// --sort option: sorting type notations
const char* fsort[] = { "RGN", "CC" };		// corresponds to eRS; rsOFF is hidden
// --info option: types of info notations
const char* infos[] = { "LAC", "NM", "CNT", "STAT" };	// corresponds to eInfo; iNONE is hidden

const char* ForAligns = "For the alignments only";
const char* IgnoreBed = "Ignored for the ordinary beds";

// *** Options definition

enum eOptGroup { gINPUT, gTREAT, gTREAT_R, gOUTPUT, gOTHER };	// gOTHER should be the last 
const char* Options::OptGroups[] = {
	"Input", "Processing", "Region processing", "Output", "Other"
};
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

//	{ char,	str,	Signs,	type,	group,	defVal,	minVal,	maxVal,	strVal,	descr, addDescr }
// field 7: vUNDEF if value is prohibited
// field 6: vUNDEF if no default value should be printed
Options::Option Options::List[] = {
	{ 'a', "align",	fNone,	tENUM,	gINPUT,	FALSE,	vUNDEF, 2, NULL,
	"input bed files are alignments. Ignored for bam and wig", NULL },
	{ 'g', sGen,	fNone,	tNAME,	gINPUT, vUNDEF, 0, 0, NULL,
	"chromosome sizes file", NULL },
	{ HPH,"gap-len",fHidden,tINT,	gINPUT,	1000, 50, 1e5, NULL,
	"minimal length of undefined nucleotide region in genome\nwhich is declared as a gap.\nIgnored for the chromosome sizes file and for the ordinary beds", NULL },
	{ 'd', "dupl",	fHidden,tENUM,	gINPUT, TRUE,	0, 2, (char*)Options::Booleans,
	"accept duplicate reads.", ForAligns },
	{ 'l', "list",	fNone,	tNAME,	gINPUT, vUNDEF, 0, 0, NULL,
	"list of multiple input files.\nFirst (primary) file in list is comparing with others (secondary)", NULL },
	{ 'c', Chrom::Abbr,	fNone,	tNAME,	gTREAT, vUNDEF, 0, 0, NULL,	"treat specified chromosome only", NULL },
	{ 'r', "cc",		fNone,	tCOMB,	gTREAT,	CC::ccP, CC::ccP, CC::ccS, (char*)CCs,
	"correlation coefficient, in any order: ? - Pearson, ? - signal", NULL },
	{ 'f', "fbed",		fNone,	tNAME,	gTREAT_R, vUNDEF,	0, 0, NULL,
	"'template' ordinary bed file which features define compared regions.\n", IgnoreBed},
	{ 'e', "ext-len",	fNone,	tINT,	gTREAT_R,0, 0, 1e4, NULL,
	"length by which the features in primary file (for ordinary beds) or in\n'template' (for alignments and wigs) will be extended in both directions\nbefore treatment", NULL },
	{ 's', "ext-step",	fNone,	tINT,	gTREAT_R,0, 0, 500, NULL,
	"step of extending features in primary bed file;\nif 0 then no step calculation. For the ordinary beds only", NULL },
	{ 'R',	"pr-cc",	fNone,	tCOMB,	gOUTPUT, PrintMngr::pLOC, PrintMngr::pLOC, PrintMngr::pTOT, (char*)prCCs,
	"print coefficient, in any order:\n? - for each chromosome, ? - total", NULL },
	{ 'B', "bin-width",	fNone,	tFLOAT,	gOUTPUT,0, 0, 1.0F, NULL,
	"print frequency histogram with given bin width", NULL },
	//{ 'F', "fcc-sort ",	fOptnal,	tENUM,	gOUTPUT, rsOFF,	rsR, rsC, (char*)fsort,
	{ 'F', "fcc-sort",	fOptnal,	tENUM,	gOUTPUT, rsC,	rsR, rsC, (char*)fsort,
	"print region coefficients, sorted by: ? - regions, ? - coefficients", NULL },
	{ 'V', "verbose ",	fNone,	tENUM, gOUTPUT,	float(Obj::eInfo::NM), float(Obj::eInfo::LAC), float(Obj::eInfo::STAT), (char*)infos,
	"set verbose level:\n?  - laconic\n?   - file names\n?  - file names and number of items\n? - file names and statistics", NULL },
	{ 'w', "warn",	fHidden,tENUM, gOUTPUT,	FALSE,	vUNDEF, 2, NULL,
	"print each file's item ambiguity, if they exist.", NULL },
	{ 'o', sOutput,	fNone,	tENUM,	gOUTPUT,FALSE,	vUNDEF, 2, NULL, HelpOutFile.c_str(), NULL },
	{ 't', sTime,	fNone,	tENUM,	gOTHER,	FALSE,	vUNDEF, 2, NULL, sPrTime, NULL },
	{ HPH, sSumm,	fHidden,tSUMM,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, sPrSummary, NULL },
	{ 'v', sVers,	fNone,	tVERS,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, sPrVersion, NULL },
	{ 'h', sHelp,	fNone,	tHELP,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, sPrUsage, NULL }
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
	const FileList* fList = NULL;

	CC::Set(CC::eCC(Options::GetIVal(oCC)));
	PrintMngr::Init(Options::GetIVal(oPRCC), Obj::eInfo(Options::GetIVal(oVERB)));

	Chrom::SetCustomOption(oCHROM);
	if (Options::GetBVal(oOUTFILE))	dout.OpenFile(OutFile);
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		char** inFiles;
		short cntInFiles;
		const char* gName = Options::GetSVal(oGEN);		// genome

		// check file-list is set & exist
		const char* fListName = Options::GetSVal(oFILE_LIST);
		if (fListName) {
			fList = new FileList(FS::CheckedFileName(fListName));
			inFiles = fList->Files();
			if (!(cntInFiles = fList->Count()))
				Err(Err::MISSED, NULL, InFiles + " (no uncommented line in " +
					string(fListName) + ")").Throw();
		}
		else	inFiles = argv + fileInd, cntInFiles = argc - fileInd;
		if (cntInFiles < 2)			// check input files
			Err(Err::MISSED, NULL, cntInFiles ? "secondary " + InFiles : InFiles).Throw();

		if (!gName && !FS::HasExt(inFiles[0], FT::Ext(FT::eType::BAM)))
			Err(Options::OptionToStr(oGEN) + " is required for all the file types except BAM",
				inFiles[0]).Throw();

		ChromSizes cSizes(gName, NULL, true);
		DefRegions gRgn(cSizes, Options::GetIVal(oGAPLEN));
		CorrPair cPair(inFiles[0], gRgn, Options::GetSVal(oFBED), cntInFiles > 2);
		for (short i = 1; i < cntInFiles; i++)
			cPair.CalcCC(inFiles[i]);
	}
	catch (const Err & e) { ret = 1;	cerr << e.what() << LF; }
	catch (const exception & e) { ret = 1;	cerr << SPACE << e.what() << LF; }
	if (fList)	delete fList;
	timer.Stop(false);
	return ret;
}


/************************ class FileList ************************/
#ifdef OS_Windows
// Returns true if 'name' is file's pattern
bool	IsFilePattern(const char* name)
{
	return strchr(name, '*') != NULL || strchr(name, '?') != NULL;
}

string GetPath(const LPCTSTR name)
{
	const char* pch = strrchr(name, '/');
	return pch ? string(name, pch - name + 1) : strEmpty;
}

//// Fills vector by filenames according name template.
//// Works only if _UNICODE is not defined.
////	@files: vector of file names
////	@templ: name template, which can include '*' and '?' marks
//void FillFilesByTemplate(vector<string>& files, const LPCTSTR templ)
//{
//	string path = GetPath(templ);
//	WIN32_FIND_DATA ffd;
//
//	// check directory and estimate listFileNames capacity
//	HANDLE hFind = FindFirstFile(templ, &ffd);
//	if (hFind == INVALID_HANDLE_VALUE)
//		Err("bad file or content", templ).Throw();
//	if (files.capacity() == 0) {		// count files to reserve capacity
//		short count = 1;
//		for (; FindNextFile(hFind, &ffd); count++);
//		files.reserve(count);
//		hFind = FindFirstFile(templ, &ffd);
//	}
//	// fill the list
//	do	files.push_back(path + string(ffd.cFileName));	//  works only if _UNICODE isn't defined
//	while (FindNextFile(hFind, &ffd));
//	FindClose(hFind);
//}
#endif	// OS_Windows

//FileList::FileList(char* files[], short cntFiles) : _files(NULL), _memRelease(true)
//{
//#ifdef OS_Windows
//	short i;
//	bool hasTemplate = false;
//	for (i = 0; i < cntFiles; i++)
//		if (hasTemplate = IsFilePattern(files[i]))
//			break;
//	if (hasTemplate) {
//		// First we fill vector of file names, second initialize _files by this vector.
//		// It needs because files[] may contain file names and file template (pattern) as well,
//		// so in case of Windows we don't know the total amount of files beforehand.
//		vector<string> tmpFiles;
//		if (cntFiles > 1)
//			tmpFiles.reserve(cntFiles);
//		// else if it's a single template name, capacity will be reserved in FillFilesByTemplate(),
//		// or if it's a single common name, capacity will not be reserved at all
//
//		for (i = 0; i < cntFiles; i++)
//			if (IsFilePattern(files[i])) 
//				FillFilesByTemplate(tmpFiles, files[i]);
//			else		// real list of files
//				tmpFiles.push_back(files[i]);
//
//		_files = new char* [_count = short(tmpFiles.size())];
//		char* fname;
//		size_t len;
//		for (i = 0; i < _count; i++) {
//			_files[i] = fname = (char*)malloc(len = tmpFiles[i].length() + 1);	// free in ~FileList()
//			memcpy(fname, tmpFiles[i].c_str(), len);
//			//strcpy_s(fname, len, tmpFiles[i].c_str());
//		}
//	}
//	else
//#endif	// OS_Windows
//	{
//		_files = files;
//		_count = cntFiles;
//		_memRelease = false;
//	}
//}

FileList::FileList(const char* fName) : _count(0), _files(NULL)//, _memRelease(true)
{
	TabFile file(fName);
	vector<char*> tmpFiles;	// temporary vector because we don't know the proof capacity

	//== fill tmpFiles
	tmpFiles.reserve(file.EstCount());
	while (file.GetNextLine()) {
		const char* src = file.StrField(0);
		size_t len = 1 + strlen(src);
		char* dst = new char[len];
		memcpy(dst, src, len);
		tmpFiles.push_back(dst);
		_count++;
	}
	//== copy to _files
	_files = new char* [_count];
	memcpy(_files, tmpFiles.data(), _count*sizeof(char*));
	//== free tmpFiles
	//for (char* item : tmpFiles)	delete item, item = NULL;

}

/************************ end of class FileList ************************/
