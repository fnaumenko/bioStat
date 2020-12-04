/************************************************************************************
fragDist - fragment's distribution - is designed to calculate the size distribution parameters of
paired-end fragments OR reads

The frequency distribution is displayedin the form of lines,
each of which contains one pair <fragment length><TAB><frequency>.
Then this output is transferred to the Excel, which allows to plot it.

Copyright (C) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 1.12.2020
-------------------------

This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See GNU General Public License for more details.
************************************************************************************/

#include "fragDist.h"
#include <memory>		// unique_ptr

using namespace std;

const string Product::Title = "fragDist";
const string Product::Version = "1.0";
const string Product::Descr = "PE-fragment-size/read-length distribution parameters caller";

const char* ProgParam = "<in-file>";	// program parameter tip
const string OutFileExt = FT::Ext(FT::eType::DIST);
const string HelpOutFile = sFileDuplBegin + string(ProgParam) + OutFileExt + sFileDuplEnd;

enum eOptGroup { oOTHER };						// the only member - no gropus in help
const char* Options::OptGroups[] = { NULL };	// no gropus in help
const BYTE Options::GroupCount = sizeof(Options::OptGroups) / sizeof(char*);

// --input option
const char* sinputs[] = { "FRAG","READ" };					// corresponds to Inp


// { char, str, Signs (8: hidden), type, group, 
//	defVal (if NO_DEF then no default value printed),
//	minVal (if NO_VAL then value is prohibited), maxVal, strVal, descr, addDescr }
Options::Option Options::List[] = {
	//{ 'c',Chrom::Abbr,fNone,tNAME,	oOTHER,	NO_DEF, 0, 0, NULL, "treat specified chromosome only", NULL },
	{ 'd', "no-dup",fHidden,tENUM,	oOTHER,	FALSE,	NO_VAL, 0, NULL, "reject duplicate reads", NULL },
	{ 'I', "info",	fHidden,tENUM,	oOTHER,	FALSE,	NO_VAL,	0, NULL,
	"print statistical info about alignment", NULL },
	{ 'w', "warn",	fHidden,tENUM,	oOTHER, FALSE,	NO_VAL, 0, NULL,
	"print each read ambiguity, if they exist" },
	{ 'i', "inp",	fNone,	tENUM,	oOTHER, float(Inp::FRAG), float(Inp::FRAG), float(Inp::UNDEF), (char*)sinputs,
	"input data to call distribution: ? - fragments, ? - reads", NULL },
	{ 'n', "norm",	fNone,	tENUM,	oOTHER,	FALSE,	NO_VAL, 0, NULL, "call normal distribution parameters anyway", NULL },
	{ 'D', "dist",	fNone,	tENUM,	oOTHER,	FALSE,	NO_VAL, 0, NULL, "print frequency distribution", NULL },
	{ 'o', sOutput,	fOptnal,tNAME,	oOTHER,	NO_DEF,	0,	0, NULL, HelpOutFile.c_str(), NULL },
	{ 't',	sTime,	fNone,	tENUM,	oOTHER,	FALSE,	NO_VAL, 0, NULL, sPrTime, NULL },
	{ HPH,	sSumm,	fHidden,tSUMM,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrSummary, NULL },
	{ 'v',	sVers,	fNone,	tVERS,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrVersion, NULL },
	{ 'h',	sHelp,	fNone,	tHELP,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrUsage, NULL }
};
const BYTE	Options::OptCount = sizeof(Options::List) / sizeof(Options::Option);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, ProgParam, true, "paired-end alignment in bam/bed format OR reads in fq/bam/bed format" }
};
const BYTE Options::UsageCount = sizeof(Options::Usages) / sizeof(Options::Usage);

dostream dout;	// stream's duplicator

/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv, ProgParam);
	if (fileInd < 0)	return 1;		// wrong option or tip output
	int ret = 0;						// main() return code

	//Chrom::SetCustomOption(oCHROM, true);
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* iName = FS::CheckedFileName(argv[fileInd]);		// input
		const char* oName = Options::GetSVal(oOUTFILE);				// output
		FT::eType type = FT::GetType(iName);
		Inp input = Inp(Options::GetIVal(oINPUT));
		bool prDist = Options::GetBVal(oPR_DIST);

		// while input file is distribution, add "_out" to output file name
		// only if isn't specified or the same as input
		if (Options::Assigned(oOUTFILE)
		&& !dout.OpenFile(Options::GetFileName(oOUTFILE, iName,
			(type == FT::eType::DIST && 
				(!oName || !strcmp(FS::FileNameWithoutExt(iName).c_str(), oName)) ?
				"_out" : strEmpty)
			+ OutFileExt))
			)
			Err(Err::FailOpenOFile).Throw();

		// take distribution
		if (type == FT::eType::BED || type == FT::eType::BAM) {
			if (input == Inp::FRAG) {
				ChromSizes cSizes;
				Reads reads(ProgParam, iName, cSizes,
					Obj::eInfo(Options::GetBVal(oINFO) ? Obj::eInfo::STAT : Obj::eInfo::STD),
					true, true, Options::GetBVal(oALARM),
					true, !Options::GetBVal(oNO_DUPL), -1);

				FragDist(reads).PrintDist(prDist);
			}
			else {			// reads
				unique_ptr<DataInFile> file;
				if (type == FT::eType::BED)
					file.reset(new BedInFile(string(iName), FT::eType::ABED, 0, true, false));
				else
					file.reset(new BamInFile(string(iName), false));
				ReadDist(*file).PrintDist(true, prDist);
			}
		}
		else if (type == FT::eType::FQ)
			if (input == Inp::READ) {
				FqFile file(iName);
				ReadDist(file).PrintDist(true, prDist);
			}
			else
				Err(string(iName) + ": wrong format for fragment distribution; should be BED or BAM").Throw();
		else if (type == FT::eType::DIST) {
			//dout << LF << iName << COLON << LF;
			LenFreq(iName).Print(dout, LenFreq::eType::UNDEF, Options::GetBVal(oNORM), prDist);
		}
		else
			Err(string(iName) + ": wrong format; should be BED or BAM or FQ").Throw();
	}
	catch (const Err& e)		{ ret = 1; cerr << e.what() << endl; }
	catch (const exception& e)	{ ret = 1; cerr << e.what() << endl; }
	timer.Stop("wall-clock: ", false, true);
	return ret;
}


// *********************** FragDist *********************************

// Add read to statistics if its mate is waiting already, otherwhise add it to the waiting list
//	@rit: read's iterator
void FragDist::Add(const Reads::cItemsIter& rit)
{
	auto mit = _waits.find(rit->Numb);	// look for the read with given Numb

	rit->Pos;
	if (mit == _waits.end())							// is read not on the waiting list?
		_waits[rit->Numb] = rit;						// add read the waiting list
	else if (mit->second->Strand != rit->Strand) {		// is it mate, not duplicate read?
		IncrStat(mit->second, rit);						// insert fragment into statistics
		_waits.erase(rit->Numb);						// remove read from the waiting list
	}
}

// Constructor by alignment; the instance is initialized according to the reads
//	@test: paired-end reads collection
FragDist::FragDist(Reads& test)
{
	auto cit = test.cBegin();	// chrom iterator
	Read::FixedLen = test.ReadLen();

	_waits.rehash(test.ItemsCount(CID(cit)) / 2);	// mem reserve to avoid excessive rehash
	for (; cit != test.cEnd(); cit++) {
		const auto end = test.ItemsEnd(cit);		// curr chrom last read iterator
		for (auto rit = test.ItemsBegin(cit); rit != end; rit++)	// read iteration through chrom 
			Add(rit);
		Clear();
	}
}
// *********************** end of FragDist *********************************
