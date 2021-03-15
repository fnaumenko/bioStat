/************************************************************************************
callDist is designed to calculate the size distribution parameters of
paired-end fragments OR reads

The frequency distribution is displayedin the form of lines,
each of which contains one pair <fragment length><TAB><frequency>.
Then this output is transferred to the Excel, which allows to plot it.

Copyright (C) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 15.03.2021
-------------------------

This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See GNU General Public License for more details.
************************************************************************************/

#include "callDist.h"
#include <memory>		// unique_ptr

using namespace std;

const string Product::Title = "callDist";
const string Product::Version = "1.0";
const string Product::Descr = "PE-fragment-size/read-length distribution parameters caller";

const char* ProgParam = "<in-file>";	// program parameter tip
const string OutFileExt = FT::Ext(FT::eType::DIST);
const string HelpOutFile = sFileDuplBegin + string(ProgParam) + OutFileExt + sFileDuplEnd;

const char* inputs[] = { "FRAG","READ" };	// input option; corresponds to Inp
const char* dTypes[] = { "N","LN","G" };	// input distrib option; corresponds to InType	

// *** Options definition

enum eOptGroup { gOTHER };						// the only member - no gropus in help
const char* Options::OptGroups[] = { NULL };	// no gropus in help
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

// { char, str, Signs (8: hidden), type, group, 
//	defVal (if NO_DEF then no default value printed),
//	minVal (if NO_VAL then value is prohibited), maxVal, strVal, descr, addDescr }
Options::Option Options::List[] = {
	//{ 'c',Chrom::Abbr,fNone,tNAME,	gOTHER,	NO_DEF, 0, 0, NULL, "treat specified chromosome only", NULL },
	{ 'd', "no-dup",fHidden,tENUM,	gOTHER,	FALSE,	NO_VAL, 0, NULL, "reject duplicate reads", NULL },
	{ 'I', "info",	fHidden,tENUM,	gOTHER,	FALSE,	NO_VAL,	0, NULL,
	"print statistical info about alignment", NULL },
	{ 'w', "warn",	fHidden,tENUM,	gOTHER, FALSE,	NO_VAL, 0, NULL,
	"print each read ambiguity, if they exist" },
	{ 'i', "inp",	fNone,	tENUM,	gOTHER, float(InpType::FRAG), float(InpType::FRAG), ArrCnt(inputs),
	 (char*)inputs, "input data to call distribution: ? - fragments, ? - reads", NULL },
	{ 'D',"dist",	fNone,	tCOMB,	gOTHER, float(LenFreq::eCType::LNORM), float(LenFreq::eCType::NORM), 
	ArrCnt(dTypes),	(char*)dTypes,
	"called distribution, in any order: ? - normal, ? - lognormal, ? - Gamma", NULL },
	{ 'p', "pr-dist",fNone,	tENUM,	gOTHER,	FALSE,	NO_VAL, 0, NULL, "print original frequency distribution", NULL },
	{ 'o', sOutput,	fOptnal,tNAME,	gOTHER,	NO_DEF,	0,	0, NULL, HelpOutFile.c_str(), NULL },
	{ 't',	sTime,	fNone,	tENUM,	gOTHER,	FALSE,	NO_VAL, 0, NULL, sPrTime, NULL },
	{ HPH,	sSumm,	fHidden,tSUMM,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrSummary, NULL },
	{ 'v',	sVers,	fNone,	tVERS,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrVersion, NULL },
	{ 'h',	sHelp,	fNone,	tHELP,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrUsage, NULL }
};
const BYTE Options::OptCount = ArrCnt(Options::List);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, ProgParam, true, "paired-end alignment in bam/bed format OR reads in fq/bam/bed format" }
};
const BYTE Options::UsageCount = ArrCnt(Options::Usages);

dostream dout;	// stream's duplicator

// Returns explicitly defined by user combo type, otherwise default combo type
LenFreq::eCType GetType(LenFreq::eCType defType)
{
	return Options::Assigned(oDTYPE) ? LenFreq::eCType(Options::GetIVal(oDTYPE)) : defType;
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
		FT::eType ftype = FT::GetType(iName);
		InpType inpType = InpType(Options::GetIVal(oINPUT));
		bool prDist = Options::GetBVal(oPR_DIST);

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
		dout << iName << SepCl;	cout.flush();
		switch (ftype) {
		case FT::eType::BED:
		case FT::eType::BAM:
			if (inpType == InpType::FRAG) {
				ChromSizes cSizes;
				Reads reads(NULL, iName, cSizes,
					Obj::eInfo(Options::GetBVal(oINFO) ? Obj::eInfo::STAT : Obj::eInfo::NONE),
					false, true, Options::GetBVal(oALARM),
					true, !Options::GetBVal(oNO_DUPL), -1);
				FragDist(reads).Print(GetType(LenFreq::eCType::LNORM), prDist);
			}
			else {			// reads
				unique_ptr<DataInFile> file;
				if (ftype == FT::eType::BED)
					file.reset(new BedInFile(iName, FT::eType::ABED, 0, true, false));
				else
					file.reset(new BamInFile(iName, false));
				ReadDist(*file).Print(GetType(LenFreq::eCType::NORM), prDist);
			}
			break;
		case FT::eType::FQ:
			if (inpType == InpType::FRAG && Options::Assigned(oINPUT))
				Err("wrong format for fragment distribution; should be BED or BAM").Throw();
			else {
				FqFile file(iName);
				ReadDist(file).Print(GetType(LenFreq::eCType::NORM), prDist);
			}
			break;
		case FT::eType::DIST:
			LenFreq(iName).Print(dout, LenFreq::eCType(Options::GetIVal(oDTYPE)), prDist);
			break;
		default:
			Err("wrong format; should be BED or BAM or FQ").Throw();
			break;
		}
	}
	catch (const Err& e) { ret = 1; cerr << e.what() << endl; }
	catch (const exception& e) { ret = 1; cerr << e.what() << endl; }
	timer.Stop("wall-clock: ", false, true);
	return ret;
}


const char* LenDist::ItemTitles[] = { "fragments", "reads" };

// *********************** FragDist *********************************

// Add read to statistics if its mate is waiting already, otherwhise put it on the waiting list
//	@rit: read's iterator
void FragDist::PutRead(const Reads::cItemsIter& rit)
{
	auto mit = _waits.find(rit->Numb);		// look for the read with given Numb

	if (mit == _waits.end())							// is read not on the waiting list?
		_waits[rit->Numb] = rit;						// add read the waiting list
	else if (mit->second->Strand != rit->Strand) {		// is it mate, not duplicate read?
		AddFrag(mit->second, rit);						// insert fragment into statistics
		_waits.erase(rit->Numb);						// remove read from the waiting list
	}
}

// Constructor by alignment; the instance is initialized according to the reads
//	@test: paired-end reads collection
FragDist::FragDist(Reads& test) : LenDist(InpType::FRAG)
{
	auto cit = test.cBegin();	// chrom iterator
	Read::FixedLen = test.ReadLen();

	_waits.rehash(test.ItemsCount(CID(cit)) / 2);	// mem reserve to avoid excessive rehash
	for (; cit != test.cEnd(); cit++) {
		const auto end = test.ItemsEnd(cit);		// curr chrom last read iterator
		for (auto rit = test.ItemsBegin(cit); rit != end; rit++)	// read iteration through chrom 
			PutRead(rit);
		Clear();
	}
}
// *********************** end of FragDist *********************************
