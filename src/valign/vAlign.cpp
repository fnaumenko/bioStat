/************************************************************************************
vAlign is a fast verifier of reads of aligned DNA sequence,
recalled from initial artificial FastQ sequence.
It compares the original and actual coordinates of each read
and prints statistics of right and wrong mappings.

Copyright (C) 2017 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 1.12.2020
-------------------------

This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the	GNU General Public License for more details.
************************************************************************************/

#include "Data.h"
#include "vAlign.h"

using namespace std;

const string Product::Title = "vAlign";
const string Product::Version = "1.0";
const string Product::Descr = "alignment verifier";

const char* ProgParam = "<in-file>";	// program parameter tip
const string OutFileSuff = "_valign.txt";
const string HelpOutFile = sFileDuplBegin + string(ProgParam) + OutFileSuff + sFileDuplEnd;

// --info option: types of info notations
const char* infos[] = { "NM", "CNT", "STAT" };	// corresponds to eInfo; iNONE and iLAC are hidden

// *** Options definition

enum eOptGroup { oTREAT, oOTHER };	// oOTHER should be the last 
const char* Options::OptGroups[] = { "Treatment", "Other" };
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

//	{ char,	str,	Signs,	type,	group,	defVal,	minVal,	maxVal,	strVal,	descr, addDescr }
// field 7: vUNDEF if value is prohibited
// field 6: vUNDEF if no default value should be printed
Options::Option Options::List[] = {
	{ 'g', sGen,	fOblig, tNAME, oTREAT, vUNDEF, 0, 0, NULL,
	"reference genome library or single nucleotide sequence.", NULL },
	{ 'c',Chrom::Abbr,fNone,tNAME,	oTREAT, vUNDEF, 0, 0, NULL,
	"treat specified chromosome only. For reference genome only", NULL },
	{ HPH,"min-scr",fNone,	tINT,	oTREAT, vUNDEF, 0, 1000, NULL, "score threshold for treated reads", NULL },
	{ HPH,"char-case",fNone,tENUM,	oTREAT, FALSE,	0, 2, (char*)Options::Booleans,
	"recognize uppercase and lowercase characters in template and test\nas different", NULL },
	{ 'i', "info",	fHidden,tENUM, oTREAT,	float(Obj::eInfo::STAT), float(Obj::eInfo::NM), float(Obj::eInfo::STAT), (char*)infos,
	"print information about file:\n? - name only, ? - number of reads, ? - statistics", NULL },
	{ 'w', "warn",	fHidden,tENUM,	oTREAT, FALSE,	vUNDEF, 2, NULL,
	"print each read ambiguity, if they exist", NULL },
	{ 'o', sOutput,	fOptnal,tNAME,	oTREAT,NO_DEF,	0,	0, NULL, HelpOutFile.c_str(), NULL },
	{ 't', sTime,	fNone,	tENUM,	oOTHER,	FALSE,	NO_VAL, 2, NULL, sPrTime, NULL },
	{ HPH, sSumm,	fHidden,tSUMM,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrSummary, NULL },
	{ 'v', sVers,	fNone,	tVERS,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrVersion, NULL },
	{ 'h', sHelp,	fNone,	tHELP,	oOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrUsage, NULL }
};
const BYTE	Options::OptCount = ArrCnt(Options::List);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, ProgParam, true, "alignment in bam/bed format" }
};
const BYTE Options::UsageCount = ArrCnt(Options::Usages);

dostream dout;	// stream's duplicator

/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv, ProgParam);
	if (fileInd < 0)	return 1;		// wrong option
	int ret = 0;						// main() return code

	Chrom::SetCustomOption(oCHROM);
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* iName = FS::CheckedFileName(argv[fileInd]);	// input alignment
		const char* oName = Options::GetSVal(oOUTFILE);			// output

		if (Options::Assigned(oOUTFILE)
		&& (!dout.OpenFile(Options::GetFileName(oOUTFILE, iName, OutFileSuff))))
			Err(Err::FailOpenOFile).Throw();

		ChromSizes cSizes(Options::GetSVal(oGEN), NULL, true);
		Reads test(ProgParam, iName, cSizes,
			Obj::eInfo(Options::GetIVal(oINFO)),
			true, true, Options::GetBVal(oALARM),
			true, Options::GetIVal(oMINSCR));

		vAlign(cSizes, test);
	}
	catch (Err & e) { ret = 1;	cerr << e.what() << LF; }
	catch (exception & e) { ret = 1;	cerr << e.what() << LF; }
	//catch(...)			{ ret = 1;	cerr << "Unregistered error\n"; }

	timer.Stop();
	//#ifdef OS_Windows
	//	system("pause");
	//#endif
	return ret;
}

/************************ class vAlign ************************/

vAlign::vAlign(const ChromSizes& cSizes, Reads& reads) :
	_caseDiff(Options::GetBVal(oCCASE)),
	_maxScore(reads.MaxScore())
{
	Reads::cItemsIter rit, ritend;
	chrid	cID;
	//bool	isPE = reads.IsPE();

	Read::FixedLen = reads.ReadLen();
	_mismAccums.Reserve(Read::FixedLen + 1);
	for (Reads::cIter cit = reads.cBegin(); cit != reads.cEnd(); cit++) {
		cID = CID(cit);
		if (!Chrom::NoCustom() && Chrom::CustomID() != cID)	continue;
		dout << Chrom::AbbrName(cID) << LF;
		_preciseAccum.Reset();
		_mismAccums.Clear();
		RefSeq seq(cID, cSizes);
		ritend = reads.ItemsEnd(cit);
		for (rit = reads.ItemsBegin(cit); rit != ritend; rit++)
			if (rit->InitCID == cID)
				if (rit->Pos == rit->RecPos)	// second number in PE Read name is the end of fragment
					_preciseAccum.AddRead(rit->Score);
				else
					_mismAccums[VerifyRead(seq, chrlen(rit->RecPos), rit->Pos)].AddRead(rit->Score);
		PrintStats(cID, reads.ItemsCount(cID));
	}
}

// Gets count of mismatches for tested Read
//	@seq: chromosome sequence
//	@templPos: template Read start position 
//	@testPos: tested Read start position
//	return: count of testet Read's mismatches in comparison with template Read
readlen vAlign::VerifyRead(const RefSeq& seq, chrlen templPos, chrlen testPos)
{
	const char* pos1 = seq.Read(templPos);
	const char* pos2 = seq.Read(testPos);
	readlen cnt = 0;

	for (readlen i = 0; i < Read::FixedLen; i++)
		if (_caseDiff) {
			if (*(pos1 + i) != *(pos2 + i))	cnt++;
		}
		else if (toupper(*(pos1 + i)) != toupper(*(pos2 + i)))	cnt++;
	return cnt;
}

// Prints statistic for given chrom
//	@cID: chrom ID
//	@rCnt: total count of Reads for given chrom
void vAlign::PrintStats(chrid cID, size_t rCnt)
{
	dout << "mism\treadCnt\tquality\n";
	dout << "precise\t";	_preciseAccum.Print(_maxScore);
	//dout << "_lowScoreCnt = " << _lowScoreCnt << LF;
	size_t cnt = _preciseAccum.Count();		// count of Reads at given chrom
	for (readlen i = 0; i <= Read::FixedLen; i++) {
		dout << int(i) << TAB;
		_mismAccums[i].Print(_maxScore);
		cnt += _mismAccums[i].Count();
	}
	dout << "total reads per chrom " << Chrom::Mark(cID) << ":\t" << cnt << TAB
		<< sPercent(cnt, rCnt, 0, 0, false) //<< LF;
		<< TAB << rCnt << LF;
	cnt = rCnt - cnt;
	dout << "reads per different chroms:\t" << cnt << TAB << sPercent(cnt, rCnt, 0, 0, false) << LF;
	fflush(stdout);		// when called from a package
}
/************************ end of class vAlign ************************/
