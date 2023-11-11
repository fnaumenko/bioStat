/************************************************************************************
vAlign is a fast verifier of reads of aligned DNA sequence,
recalled from initial artificial FastQ sequence.
It compares the original and actual coordinates of each read
and prints statistics of right and wrong mappings.

Copyright (C) 2017 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 11.12.2021
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
const string Product::Version = "2.0";
const string Product::Descr = "alignment verifier";

const char* ProgParam = "<in-file>";	// program parameter tip
const string OutFileSuff = "_valign.txt";
const string HelpOutFile = sFileDuplBegin + string(ProgParam) + OutFileSuff + sFileDuplEnd;

// *** Options definition
const char* verbs[] = { "TOT", "LAC","DET" };	// verbose option; corresponds to Inp

// --info option: types of info notations
enum tOptGroup { gTREAT, gOUTPUT, gOTHER };	// gOTHER should be the last 
const char* Options::OptGroups[] = { "Processing", "Output", "Other" };
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

const BYTE Options::Option::IndentInTabs = 3;

//	{ char,	str,	Signs,	type,	group,	defVal,	minVal,	maxVal,	strVal,	descr, addDescr }
// field 7: vUNDEF if value is prohibited
// field 6: vUNDEF if no default value should be printed
Options::Option Options::List[] = {
	{ 'g', sGen,	tOpt::OBLIG, tNAME,	gTREAT, vUNDEF, 0, 0, NULL,
	"reference genome library or single nucleotide sequence.", NULL },
	{ 'c',Chrom::Abbr,tOpt::NONE,tNAME,	gTREAT, vUNDEF, 0, 0, NULL,
	"treat specified chromosome only. For reference genome only", NULL },
	{ HPH,"min-scr",  tOpt::NONE,tINT,	gTREAT, 0, 0, 1000, NULL, "score threshold for treated reads", NULL },
	{ HPH,"char-case",tOpt::NONE,tENUM,	gTREAT, FALSE,	0, 2, (char*)Options::Booleans,
	"recognize uppercase and lowercase characters in template and test\nas different", NULL },
	{ 'o', sOutput,	tOpt::OBLIG,tNAME,	gOUTPUT,NO_DEF,	0,	0, NULL, HelpOutFile.c_str(), NULL },
	{ 'T', "sep",	tOpt::NONE,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL, "use 1000 separator in output", NULL },
	{ 'V', "verbose",tOpt::NONE, tENUM,	gOUTPUT, float(eVerb::LAC), float(eVerb::TOT), ArrCnt(verbs), (char*)verbs,
	 "set output verbose level:\n? - only total detailed,\n? - laconic for each chromosome and total detailed,\n? - detailed for each chromosome", NULL },
	{ 't', sTime,	tOpt::NONE,	tENUM,	gOTHER,	FALSE,	NO_VAL, 2, NULL, sPrTime,	NULL },
	{ HPH, sSumm,	tOpt::HIDDEN,tSUMM,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrSummary,NULL },
	{ 'v', sVers,	tOpt::NONE,	tVERS,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrVersion,NULL },
	{ 'h', sHelp,	tOpt::NONE,	tHELP,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrUsage,	NULL }
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

	if (Options::GetBVal(oLOCALE))	cout.imbue(locale(LOCALE_ENG));

	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* iName = FS::CheckedFileName(argv[fileInd]);	// input alignment
		const char* oName = Options::GetSVal(oOUTFILE);			// output

		if (Options::Assigned(oOUTFILE)
		&& (!dout.OpenFile(Options::GetFileName(oOUTFILE, oName, OutFileSuff))))
			Err(Err::FailOpenOFile).Throw();

		dout << iName << LF;	cout.flush();
		ChromSizes cSizes(Options::GetSVal(oGEN), oCHROM, true);
		vAlign(iName, cSizes);
	}
	catch (Err & e) { ret = 1;	cerr << e.what() << LF; }
	catch (exception & e) { ret = 1;	cerr << e.what() << LF; }

	timer.Stop(0, false, true);
	//#ifdef OS_Windows
	//	system("pause");
	//#endif
	return ret;
}

/************************ class vAlign ************************/

// Adds chrom statistisc to total one
void vAlign::Stat::Add(const Stat& stat, readlen rLen)
{
	_lowScoreCnt += stat._lowScoreCnt;
	_preciseAccum.Add(stat._preciseAccum);
	SetMaxScore(stat._maxScore);
	for (readlen i = 0; i <= rLen; i++) {
		const bool found = stat._mismAccum.find(i) != stat._mismAccum.end();
		if (_mismAccum.find(i) == _mismAccum.end()) {
			if (found)
				_mismAccum[i] = stat._mismAccum.at(i);
		}
		else if (found)
			_mismAccum[i].Add(stat._mismAccum.at(i));
	}
}

// Prints count and percentage of total
void PrintCount(const string& title, int sWidth, size_t cnt, ULONG totalCnt, int dWidth)
{
	if(cnt)
		dout<< left << setw(sWidth) << title
			<< right << setw(dWidth) << cnt << SPACE << SPACE << sPercent(ULONG(cnt), totalCnt, 4, 0, false) << LF;
}

// Prints statistic for given chrom
//	@cID: chrom ID
//	@cnt: total count of Reads for given chrom
//	@duplCnt: number of duplicates for given chrom
//	@prMismDist: if TRUE then print mismatches distribution
void vAlign::Stat::Print(chrid cID, ULONG cnt, size_t duplCnt, bool prMismDist) const
{
	const bool isTotal = cID == Chrom::UnID;

	if (isTotal)		dout << "TOTAL\n";
	int wd;			// wigth of digital field
	// *** print mismathes distribution
	ULONG rCnt = _preciseAccum.Count();
	ULONG rPrecCnt = rCnt;
	if (prMismDist) {
		//dout << "mismCnt\treadCnt\tAvrQual\n";
		//wd = 3 * 8;		// wigth of digital field
		dout << "mismCnt\treadCnt\n";
		wd = 2 * 8;		// wigth of digital field
		PrintHorLine(wd);
		dout << "precise\t";	_preciseAccum.Print(_maxScore);
	}
	for (const auto& acc : _mismAccum) {
		if (prMismDist) {
			dout << acc.first << TAB;
			acc.second.Print(_maxScore);
		}
		if (!acc.first)	rPrecCnt += acc.second.Count();
		rCnt += acc.second.Count();
	}
	if (prMismDist)	PrintHorLine(wd);

	// *** print stats
	wd = DigitsCount(cnt, Options::GetBVal(oLOCALE));
	const string reads = "reads ";
	const string mapped = "mapped ";
	const string misms = " mismathes:";
	const string title = reads + mapped + " to the correct chrom";
	const int ws = int(title.length() + 3 + (isTotal ? 0 : Chrom::MarkLength(cID)));

	// total reads
	dout << left << setw(ws) << (reads + "total" + (isTotal ? strEmpty : (" per chrom " + Chrom::Mark(cID))) + COLON);
	dout << right << setw(wd) << cnt;
	if (duplCnt)
		dout << "  (including " << duplCnt << sPercent(ULONG(duplCnt), cnt, 3, 0, true) << " duplicates)";
	dout << LF;
	// reads discarded due to low score
	PrintCount(reads + "discarded due to low score:", ws, _lowScoreCnt, cnt, wd);
	// reads mapped correctly to the chrom
	PrintCount(title + (isTotal ? strEmpty : (SPACE + Chrom::Mark(cID))) + COLON, ws, rCnt, cnt, wd);
	dout << "from wich:\n";
	// reads mapped reads mapped without mismathes
	PrintCount(string(reads.length(), SPACE) + mapped + "without" + misms, ws, rPrecCnt, cnt, wd);
	// reads mapped reads mapped with mismathes
	PrintCount(string(reads.length(), SPACE) + mapped + "with" + misms, ws, rCnt - rPrecCnt, cnt, wd);
	// reads mapped to different chroms
	PrintCount(reads + mapped + "to different chroms:", ws, cnt - rCnt - _lowScoreCnt, cnt, wd);
	fflush(stdout);				// when called from a package
}

// Gets count of mismatches for tested Read
//	@seq: chromosome sequence
//	@r: tested Read 
//	return: count of testet Read's mismatches in comparison with template pattern
readlen vAlign::VerifyRead(const RefSeq& seq, const Read& r)
{
	const char* pos1 = seq.Seq(r.RecPos) - 1;
	const char* pos2 = seq.Seq(r.Pos) - 1;
	readlen cnt = 0;

	if (_caseDiff)
		for (readlen i = 0; i < r.Len; i++)
			cnt += *++pos1 != *++pos2;
	else
		for (readlen i = 0; i < r.Len; i++)
			cnt += toupper(*++pos1) != toupper(*++pos2);

	return cnt;
}

// Treats current read
bool vAlign::operator()()
{
	Read r(*_file);
	if (r.Score < _minScore)	_chrStat.IncrLowScoreCnt();
	else {
		_chrStat.SetMaxScore(r.Score);
		if (r.RecCID == _cID)			// mapped to the same chrom
			if (r.Pos == r.RecPos)		// second number in PE Read name is the end of fragment
				_chrStat.IncrPrecise(r.Score);
			else
				_chrStat.IncrMism(VerifyRead(*_seq, r), r.Score);
	}
	return true;
}

// Closes current chrom, open next one
//	@cID: current chrom ID
//	@cnt: current chrom items count
//	@nextcID: next chrom ID
void vAlign::operator()(chrid cID, chrlen, size_t cnt, chrid nextcID)
{
	if (cID != Chrom::UnID)	// not pre-first chrom
		CloseChromStat(cID, cnt, _file->DuplCount());
	if (_verb >= eVerb::LAC) { dout << Chrom::ShortName(nextcID) << LF;	fflush(stdout); }
	_seq.reset(new RefSeq(_cID = nextcID, _cs));
	_chrStat.Clear();
}

// Closes last chrom
//	@cID: last chrom ID
//	@cnt: last chrom items count
//	@tCnt: total items count
void vAlign::operator()(chrid cID, chrlen, size_t cCnt, ULONG tCnt)
{
	CloseChromStat(cID, cCnt, _file->DuplCount());
	if (_file->ReadedChromCount() > 1)
		_totStat.Print(Chrom::UnID, tCnt, _file->DuplTotalCount(), _verb >= eVerb::TOT);
}

/************************ end of class vAlign ************************/
