/**********************************************************
vAlign.h
Provides option emum and main functionality
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 07/28/2024
***********************************************************/
#pragma once

#include "Options.h"
#include "ChromSeq.h"
#include "DataReader.h"

enum optValue {
	oGEN,
	oCHROM,
	oMIN_SCORE,
	oCCASE,
	oDOUT_FILE,
	oLOCALE,
	oVERBOSE,
	oTIME,
	oSUMM,
	oVERSION,
	oHELP
};

enum class eVerb {	// verbose level
	TOT,	// total laconic output
	LAC,	// laconic for each chrom output
	DET		// detailed for each chrom output
};

// 'vAlign' represents Read's mismatches frequency and its output
class vAlign
{
	// Statistics accumulater
	class Stat
	{
		// 'ReadAccum' - accumulator of Read's count and average score 
		class ReadAccum
		{
			size_t	_count = 0;		// count of Reads
			double	_avrScore = 0;	// average score

		public:
			// Gets count of mismatches
			size_t Count() const { return _count; }

			void Clear() { _avrScore = 0, _count = 0; }

			// Add Read's score
			//	@param score: Read's score
			void AddRead(float score);

			// Adds chrom read accum to total one
			//	@param rAcc: added Read's accumulator
			void Add(const ReadAccum& rAcc);

			void Print(float maxScore) const {
				dout << _count /*<< TAB << (_avrScore / maxScore)*/ << LF;
			}
		};

		float	_maxScore = 0;
		size_t	_lowScoreCnt = 0;
		size_t	_duplCnt;
		ReadAccum	_preciseAccum;			// accumulator for exactly matched Reads
		map<readlen,ReadAccum>	_mismAccum;	// mismatches accumulator with count of mismatches as a key

	public:
		void SetMaxScore(float score) { if (score > _maxScore)	_maxScore = score; }

		// Increments count of too low scored reads
		void IncrLowScoreCnt() { _lowScoreCnt++; }

		// Increments count of precisely mapped reads
		//	@score: read's score
		void IncrPrecise(float score) { _preciseAccum.AddRead(score); }

		// Increments count of with mismatches mapped reads
		//	@mCnt: number of mismatches
		//	@score: read's score
		void IncrMism(readlen mCnt, float score) { _mismAccum[mCnt].AddRead(score); }

		// Adds chrom statistisc to total one
		void Add(const Stat& stat, readlen rLen);

		void Clear() { _lowScoreCnt = 0; _preciseAccum.Clear(); _mismAccum.clear(); }

		// Prints statistic for chrom
		//	@param cID: chrom ID
		//	@param rCnt: total count of Reads for given chrom
		//	@param prMismDist: if TRUE then print mismatches distribution
		void Print(chrid cID, size_t cnt, size_t duplCnt, bool prMismDist) const;
	};

	const bool _caseDiff;		// differ uppercase and lowercase characters
	const bool	_multy;			// if TRUE then output of statistics not only for one chrom
	const eVerb	_verb;
	const float _minScore;		// min possible score
	chrid	_cID;				// current chrom ID
	unique_ptr<ChromSeq> _seq;	// ref sequence
	const ChromSizes&	 _cs;
	RBedReader* _file;			// valid in constructor only!
	Stat	_chrStat;			// current chrom statistics
	Stat	_totStat;			// total statistics

	// Gets count of mismatches for tested Read
	//	@seq: chromosome sequence
	//	@r: tested Read
	//	@duplCnt: number of duplicates
	//	return: count of testet Read's mismatches in comparison with template pattern
	readlen VerifyRead(const ChromSeq& seq, const Read& r);

	void CloseChromStat(chrid cID, size_t cnt, size_t duplCnt)
	{
		if (_verb >= eVerb::LAC || _file->ReadedChromCount() == 1)
			_chrStat.Print(cID, ULONG(cnt), duplCnt, _verb == eVerb::DET);
		if (_multy)		_totStat.Add(_chrStat, _file->ReadLength());
	}

public:
	vAlign(const char* fname, ChromSizes& cs) :
		_caseDiff(Options::GetBVal(oCCASE)),
		_multy(Chrom::UserCID() == Chrom::UnID),
		_verb(eVerb(Options::GetIVal(oVERBOSE))),
		_minScore(Options::GetFVal(oMIN_SCORE)),
		_cs(cs)
	{
		RBedReader file(fname, &cs, BYTE_UNDEF, eOInfo::NM);
		_file = &file;
		file.Pass(*this);
		_file = nullptr;
	}

	// Treats current read
	bool operator()();

	// Closes current chrom, open next one
	//	@param cID: current chrom ID
	//	@param cLen: current chrom length
	//	@param cnt: current chrom items count
	//	@param nextcID: next chrom ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID);

	// Closes last chrom
	//	@param cID: last chrom ID
	//	@param cLen: current chrom length
	//	@param cnt: last chrom items count
	//	@param tCnt: total items count
	void operator()(chrid cID, chrlen cLen, size_t cnt, size_t tCnt);
};