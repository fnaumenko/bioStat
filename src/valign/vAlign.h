/**********************************************************
vAlign.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 1.1.2020
-------------------------
Provides option emum and main functionality
***********************************************************/
#pragma once

enum optValue {
	oGEN,
	//oTEMPL,
	oCHROM,
	oMINSCR,
	oCCASE,
	oINFO,
	oALARM,
	oOUTFILE,
	oTIME,
	oSUMM,
	oVERSION,
	oHELP
};

class vAlign
/*
 * 'vAlign' represents Read's mismatches density
 */
{
private:
	struct ReadAccum
	/*
	 * 'ReadAccum' - accumulator of Read's count and average score 
	 */
	{
	private:
		size_t	_count;		// count of Reads
		double	_avrScore;	// average score

	public:
		ReadAccum() : _count(0), _avrScore(0) {}

		// Gets count of mismatches
		inline size_t	Count() const { return _count; }

		inline void	Reset() { _count = 0; _avrScore = 0.0; }

		// Add Read's score
		//	@score: Read's score
		void AddRead(float score) {
			_avrScore = _count ? 
				(_avrScore*(_count-1) + score)/_count:	// rolling average
				score;
			_count++;
		}

		void Print(float maxScore) {
			dout << _count;
			if( _count )	dout << TAB << (_avrScore/maxScore);
			dout << LF;
		}
	};

	bool _caseDiff;		// differ uppercase and lowercase characters
	float _maxScore;
	//int _lowScoreCnt;

	// container of mismatches Reads:
	// index - count of mismatches,
	// value by index - count of Reads and mean score by given amount of mismatches
	Array<ReadAccum> _mismAccums;
	ReadAccum	_preciseAccum;		// accumulator for exactly matched Reads

	// Gets count of mismatches for tested Read
	//	@seq: chromosome sequence
	//	@templPos: template Read start position 
	//	@testPos: tested Read start position
	//	return: count of testet Read's mismatches in comparison with template Read
	readlen VerifyRead(const RefSeq& seq, chrlen templPos, chrlen testPos);

public:
	vAlign(const ChromSizes& cSizes, Reads &reads);

	// Prints statistic for chrom
	//	@cID: chrom ID
	//	@rCnt: total count of Reads for given chrom
	void PrintStats(chrid cID, size_t rCnt);
};