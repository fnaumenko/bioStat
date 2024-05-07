/**********************************************************
callDist.h (c) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 05/07/2024
-------------------------
Provides main functionality
***********************************************************/

#pragma once

#include "Options.h"
#include "ChromData.h"
#include "Distrib.h"
#include "FqReader.h"

enum optValue {		// options id
	oINPUT,
	oCHROM,
	oDTYPE,
	oDUPL,
	oPR_DIST,
	oPR_STATS,
	oOUTFILE,
	oTIME,
	oSUMM,
	oVERSION,
	oHELP,
};

// Input data type
enum class InpType { FRAG, READ };

// Base length distribution class 
class LenDist
{
	Distrib	_freq;					// length frequency statistics
	RBedReader* _file = nullptr;	// valid in constructor only!

protected:
	// pass through file records
	template<typename T>
	void Pass(T* obj, RBedReader& file) {
		_file = &file;
		file.Pass(*obj);
		_file = nullptr;
	}

	// Gets the file being read
	inline const RBedReader& File() { return *_file; }

	// Adds frag length to the frequency distribution
	inline void AddLen(fraglen len) { _freq.AddVal(len); }

public:
	// Print actual frequency distribution on a new line
	//	@param ctype: combined type of distribution
	//	@param prDistr: if true then print distribution additionally
	void Print(Distrib::eCType ctype, bool prDistr)
	{
		// empty input is checked already in the 'UniBedReader' constructor
		if (_freq.Size())
			_freq.Print(dout, ctype, prDistr);
	}
};

// 'FragDist' represents fragment's length frequency statistics ('fragment distribution')
class FragDist : public LenDist
{
	FragIdent _fIdent;
	bool	_dupl;						// if TRUE if duplicate frags are allowed
	bool	_checkedPE = false;			// if TRUE if reads have been checked for PE

public:
	FragDist(const char* fname, bool prStats) : _fIdent(_dupl = Options::GetBVal(oDUPL))
	{
		vector<UniBedReader::Issue> issues = { "duplicates" };
		RBedReader file(fname, nullptr, BYTE_UNDEF, eOInfo::NONE, false);

		Pass(this, file);

		size_t cnt = _fIdent.Count();
		issues[0].Cnt = _fIdent.DuplCount();
		UniBedReader::PrintItemCount(cnt, "fragments");
		if (issues[0].Cnt) {
			if (_dupl)	issues[0].Action = UniBedReader::eAction::ACCEPT;
			UniBedReader::PrintStats(cnt, issues[0].Cnt, issues, prStats);
		}
	}

	// treats current read
	bool operator()();

	// Closes current chrom, open next one
	inline void operator()(chrid, chrlen, size_t, chrid) {}

	// Closes last chrom
	inline void operator()(chrid, chrlen, size_t, size_t) {
#ifdef MY_DEBUG
		printf(">>max size: %d ", _maxSize);
#endif
	}
};

// 'ReadDist' represents Read's length frequency statistics ('Read distribution')
class ReadDist : public LenDist
{
public:
	// Constructor by BAM/BED file
	ReadDist(const char* fname, bool prStats) {
		RBedReader file(fname, nullptr,
			BYTE(Options::GetRDuplLevel(oDUPL)),
			prStats ? eOInfo::STAT : eOInfo::STD,
			false);

		Pass(this, file);
	}

	// treats current read
	inline bool operator()() { AddLen(File().ItemRegion().Length()); return true; }

	// Closes current chrom, open next one
	inline void operator()(chrid, chrlen, size_t, chrid) {}

	// Closes last chrom
	inline void operator()(chrid, chrlen, size_t, size_t) {}
};

// 'FqReadDist' represents 'row' (fastq) Read's length frequency statistics ('Read distribution')
class FqReadDist : public LenDist
{
public:
	// Constructor by FastQ file
	FqReadDist(const char* fname)  {
		ULONG	cnt = 0;		// count of reads

		for (FqReader file(fname); file.GetSequence(); cnt++)
			AddLen(file.ReadLength());
		UniBedReader::PrintItemCount(cnt, FT::ItemTitle(FT::eType::FQ, cnt > 1));
		dout << LF;
	}
};