/**********************************************************
callDist.h (c) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 02/05/25
-------------------------
Provides main functionality
***********************************************************/

#pragma once

#include "Options.h"
#include "ChromData.h"
#include "Distrib.h"
#include "DataReader.h"
#include "FqReader.h"

enum optValue {		// options id
	oINPUT,
	oCHROM,
	oDTYPE,
	oDUPL,
	oPR_DIST,
	oPR_STATS,
	oDOUT_FILE,
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
	void Pass(T* obj, RBedReader& file, bool prLF = true) {
		_file = &file;
		file.Pass(*obj, prLF);
		_file = nullptr;
	}

	// Gets the file being read
	const RBedReader& File() { return *_file; }

	// Adds frag length to the frequency distribution
	void AddLen(fraglen len) { _freq.IncrFreq(len); }

public:
	// Print actual frequency distribution on a new line
	//	@param dtype: combined type of distribution
	//	@param prDistr: if true then print distribution additionally
	void Print(Distrib::eDType dtype, bool prDistr)
	{
		// empty input is checked already in the 'UniBedReader' constructor
		if (_freq.Size()) {
			_freq.CalcADParams(dtype);
			_freq.Print(dout, true, prDistr);
		}
	}
};

// 'FragDist' represents fragment's length frequency statistics ('fragment distribution')
class FragDist : public LenDist
{
	FragIdent _fIdent;
	bool	_duplAccept;			// if TRUE if duplicate frags are allowed

public:
	FragDist(const char* fname, bool prStats) : _fIdent(_duplAccept = Options::GetBVal(oDUPL))
	{
		RBedReader file(
			fname,
			nullptr,	// chrom sizes
			0,			// no read duplicates control; fragment duplicates is controlled by _fIdent
			eOInfo::STD,
			true,		// abort invalid
			true		// first line will be pre-read
		);

		// pre-read first item to check for PE sequence
		file.GetNextItem();		// no need to check for empty sequence
		if (!file.IsPaired())
			Err("only paired-end reads are acceptable to call fragments distribution",
				file.CondFileName()).Throw();

		Pass(this, file, false);

		// print statistic
		size_t cnt = _fIdent.Count();
		dout << SepSCl << cnt << sPercent(cnt, file.ReadedItemCount()/2, 2, 0, true) << " identified fragments";

		if (prStats) {
			vector<UniBedReader::Issue> issues = { "unacceptably short","duplicates" };

			issues[0].Cnt = _fIdent.ShortCount();
			issues[1].Cnt = _fIdent.DuplCount();
			if (_duplAccept)	issues[1].Action = UniBedReader::eAction::ACCEPT;
			file.PrintStats(cnt, issues);
		}
		dout << LF;
	}

	// treats current read
	bool operator()(bool)
	{
		Region frag;
		const Read read(File());

		if (_fIdent(read, File().ReadLength(), frag))
			AddLen(frag.Length());
		return true;
	}

	// Closes current chrom, open next one
	void operator()(chrid, chrlen, size_t, chrid) {}

	// Closes last chrom
	void operator()(chrid, chrlen, size_t, size_t) {
#ifdef MY_DEBUG
		printf(">>max size: %d ", _fIdent.MaxMapSize());
#endif
	}
};

// 'ReadDist' represents Read's length frequency statistics ('Read distribution')
class ReadDist : public LenDist
{
public:
	// Constructor by BAM/BED file
	ReadDist(const char* fname, bool prStats) {
		RBedReader file(
			fname,
			nullptr,
			Options::GetRDuplPermit(oDUPL),
			prStats ? eOInfo::STAT : eOInfo::STD
		);

		Pass(this, file);
	}

	// treats current read
	bool operator()(bool) { AddLen(File().ItemRegion().Length()); return true; }

	// Closes current chrom, open next one
	void operator()(chrid, chrlen, size_t, chrid) {}

	// Closes last chrom
	void operator()(chrid, chrlen, size_t, size_t) {}
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