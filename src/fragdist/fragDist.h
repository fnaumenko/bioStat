/**********************************************************
fragDist.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 18.11.2020
-------------------------
Provides main functionality
***********************************************************/

#pragma once

#include "Data.h"
#ifdef _WIN32
#include <unordered_map>
#else
#include <tr1/unordered_map>
using namespace std::tr1;
#endif

enum optValue {		// options id
	//oCHROM,
	oNO_DUPL,
	oINFO,
	oALARM,
	oINPUT,
	oNORM,
	oPR_DIST,
	oOUTFILE,
	oTIME,
	oSUMM,
	oVERSION,
	oHELP,
};

enum class Inp {
	FRAG,
	READ,
	UNDEF
};

// Base length distribution class 
class LenDist
{
protected:
	LenFreq	_freq;					// length frequency statistics

public:
	// Print actual  distribution
	inline void Print(LenFreq::eType type, bool callNorm, bool prDistr) { _freq.Print(dout, type, callNorm, prDistr); }
};

// 'FragDist' represents fragment's length frequency statistics ('fragment distribution')
class FragDist : LenDist
{
	unordered_map<ULLONG, Reads::cItemsIter> _waits;	// 'waiting list' - pair mate candidate's collection
	
	// Increments statistics
	//	@rit1: iter to the first read in a pair
	//	@rit2: iter to the second read in a pair
	void IncrStat(const Reads::cItemsIter& rit1, const Reads::cItemsIter& rit2) {
		_freq.AddLen(rit1->Strand ?
			rit2->Pos + rit2->Len - rit1->Pos :
			rit1->Pos + rit1->Len - rit2->Pos
		);
	}

	// Add read to statistics if its mate is waiting already, otherwhise add it to the waiting list
	//	@rit: read's iterator
	void Add(const Reads::cItemsIter& rit);

	// clear waiting list to count the next chrom
	inline void Clear() { _waits.clear(); }

public:
	// Constructor by alignment; the instance is initialized according to the reads
	//	@test: paired-end reads collection
	FragDist(Reads& test);

	// Print actual frequency distribution
	inline void PrintDist(bool prDistr) { Print(LenFreq::eType::UNDEF, false, prDistr); }
};

// 'ReadDist' represents Read's length frequency statistics ('Read distribution')
class ReadDist : LenDist
{
	ULONG _cnt = 0;	// count of Reads

public:
	// Constructor by FastQ file; the instance is initialized according to the reads
	ReadDist(FqFile& file) {
		for (; file.GetSequence(); _cnt++)	
			_freq.AddLen(file.ReadLength());
	}

	ReadDist(DataInFile& file) {
		for (; file.GetNextItem(); _cnt++)
			_freq.AddLen(file.ItemLength());
	}

	// Print actual frequency distribution
	inline void PrintDist(bool callNorm, bool prDistr) { 
		dout << _cnt << " items\n";
		Print(LenFreq::eType::NORM, callNorm, prDistr);
	}
};