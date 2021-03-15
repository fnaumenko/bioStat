/**********************************************************
callDist.h (c) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 14.03.2021
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
	oDTYPE,
	//oNORM,
	oPR_DIST,
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
	static const char* ItemTitles[];
	int	_itemInd;		// index of item title

protected:
	LenFreq	_freq;		// length frequency statistics
	ULONG _cnt = 0;		// count of items

	inline LenDist(InpType iType) : _itemInd(int(iType)) {}

public:
	// Print actual frequency distribution
	//	@ctype: combined type of distribution
	//	@prDistr: if true then print distribution additionally
	void Print(LenFreq::eCType ctype, bool prDistr) {
		dout << _cnt << SPACE << ItemTitles[_itemInd] << LF;
		_freq.Print(dout, ctype, prDistr);
	}

};

// 'FragDist' represents fragment's length frequency statistics ('fragment distribution')
class FragDist : public LenDist
{
	unordered_map<ULLONG, Reads::cItemsIter> _waits;	// 'waiting list' - pair mate candidate's collection

	// Adds frag to the freq distribution
	//	@rit1: iter to the first read in a pair
	//	@rit2: iter to the second read in a pair
	void AddFrag(const Reads::cItemsIter& rit1, const Reads::cItemsIter& rit2) {
		_cnt++;
		_freq.AddLen(rit1->Strand ?
			rit2->Pos + rit2->Len - rit1->Pos :
			rit1->Pos + rit1->Len - rit2->Pos
		);
	}

	// Add read to statistics if its mate is waiting already, otherwhise put it on the waiting list
	//	@rit: read's iterator
	void PutRead(const Reads::cItemsIter& rit);

	// clear waiting list to count the next chrom
	inline void Clear() { _waits.clear(); }

public:
	// Constructor by alignment; the instance is initialized according to the reads
	//	@test: paired-end reads collection
	FragDist(Reads& test);
};

// 'ReadDist' represents Read's length frequency statistics ('Read distribution')
class ReadDist : public LenDist
{
public:
	// Constructor by FastQ file
	ReadDist(FqFile& file) : LenDist(InpType::READ) {
		for (; file.GetSequence(); _cnt++)
			_freq.AddLen(file.ReadLength());
	}

	// Constructor by BAM/BED file
	ReadDist(DataInFile& file) : LenDist(InpType::READ) {
		for (; file.GetNextItem(); _cnt++)
			_freq.AddLen(file.ItemLength());
	}
};