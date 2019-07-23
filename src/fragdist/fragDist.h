/**********************************************************
fragDist.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 18.06.2019
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


enum optValue {
	oCHROM,
	oNO_DUPL,
	oINFO,
	oALARM,
	oPR_DIST,
	oOUTFILE,
	oTIME,
	oVERSION,
	oSUMM,
	oHELP
};


class FragDist
{
	struct Waiting {
		Reads::cItemsIter It;
		bool	Free;

		inline Waiting(Reads::cItemsIter& it) : It(it), Free(false) {}

		inline void Replace(Reads::cItemsIter& it) { It = it; Free = false; }
	};

	FragFreq	_freq;
	unordered_map<chrlen, Reads::cItemsIter> _waits;
	
	// Increments statistics
	//	@rit1: first read in pair
	//	@rit2: second read in pair
	void IncrStat(const Reads::cItemsIter& rit1, const Reads::cItemsIter& rit2) {
		// get frag len
		fraglen len = rit1->Strand ?
			rit2->Pos + rit2->Len - rit1->Pos:
			rit1->Pos + rit1->Len - rit2->Pos;
		_freq.AddFrag(len);
	}

	// Check read mate in waiting list, add to statistics if founded
	//	@rit: read iterator
	void Add(const Reads::cItemsIter& rit) {
		Reads::cItemsIter rit0 = _waits[rit->Numb];	// iterator to pair mate
		//if(rit0._Ptr) {		// item exists; doesn't compiled with -std=c++11
		if(&(*rit0)) {			// item exists
			if(rit0->Strand != rit->Strand) {
				IncrStat(rit0, rit);
				_waits.erase(rit->Numb);
			}
			//else	cout << rit->Numb << ": repeated strand\n";
		}
		else	
			_waits[rit->Numb] = rit;	// add item
	}

	// clear counter before new chrom
	inline void Clear() { _waits.clear(); }

public:
	FragDist(Reads& test);

	// Print actual frequency distribution
	inline void PrintDist(bool prDistr) const { _freq.Print(dout, prDistr); }
};
