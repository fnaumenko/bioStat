/**********************************************************
DefRegions.h  2023 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 11/21/2023
-------------------------
***********************************************************/
#pragma once

#include "Data.h"

// 'DefRegions' represents chrom's defined regions,
// initialized from ChromSizes (statically, at once)
// or from .fa files (dynamically, by request).
class DefRegions : public Chroms<Regions>
{
	ChromSizes& _cSizes;
	const chrlen	_minGapLen;	// minimal allowed length of gap
#ifdef _BIOCC
	const bool		_singleRgn = true;	// true if this instance has single Region for each chromosome
#endif


public:
	// Creates an instance by genome name, from chrom sizes file or genome.
	//	@param cSizes: chrom sizes
	//	@param minGapLen: minimal length which defines gap as a real gap
	DefRegions(ChromSizes& cSizes, chrlen minGapLen)
		: _cSizes(cSizes), _minGapLen(minGapLen) { Init(); }

	void Init();

	// Returns true if regions are not initialized
	bool IsEmpty() const { return !Count(); }

	ChromSizes& ChrSizes() { return _cSizes; }

	// Returns chrom's size by chrom's iterator
	chrlen Size(cIter it)	const { return Data(it).LastEnd(); }

	// Returns chrom's size by chrom's ID
	chrlen Size(chrid cID)	const { return At(cID).Data.LastEnd(); }

	// Gets chrom regions by chrom ID; lazy for real chrom regions
	const Regions& operator[] (chrid cID);

	// Copying constructor: creates empty copy!
	DefRegions(const DefRegions& gRgns) :
		_cSizes(gRgns._cSizes), _minGapLen(gRgns._minGapLen), _singleRgn(gRgns._singleRgn) {}

	// Returns miminal size of chromosome: for represented chromosomes only
	chrlen MinSize() const;

	// Returns total genome's size: for represented chromosomes only
	genlen GenSize() const;

	// Adds chromosome and regions without check up
	void AddChrom(chrid cID, const Regions& rgns) { AddVal(cID, rgns); }

#ifdef MY_DEBUG
	void Print() const;
#endif
};
