/**********************************************************
PDTest - Peak Detectors test
Copyright (C) 2024 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07/19/2024
-------------------------
***********************************************************/

#pragma once
#include "Features.h"

enum optValue {		// options id
	//oGEN,
	//oCHROM,
	oTEMPL,
	oMIN_SCORE,
	oALARM,
	oOUTFILE,
	oTIME,
	oVERSION,
	oSUMM,
	oHELP,
};

const fraglen ROI_ext = 500;

// 'IGVlocus' is designed to print locus that can be pasted into the address bar of Integrative Genomics Viewer
class IGVlocus
{
	string _chrom = Chrom::AbbrName(0);
	mutable char _buf[2 * 10 + 5 + 2]{ 0 };	// 2 * max position length + Chrom::MaxAbbrNameLength + 2 separators

	// Prints IGV locus to inner buffer
	void NPrint(chrlen start, chrlen end) const
	{
		sprintf(_buf, "%s:%d-%d", _chrom.c_str(), start - ROI_ext, end + ROI_ext);
	}

public:
	//IGVlocus(chrid cID) : _chrom(Chrom::AbbrName(cID)) {}

	void SetChrom(chrid cID) { _chrom = Chrom::AbbrName(cID); }

	const char* ChromAbbrName() const { return _chrom.c_str(); }

	// Prints IGV locus to inner buffer
	//	@param start: feature's start position
	//	@param end: feature's end position
	//	@returns: inner buffer
	const char* Print(chrlen start, chrlen end) const { NPrint(start, end); return _buf; }

	// Prints IGV locus to inner buffer
	//	@param pos: feature's start position
	//	@returns: inner buffer
	const char* Print(chrlen pos) const { return Print(pos, pos); }
};

class TxtOutFile
{
	FILE* _file;
public:
	TxtOutFile(const char* name) { _file = fopen(name, "w"); }
	~TxtOutFile() { fclose(_file); }

	template<typename... Args>
	void Write(const char* format, Args ... args) { fprintf(_file, format, args...); }
};


// Statistical Binary Classifiers
static class BC
{
	static const char* titles[2];

public:
	enum eBC {
		FP,	// false positive error
		FN,	// false negative error
	};

	static const char* Title(eBC bc) { return titles[bc]; }
} bc;


// Features wrapper for manipulating iterators for a given chromosome 
class ChromFeaturesIters
{
public:
	using iterator = Features::cItemsIter;

	static void SetChrom(chrid cID)	{ if (locus)	locus->SetChrom(cID); }

	ChromFeaturesIters(BC::eBC bc, const Features& fs, Features::cIter cIt, float minScore = 0)
		: _bc(bc), _minScore(minScore)
	{
		auto& data = fs.Data(cIt);
		_beginII = fs.ItemsBegin(data);
		_endII = fs.ItemsEnd(data);
		if (!locus)	locus = new IGVlocus;
	}

	~ChromFeaturesIters() {
		if (oFile) { delete oFile; oFile = nullptr;	}
		if (locus) { delete locus; locus = nullptr;	}
	}

	static void SetOutFile(const char* fname)
	{
		if (fname) {
			oFile = new TxtOutFile(fname);
			oFile->Write("chrom\tBC start\tend  \tscore\tlocus\n");
		}
	}

	void PrBcCount(char sep) const { cout << BC::Title(_bc) << SepCl << _bcCount << sep; }

	iterator& begin()	{ return _beginII; }
	iterator& end()		{ return _endII; }

	static chrlen	Start(iterator it) { return it->Start; }
	static chrlen	End(iterator it) { return it->End; }
	static bool		IsWeak(iterator it) { return false; }
	void			Discard(iterator it) { 
		if (it->Value < _minScore)	return;
		_bcCount++;
		if (oFile) {
			float score = it->Value;
			if (score > 1)	score /= 1000;
			oFile->Write("%s\t%-3s%d\t%d\t%.2f\t%s\n",
				locus->ChromAbbrName(),
				BC::Title(_bc),
				it->Start,
				it->End,
				score,
				locus->Print(it->Start, it->End)
			);
		}
	}

private:
	static IGVlocus* locus;
	static TxtOutFile* oFile;

	BC::eBC	 _bc;
	USHORT	 _bcCount = 0;
	float	 _minScore = 0;
	iterator _beginII;
	iterator _endII;
};


template<typename T>
void DiscardNonOverlapRegions(T rgns[2], fraglen minOverlapLen)
{
	const typename T::iterator itEnd[2]{ rgns[0].end(), rgns[1].end() };
	typename T::iterator it[2]{ rgns[0].begin(), rgns[1].begin() };
	BYTE s;		// index of the left started region: 0 - direct, 1 - reverse
	BYTE suspended = 0;	// if 1st|2nd bit is set then direct|reverse region is suspended and should be analyzed in the next pass

	auto start = [&it](BYTE s) -> chrlen { return T::Start(it[s]); };
	auto end = [&it](BYTE s) -> chrlen { return T::End(it[s]); };

	// unconditional discarde the region
	auto discardeRgn = [&](BYTE s) {
		//T::Discard(it[s], s);
		rgns[s].Discard(it[s]);
		suspended &= ~(1 << s);	// reset suspended s
		it[s]++;
	};
	// unconditional discarde the regions; always applied to the left region
	auto discardeLastRgns = [&it, &itEnd, &discardeRgn](BYTE s) {
		while (it[s] != itEnd[s])
			discardeRgn(s);
	};

	while (it[0] != itEnd[0] && it[1] != itEnd[1]) {
		s = start(0) > start(1);
		if (end(s) > start(!s) + minOverlapLen) {	// strong intersection
			/***************************************
				-----    ?????	more regions?		left overlaps
			+++++++++++++		suspended...		left overlaps
				---------		suspended...		right overlaps
			+++++++	  ?????		more regions?		right overlaps
			***************************************/
			bool valid = true;
			if (T::IsWeak(it[s]))	discardeRgn(s), valid = false;
			if (T::IsWeak(it[!s]))	discardeRgn(!s), valid = false;
			if (valid) {
				s ^= end(s) < end(!s);	// flip by condition
				it[!s]++;
				suspended |= 1 << s;	// set suspended s
			}
		}
		else {							// weak intersection or no one
			// conditional close of the region
			if (suspended)	suspended = 0;
			else			
				//T::Discard(it[s], s);
				rgns[s].Discard(it[s]);

			// close remaining complementary regions
			if (++it[s] == rgns[s].end()) {
				discardeLastRgns(!s);
				break;			// no need to check in while()
			}
		}
	}

	// close 'out of scope' regions
	if ((s = it[1] != itEnd[1]) || it[0] != itEnd[0]) {
		it[s]++;
		discardeLastRgns(s);
	}
}