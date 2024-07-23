/**********************************************************
PDTest - Peak Detectors test
Copyright (C) 2024 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07/23/2024
-------------------------
***********************************************************/

#pragma once
#include "Features.h"
#include <numeric>		// std::reduce

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
				rgns[0].Accept(it);
				s ^= end(s) < end(!s);	// flip by condition
				it[!s]++;
				suspended |= 1 << s;	// set suspended s
			}
		}
		else {							// weak intersection or no one
			// conditional close of the region
			if (suspended)	suspended = 0;
			else			rgns[s].Discard(it[s]);

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

// target set for computation
enum eTarget {
	CHR,	// for the current chromosome
	ALL		// for all
};

// Keeps data and calculates Standard Deviation
class StandDev
{
	vector<short> _devs;	// deviations
	UINT	_startInd = 0;	// starting index in _devs for the current chromosome
	int		_sumDev[2]{ 0,0 };	// sum of deviation for current chromosome, for all

public:
	void Reserve(size_t capacity) { _devs.reserve(capacity); }

	// Returns Standard Deviation
	//	@param t: target set for computation
	float GetSD(eTarget t) const
	{
		UINT startInd = !t * _startInd;		// _startInd for chrom, 0 for all
		chrlen cnt = chrlen(_devs.size()) - startInd;
		auto avr = float(_sumDev[t]) / cnt;
		float sumD = 0;

		for (auto it = _devs.begin() + startInd; it != _devs.end(); it++)
			sumD += float(pow(float(*it) - avr, 2));
		return sqrt(sumD / cnt);
	}

	// Adds deviation value
	void AddDev(short dev)
	{
		_devs.push_back(dev);
		_sumDev[CHR] += dev;
	}

	// Stops counting data for the current chromosome
	void ResetChrom()
	{
		_startInd = UINT(_devs.size() - 1);
		_sumDev[ALL] += _sumDev[CHR];
		_sumDev[CHR] = 0;
	}
};

// 'Features' wrapper for manipulating iterators for a given chromosome 
class FeaturesStatData
{
public:
	using iterator = Features::cItemsIter;

private:
	static IGVlocus*	locus;
	static TxtOutFile*	oFile;

	// prints binary classifiers
	static void PrBcCounts(const chrlen count[2])
	{
		const char* delim = "  ";
		cout << setw(4) << count[0] << delim << setprecision(3) << float(count[0]) / count[1] << delim;
	}

	const Features& _fs;
	StandDev&	_sd;
	iterator	_beginII;
	iterator	_endII;
	BC::eBC		_bc;				// needed for printing to dump file
	float		_minScore;
	/*
	* valid feature's counters (per chrom and total) are only needed 
	* for autonomous calculation of FNR and FDR (when printing).
	* In FeaturesStatTuple::PrintF1() we can use _sd.size() as True Positive
	*/
	chrlen	_cntBC[2][2]{ 
		{ 0,0 },	// false, total valid feature's count for chromosome
		{ 0,0 }		// false, total valid feature's count in total
	};

public:
	static void SetOutFile(const char* fname)
	{
		if (fname) {
			oFile = new TxtOutFile(fname);
			oFile->Write("chrom\tBC start\tend  \tscore\tlocus\n");
		}
	}

	FeaturesStatData(BC::eBC bc, const Features& fs, StandDev& sd, float minScore)
		: _bc(bc), _fs(fs), _sd(sd), _minScore(minScore)
	{
		if(oFile)
			locus = new IGVlocus;
	}

	~FeaturesStatData()
	{
		if (oFile) { delete oFile; oFile = nullptr; }
		if (locus) { delete locus; locus = nullptr; }
	}

	// returns binary classifiers
	//	@param t: target set for computation
	const chrlen* BC(eTarget t) const { return _cntBC[t]; }

	// prints binary classifiers
	//	@param t: target set for computation
	void PrintStat(eTarget t) const { PrBcCounts(_cntBC[t]); }

	// sets counting local stats data (for given chromosome)
	//	@param cIt: chromosome's iterator
	void SetChrom(Features::cIter cIt)
	{
		auto& data = _fs.Data(cIt);
		_cntBC[CHR][1] = chrlen(data.ItemsCount());
		_beginII = _fs.ItemsBegin(data);
		_endII = _fs.ItemsEnd(data);
		if (locus)	locus->SetChrom(CID(cIt));
	}

	// stops counting local stats data
	void ResetChrom()
	{
		_cntBC[ALL][0] += _cntBC[CHR][0];
		_cntBC[ALL][1] += _cntBC[CHR][1];
		memset(_cntBC[CHR], 0, sizeof(_cntBC) / 2);
	}

	iterator& begin()	{ return _beginII; }
	iterator& end()		{ return _endII; }

	static chrlen	Start	(iterator it) { return it->Start; }
	static chrlen	End		(iterator it) { return it->End; }
	static bool		IsWeak	(iterator it) { return false; }	// stub
	void Accept(const iterator it[2])
	{
		_sd.AddDev(short(int(it[0]->Centre()) - it[1]->Centre()));
	}
	void Discard	(iterator it)
	{ 
		if (it->Value < _minScore)
			//_cCntBC[1]--;		// decrease chrom counter
			_cntBC[CHR][1]--;		// decrease chrom counter
		else {
			//_cCntBC[0]++;		// increase false counter
			_cntBC[CHR][0]++;		// increase false counter
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
	}
};


class FeaturesStatTuple
{
	static const USHORT titleLineLen = 42;
	static void PrintFooter()
	{
		PrintSolidLine(titleLineLen);
		cout << sTotal << COLON << TAB;
	}

	StandDev _sd;
	FeaturesStatData _data[2];
	chrid _cID = Chrom::UnID;

	// returns F1 score
	//	@param t: target set for computation
	float GetF1(eTarget t) const
	{
		auto fn = _data[0].BC(t)[0];
		auto dtp = 2 * (_data[0].BC(t)[1] - fn);		// double TP
		return float(dtp) / (dtp + _data[1].BC(t)[0] + fn);
	}

	void PrintF1(eTarget t) const { cout << GetF1(t) << "  "; }

public:
	static void PrintHeader()
	{
		cout << setw(11) << BC::Title(BC::FN) << "   FNR   " << BC::Title(BC::FP) << "   FDR    F1     SD\n";
		PrintSolidLine(titleLineLen);
	}

	FeaturesStatTuple(const Features& tmpl, const Features& test, float minScore)
		: _data{
			FeaturesStatData(BC::FN, tmpl, _sd, minScore),
			FeaturesStatData(BC::FP, test, _sd, 0)
		}
	{
		_sd.Reserve(tmpl.ItemsCount());
	}

	// sets counting local stats data (for given chromosome)
	//	@param tmplIt: template chromosome's iterator
	//	@param testIt: test chromosome's iterator
	void SetChrom(Features::cIter tmplIt, Features::cIter testIt)
	{
		_cID = CID(tmplIt);
		_data[0].SetChrom(tmplIt);
		_data[1].SetChrom(testIt);
	}

	// stops counting local stats data
	void ResetChrom()
	{
		_data[0].ResetChrom();
		_data[1].ResetChrom();
		_sd.ResetChrom();
	}

	void Treat() { DiscardNonOverlapRegions<FeaturesStatData>(_data, 1); }

	// prints BC, F1 and SD
	//	@param t: target set for computation
	void PrintStat(eTarget t) const
	{
		if (t == CHR)	cout << Chrom::AbbrName(_cID) << COLON << TAB;
		else			PrintFooter();
		_data[0].PrintStat(t);
		_data[1].PrintStat(t);
		PrintF1(t);
		cout << _sd.GetSD(t) << LF;
	}
};
