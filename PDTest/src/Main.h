/**********************************************************
PDTest - Peak Detectors test
Copyright (C) 2024 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07/25/2024
-------------------------
***********************************************************/

#pragma once
#include "Features.h"

enum optValue {		// options id
	//oGEN,
	//oCHROM,
	oTEMPL,
	oMIN_DEV,
	oMIN_SCORE,
	oALARM,
	oOUTFILE,
	oTIME,
	oVERSION,
	oSUMM,
	oHELP,
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
		const fraglen ROI_ext = 500;
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

// 'Formatted Output File' - simple output text file with formatted line writer
class FOutFile
{
	FILE* _file;

public:
	FOutFile(const char* name) { _file = fopen(name, "w"); }

	~FOutFile() { fclose(_file); }

	template<typename... Args>
	void Write(const char* format, Args ... args) { fprintf(_file, format, args...); }
};


// A tuple of two Feature files - sample and test - that calculates and prints BC, F1 and SD statistics
class FeaturesStatTuple
{
	// target set for computation
	enum eTarget {
		CHR,	// for the current chromosome
		ALL		// for all
	};

	// false or total enumeration
	enum eFT {
		FLS,	// false value
		TTL		// total values
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

	// Keeps data and calculates Standard Deviation
	class StandDev
	{
		vector<short> _devs;		// deviations
		UINT	_startInd = 0;		// starting index in _devs for the current chromosome
		int		_sumDev[2]{ 0,0 };	// sum of deviation for current chromosome, for all

	public:
		void Reserve(size_t capacity) { _devs.reserve(capacity); }

		// Returns Standard Deviation
		//	@param t: for the current chromosome or for all
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

	// False Features dump file
	class FF_OutFile : public FOutFile
	{
		// returns normalized score
		static float NormScore(const Features::cItemsIter& it)
		{
			float score = it->Value;
			if (score > 1)	score /= 1000;
			return score;
		}

		IGVlocus _locus;

	public:
		FF_OutFile(const char* fname) : FOutFile(fname)
		{
			Write("chrom\tiss start\tend  \tscore\tdev\tlocus\n");
		}

		void SetChrom(chrid cID) { _locus.SetChrom(cID); }

		// Writes BC false feature line
		void WriteFF(const Features::cItemsIter& it, BC::eBC bc)
		{
			Write("%s\t%-4s%d\t%d\t%.2f\t \t%s\n",
				_locus.ChromAbbrName(),
				BC::Title(bc),
				it->Start,
				it->End,
				NormScore(it),
				_locus.Print(it->Start, it->End)
			);
		}

		// Writes deviation false feature line
		void WriteFF(const Features::cItemsIter& it, short dev)
		{
			Write("%s\tdev %d\t%d\t%.2f\t%d\t%s\n",
				_locus.ChromAbbrName(),
				it->Start,
				it->End,
				NormScore(it),
				dev,
				_locus.Print(it->Start, it->End)
			);
		}
	};

	// 'Features' wrapper for manipulating Features iterators for a given chromosome 
	class FeaturesStatData
	{
	public:
		using iterator = Features::cItemsIter;

	private:
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
		short		_minDev;
		float		_minScore;
		unique_ptr<FF_OutFile>& _oFile;
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
		FeaturesStatData(BC::eBC bc, const Features& fs, StandDev& sd, float minScore, short minDev, unique_ptr<FF_OutFile>& oFile)
			: _bc(bc), _fs(fs), _sd(sd), _minScore(minScore), _minDev(minDev), _oFile(oFile) {}

		// returns binary classifiers
		//	@param t: for the current chromosome or for all
		const chrlen* BC(eTarget t) const { return _cntBC[t]; }

		// prints binary classifiers
		//	@param t: for the current chromosome or for all
		void PrintStat(eTarget t) const { PrBcCounts(_cntBC[t]); }

		// sets counting local stats data (for given chromosome)
		//	@param cIt: chromosome's iterator
		void SetChrom(Features::cIter cIt)
		{
			auto& data = _fs.Data(cIt);
			_cntBC[CHR][TTL] = chrlen(data.ItemsCount());
			_beginII = _fs.ItemsBegin(data);
			_endII = _fs.ItemsEnd(data);
		}

		// stops counting local stats data
		void ResetChrom()
		{
			_cntBC[ALL][FLS] += _cntBC[CHR][FLS];
			_cntBC[ALL][TTL] += _cntBC[CHR][TTL];
			memset(_cntBC[CHR], 0, sizeof(_cntBC) / 2);
		}

		iterator& begin() { return _beginII; }
		iterator& end()	{ return _endII; }

		static chrlen	Start(iterator it) { return it->Start; }
		static chrlen	End(iterator it) { return it->End; }
		static bool		IsWeak(iterator it) { return false; }	// stub
		void Accept(const iterator it[2])
		{
			auto dev = short(int(it[0]->Centre()) - it[1]->Centre());
			_sd.AddDev(dev);
			if (_minDev && _oFile && dev > _minDev)		_oFile->WriteFF(it[1], dev);
		}
		void Discard(iterator it)
		{
			if (it->Value < _minScore)
				_cntBC[CHR][TTL]--;
			else {
				_cntBC[CHR][FLS]++;
				if (_oFile)		_oFile->WriteFF(it, _bc);
			}
		}
	};


	static const USHORT titleLineLen = 44;
	static void PrintFooter()
	{
		PrintSolidLine(titleLineLen);
		cout << sTotal << COLON << TAB;
	}

	StandDev _sd;
	FeaturesStatData _data[2];
	unique_ptr<FF_OutFile> _oFile;

	// returns F1 score
	//	@param t: for the current chromosome or for all
	float GetF1(eTarget t) const
	{
		auto FN = _data[0].BC(t)[FLS];				// False Negative
		auto dTP = 2 * (_data[0].BC(t)[TTL] - FN);	// double True Positive
		return float(dTP) / (dTP + _data[1].BC(t)[FLS] + FN);
	}

	// prints BC, F1 and SD
	//	@param t: for the current chromosome or for all
	void PrintStat(eTarget t) const
	{
		_data[0].PrintStat(t);
		_data[1].PrintStat(t);
		const char* delim = "  ";
		cout << delim << GetF1(t)
			 << delim << _sd.GetSD(t) << LF;
	}

public:
	static void PrintHeader()
	{
		cout << setw(11) << BC::Title(BC::FN) << "   FNR   " << BC::Title(BC::FP) << "   FDR      F1     SD\n";
		PrintSolidLine(titleLineLen);
	}

	FeaturesStatTuple(const Features& smpl, const Features& test, float minScore, short minDev, const char* fname)
		: _data{
			FeaturesStatData(BC::FN, smpl, _sd, minScore, minDev, _oFile),
			FeaturesStatData(BC::FP, test, _sd, 0, minDev, _oFile)
		}
	{
		_sd.Reserve(smpl.ItemsCount());
		if (fname)	_oFile.reset(new FF_OutFile(fname));
	}

	// calculates and print chromosome's statistics
	//	@param sIt: sample chromosome's iterator
	//	@param tIt: test chromosome's iterator
	void GetChromStat(Features::cIter sIt, Features::cIter tIt)
	{
		// set chrom's data
		_data[0].SetChrom(sIt);
		_data[1].SetChrom(tIt);
		if (_oFile)	_oFile->SetChrom(CID(sIt));
		// treat
		DiscardNonOverlapRegions<FeaturesStatData>(_data, 1);
		// print stats
		cout << Chrom::AbbrName(CID(sIt)) << COLON << TAB;
		PrintStat(CHR);
		// reset chrom's data
		_data[0].ResetChrom();
		_data[1].ResetChrom();
		_sd.ResetChrom();
	}

	// prints BC, F1 and SD
	void PrintTotalStat() const
	{
		PrintFooter();
		PrintStat(ALL);
	}

};
