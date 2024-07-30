/**********************************************************
Main.h for FGStest - Features Gold Standard test
2024 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 07/30/2024
-------------------------
***********************************************************/

#pragma once
#include "CrossRgns.h"
#include "Features.h"

enum optValue {		// options id
	//oGEN,
	oCHROM,
	oTEMPL,
	oMIN_DEV,
	oMIN_SCORE,
	oALARM,
	oDUMP_FILE,
	oDOUT_FILE,
	oTIME,
	oVERSION,
	oSUMM,
	oHELP,
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

	// Incorrect Features dump file
	class IncFWriter : public FormWriter
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
		IncFWriter(const char* fname) : FormWriter(fname)
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

	// 'Unified Data' holds tuple common data
	//	Data is filled by 'FeaturesStatData' instance and read by 'FeaturesStatTuple' instance
	class UniData
	{
		unique_ptr<IncFWriter> _oFile;		// dump file
		StandDev _sd;
		chrlen	 _abnormDevCnt[2]{ 0,0 };	// count of abnormal (too big) deviations for chrom, for all
		short	 _minDev = 0;				// minimum deviation for the deviation issue accounting

	public:
		// Returns count of abnormal deviations
		//	@param t: for the current chromosome or for all
		chrlen GetAbnormDevCnt(eTarget t) const { return _abnormDevCnt[t]; }

		// Returns Standard Deviation
		//	@param t: for the current chromosome or for all
		float GetSD(eTarget t) const { return _sd.GetSD(t); }

		// Initializes the instance
		//	@param capacity: standard deviation capacity
		//	@param minDev: minimum deviation for the deviation issue accounting
		//	@param fname: name of dump file or NULL
		void Init(size_t capacity, short minDev, const char* fname)
		{
			_sd.Reserve(capacity);
			_minDev = minDev;
			if (fname)	_oFile.reset(new IncFWriter(fname));
		}

		void SetChrom(chrid cID) { if (_oFile)	_oFile->SetChrom(cID); }

		void ResetChrom() {
			_sd.ResetChrom(); 
			_abnormDevCnt[ALL] += _abnormDevCnt[CHR];
			_abnormDevCnt[CHR] = 0;
		}

		// Saves False Poitive/Negative issue
		void Discard(Features::cItemsIter it, BC::eBC bc)
		{
			if (_oFile)	_oFile->WriteFF(it, bc);
		}

		// Saves abnormal deviation issue
		void AcceptDev(Features::cItemsIter it, short dev)
		{
			_sd.AddDev(dev);
			if (_minDev && dev > _minDev) {
				_abnormDevCnt[CHR]++;
				if (_oFile)	_oFile->WriteFF(it, dev);
			}
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
			dout << setw(4) << count[0] << delim << setprecision(3) << float(count[0]) / count[1] << delim;
		}

		const Features& _fs;
		UniData&	_uData;
		iterator	_beginII;
		iterator	_endII;
		BC::eBC		_bc;				// needed for printing to dump file
		float		_minScore;
		/*
		* valid feature's counters (per chrom and total) are only needed
		* for autonomous calculation of FNR and FDR (when printing).
		* In FeaturesStatTuple::PrintF1() we can use _sd.size() as True Positive\a
		*/
		chrlen	_cntBC[2][2]{
			{ 0,0 },	// false, total valid feature's count for chromosome
			{ 0,0 }		// false, total valid feature's count in total
		};

	public:
		FeaturesStatData(BC::eBC bc, const Features& fs, UniData& uData, float minScore)
			: _bc(bc), _fs(fs), _uData(uData), _minScore(minScore) {}

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
			_uData.AcceptDev(it[1], short(int(it[0]->Centre()) - it[1]->Centre()));
		}
		void Discard(iterator it)
		{
			if (it->Value < _minScore)
				_cntBC[CHR][TTL]--;
			else {
				_cntBC[CHR][FLS]++;
				_uData.Discard(it, _bc);
			}
		}
	};


	static const USHORT titleLineLen = 53;
	static void PrintFooter()
	{
		PrintSolidLine(titleLineLen);
		dout << sTotal << COLON << TAB;
	}

	FeaturesStatData _data[2];
	UniData			 _uData;

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
		const char* delim = "  ";
		dout << setw(4) << _uData.GetAbnormDevCnt(t) + *_data[0].BC(t) + *_data[1].BC(t) << delim;	// issues count
		_data[0].PrintStat(t);
		_data[1].PrintStat(t);
		dout << delim << GetF1(t) << delim << _uData.GetSD(t) << LF;
	}

public:
	static void PrintHeader()
	{
		dout << setw(15) << "issues   " << BC::Title(BC::FN) << "   FNR     " << BC::Title(BC::FP) << "   FDR      F1     SD\n";
		PrintSolidLine(titleLineLen);
	}

	FeaturesStatTuple(const Features& smpl, const Features& test, float minScore, short minDev, const char* fname)
		: _data{
			FeaturesStatData(BC::FN, smpl, _uData, minScore),
			FeaturesStatData(BC::FP, test, _uData, 0)
		}
	{ _uData.Init(smpl.ItemsCount(), minDev, fname); }

	// calculates and print chromosome's statistics
	//	@param sIt: sample chromosome's iterator
	//	@param tIt: test chromosome's iterator
	void GetChromStat(Features::cIter sIt, Features::cIter tIt)
	{
		// set chrom's data
		_data[0].SetChrom(sIt);
		_data[1].SetChrom(tIt);
		_uData.SetChrom(CID(sIt));
		// treat
		DiscardNonOverlapRegions<FeaturesStatData>(_data, 1);
		// print stats
		dout << Chrom::AbbrName(CID(sIt)) << COLON << TAB;
		PrintStat(CHR);
		// reset chrom's data
		_data[0].ResetChrom();
		_data[1].ResetChrom();
		_uData.ResetChrom();
	}

	// prints BC, F1 and SD
	void PrintTotalStat() const
	{
		PrintFooter();
		PrintStat(ALL);
	}
};
