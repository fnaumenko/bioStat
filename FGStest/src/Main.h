/**********************************************************
Main.h for FGStest - Features Gold Standard test
2024 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 08/09/2024
-------------------------
***********************************************************/

#pragma once
#include "CrossRgns.h"
#include "Features.h"

enum optValue {		// options id
	//oGEN,
	oCHROM,
	oTEMPL,
	oMIN_CDEV,
	oMIN_LDEV,
	oMIN_SCORE,
	oEXPAND,
	oALARM,
	oISSUE_FILE,
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

	// Incorrect Features issues file
	class IssBedWriter : public FormWriter
	{
		// returns normalized score
		static float NormScore(const Features::cItemsIter& it)
		{
			float score = it->Value;
			if (score > 1)	score /= 1000;
			return score;
		}

		IGVlocus _locus;

		template<typename T>
		void WriteFF(const char* format, const char* title, const Features::cItemsIter& it, T val)
		{
			Write(format,
				_locus.ChromAbbrName(),
				it->Start,
				it->End,
				title,
				NormScore(it),
				val,
				_locus.Print(it->Start, it->End)
			);
		}

	public:
		IssBedWriter(const char* fname) : FormWriter(fname)
		{
			Write("#chrom\t  start\t    end\tiss\tscore\tdev\tlocus\n");
		}

		void SetChrom(chrid cID) { _locus.SetChrom(cID); }

		// Writes BC false feature line
		void WriteFF(const Features::cItemsIter& it, BC::eBC bc)
		{
			WriteFF<char>("%s\t%d\t%d\t%s\t%.2f\t%c\t%s\n", BC::Title(bc), it, SPACE);
		}

		// Writes centre deviation false feature line
		void WriteFF(const Features::cItemsIter& it, short dev)
		{
			WriteFF<short>("%s\t%d\t%d\t%s\t%.2f\t%d\t%s\n", "cD", it, dev);
		}

		// Writes width deviation false feature line
		void WriteFF(const Features::cItemsIter& it, float dev)
		{
			WriteFF<float>("%s\t%d\t%d\t%s\t%.2f\t%.1f\t%s\n", "lD", it, dev);
		}
	};

	// Keeps data and calculates Standard Deviation
	template<typename T>
	class StandDev
	{
		vector<T>	_devs;				// deviations
		IssBedWriter*	_oFile = nullptr;
		T		_minDev = 0;			// minimum deviation for the deviation issue accounting
		UINT	_startInd = 0;			// starting index in _devs for the current chromosome
		float	_sumDev[2]{ 0.,0. };	// sum of deviation for current chromosome, for all
		chrlen	_abnormDevCnt[2]{ 0,0 };// count of abnormal (too big) deviations for chrom, for all

	public:
		void Init(size_t capacity, T minDev, IssBedWriter* oFile)
		{
			_minDev = minDev;
			_oFile = oFile;
			_devs.reserve(capacity);
		}

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
		//	@param dev: deviation value
		//	@param acceptAbnorm: if true then take into account the abnormal deviation value
		//	@param it: Features iterator to save record into file
		//	@returns: true is abnormal deviation value was taken into account
		bool AddDev(T dev, bool acceptAbnorm, const Features::cItemsIter it)
		{
			_devs.push_back(dev);
			_sumDev[CHR] += dev;
			if (acceptAbnorm && _minDev && dev > _minDev) {
				_abnormDevCnt[CHR]++;
				if (_oFile)	_oFile->WriteFF(it, dev);
				return true;
			}
			return false;
		}

		chrlen GetAbnormDevCount(eTarget t) const { return _abnormDevCnt[t]; }

		// Stops counting data for the current chromosome
		void ResetChrom()
		{
			_startInd = UINT(_devs.size() - 1);
			_sumDev[ALL] += _sumDev[CHR];
			_sumDev[CHR] = 0;
			_abnormDevCnt[ALL] += _abnormDevCnt[CHR];
			_abnormDevCnt[CHR] = 0;
		}
	};

	// 'Unified Data' holds tuple common data
	//	Data is filled by 'FeaturesStatData' instance and read by 'FeaturesStatTuple' instance
	class UniData
	{
		unique_ptr<IssBedWriter> _oFile;		// issues file
		StandDev<short> _cSD;	// centre standars deviation
		StandDev<float> _lSD;	// length standars deviation

	public:
		// Returns count of centre abnormal deviations
		//	@param t: for the current chromosome or for all
		chrlen GetCentreAbnormDevCount(eTarget t) const { return _cSD.GetAbnormDevCount(t); }

		// Returns count of length abnormal deviations
		//	@param t: for the current chromosome or for all
		chrlen GetWidthAbnormDevCount(eTarget t) const { return _lSD.GetAbnormDevCount(t); }

		// Returns centre Standard Deviation
		//	@param t: for the current chromosome or for all
		float GetCentreSD(eTarget t) const { return _cSD.GetSD(t); }

		// Returns length Standard Deviation
		//	@param t: for the current chromosome or for all
		float GetLengthSD(eTarget t) const { return _lSD.GetSD(t); }

		// Initializes the instance
		//	@param capacity: standard deviation capacity
		//	@param minCDev: minimum centre deviation for the deviation issue accounting
		//	@param minLDev: minimum length deviation for the deviation issue accounting
		//	@param fname: name of issues file or NULL
		void Init(size_t capacity, short minCDev, float minLDev, const char* fname)
		{
			if (fname)	_oFile.reset(new IssBedWriter(fname));
			_cSD.Init(capacity, minCDev, _oFile.get());
			_lSD.Init(capacity, minLDev, _oFile.get());
		}

		void SetChrom(chrid cID) { if (_oFile)	_oFile->SetChrom(cID); }

		void ResetChrom() { _cSD.ResetChrom(); _lSD.ResetChrom(); }

		// Saves False Poitive/Negative issue
		void Discard(const Features::cItemsIter it, BC::eBC bc)
		{
			if (_oFile)	_oFile->WriteFF(it, bc);
		}

		// Saves deviation issue
		void AcceptDev(const Features::cItemsIter it[2])
		{
			bool acceptAbnorm = _cSD.AddDev(
				short(int(it[0]->Centre()) - it[1]->Centre()),
				true, it[1] 
			);
			_lSD.AddDev(
				float(it[1]->Length()) / it[0]->Length(),
				!acceptAbnorm, it[1]
			);
		}
	};

	// 'Features' wrapper for manipulating Features iterators for a given chromosome 
	class FeaturesStatData
	{
	public:
		using iterator = Features::cItemsIter;

	private:
		const Features& _fs;
		UniData&	_uData;
		iterator	_beginII;
		iterator	_endII;
		BC::eBC		_bc;				// needed for printing to issues file
		USHORT		_expandVal;
		float		_minScore;
		/*
		* valid feature's counters (per chrom and total) are only needed
		* for autonomous calculation of FNR and FDR (when printing).
		* In FeaturesStatTuple::PrintF1() we can use _cSD.size() as True Positive\a
		*/
		chrlen	_cntBC[2][2]{
			{ 0,0 },	// false, total valid feature's count for chromosome
			{ 0,0 }		// false, total valid feature's count in total
		};

	public:
		FeaturesStatData(BC::eBC bc, const Features& fs, UniData& uData, USHORT	expandVal, float minScore)
			: _bc(bc), _fs(fs), _uData(uData), _expandVal(expandVal), _minScore(minScore) {}

		// Returns binary classifiers
		//	@param t: for the current chromosome or for all
		const chrlen* BC(eTarget t) const { return _cntBC[t]; }

		// Returns false rate
		//	@param t: for the current chromosome or for all
		const float Rate(eTarget t) const { return float(_cntBC[t][FLS]) / _cntBC[t][TTL]; }

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

		chrlen	Start(const iterator it)	 { return it->Start - _expandVal; }
		chrlen	End	(const iterator it)		 { return it->End + _expandVal; }
		static bool	IsWeak(const iterator it){ return false; }	// stub
		void Accept	(const iterator it[2])	 { _uData.AcceptDev(it); }
		void Discard(const iterator it)
		{
			if (it->Value < _minScore)
				_cntBC[CHR][TTL]--;
			else {
				_cntBC[CHR][FLS]++;
				_uData.Discard(it, _bc);
			}
		}
	};


	static const USHORT titleLineLen = 73;
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
		chrlen cDev	= _uData.GetCentreAbnormDevCount(t);
		chrlen lDev = _uData.GetWidthAbnormDevCount(t);
		chrlen sampleBC	= *_data[0].BC(t);
		chrlen testBC	= *_data[1].BC(t);

		dout << setprecision(3) << _data[0].Rate(t) << TAB << _data[1].Rate(t)
			<< TAB << GetF1(t) << TAB << _uData.GetCentreSD(t) << TAB << _uData.GetLengthSD(t);
		dout << TAB << '|' << setw(4) << cDev + lDev + sampleBC + testBC	// issues count
			<< setw(5) << sampleBC << setw(5) << testBC << setw(5) << cDev << setw(5) << lDev << LF;
	}

public:
	static void PrintHeader()
	{
		dout << setw(6) << SPACE << TAB << "FNR" << SPACE << TAB << "FDR" << SPACE
			<< TAB << "F1" << "  \t" << "c-SD" << TAB << "l-SD";
		dout << "\t|total" << setw(4) << BC::Title(BC::FN) << setw(5) << BC::Title(BC::FP)
			<< setw(5) << "c-D" << setw(5) << "l-D" << LF;
		PrintSolidLine(titleLineLen);
	}

	FeaturesStatTuple(
		const Features& smpl,
		const Features& test,
		USHORT expandVal,
		float minScore,
		short minCDev,
		float minLDev,
		const char* fname
	)
		: _data{
			FeaturesStatData(BC::FN, smpl, _uData, expandVal, minScore),
			FeaturesStatData(BC::FP, test, _uData, 0, 0)
		}
	{ _uData.Init(smpl.ItemsCount(), minCDev, minLDev, fname); }

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
