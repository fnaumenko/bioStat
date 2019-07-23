/**********************************************************
Calc.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 17.04.2019
-------------------------
Provides classes for calculating CC
***********************************************************/
#pragma once

#include "Data.h"
#include "bioCC.h"

typedef pair<double, double> pairDbl;
typedef pair<Regions::Iter, Regions::Iter> RegionsRange;

enum eRS {	// defines types of printed template regions (features) CC
	rsOFF = 0,		// not set: it is never pointed in command line
	rsR	= 1,		// sorted by regions
	rsC	= 2			// sorted by coefficients
};

// 'CCkey' Represents correlation coefficient constants and methods
static struct CCkey
{
private:
	static const int Undef = 666;	// undefined coefficient

public:
	enum eCC {	// defines types of CC
		ccP = 0x1,	// Pearson correlation coefficient: first bit
		ccS = 0x2,	// signal correlation coefficient: second bit
	};

	// Returns true if parameter is signal coefficient
	static inline bool IsS(eCC ecc)		{ return (ecc&ccS) != 0; }

	// Returns true if parameter is Pearson coefficient
	static inline bool IsP(eCC ecc)		{ return static_cast<bool>(ecc&ccP); }

	// Returns true if parameter is both coefficients
	static inline bool IsBoth(eCC ecc)	{ return static_cast<bool>(ecc&ccS && ecc&ccP); }

	// Returns correlation coefficient
	static inline double CalcR(double v1, double v2, double cov)
	{ return  !v1 && !v2 ? Undef : cov/(sqrt(v1) * sqrt(v2)); }	// not sqrt(v1*v2) 

	static inline void Print(double val) {
		if( val==Undef || isnan(val) )	dout << "UNDEFINED";
		else	dout << val;
	}
} cc;

// 'CC' represents a pair of correlation coeficients: Pearson(first) and signal(second),
// and provides theirs output.
struct CC : pairDbl
{
private:
	static const int Empty = -2;	// uninitialised coefficient

public:
	inline CC() { first = second = Empty; }

	inline CC(double p, double s) { first = p; second = s; }

	inline bool operator < (const CC & ccres) const { 
		return first != Empty ? first < ccres.first : second < ccres.second; 
	}

	// Sets Pearson coefficient
	inline void SetP(double val)	{ first = val;	}

	// Sets signal coefficient; undefined value by default
	inline void SetS(double val=0)	{ second = val;	}

	// Returns true if even one value in pair is setting
	inline bool NotEmpty() const { return( first != Empty || second != Empty ); }

	// returns single value
	inline double GetSingleVal() const { return first != Empty ? first : second; }

	// set single negative value to absolute value
	void SetSingleAbsVal() { 
		if( first != Empty )	first = fabs(first);
		else					second = fabs(second);
	}

	void Print() const { 
		if( first != Empty )	{
			CCkey::Print(first);
			dout << TAB;
			//if( SCNT(first) < 8 )	dout << TAB;
			if( DigitsCount(first) < 8 )	dout << TAB;
		}
		if( second != Empty )
			CCkey::Print(second);
		dout << EOL; 
	}
};

//  CCaggregate represents accumulative values needed for calculation CC for a set of arrays
struct CCaggr
{
private:
	static const UINT Undef = -1;	// uninitialised coefficient
	ULLONG check;					// needs for check-up only
	ULLONG VarS[3];								// Signal variances
	double VarP1, VarP2, CovP, Mean1, Mean2;	// Pearson variances

public:
	// Increases S value with check-up for exceeding
	//	@i: index of value on VarS
	//	@val: value to increase
	//	@cID: chrom's ID, needed for exception only
	//	@return: false if all right 
	bool IncreaseSVar(BYTE i, ULLONG val, chrid cID);

public:
	CCaggr() {
		VarP1 = VarP2 = CovP = Mean1 = Mean2 = 0;
		//VarS[0] = VarS[1] = VarS[2] = 0;
		memset(VarS, 0, 3*sizeof(ULLONG));
	}
	// Sets predefined means
	inline void SetMeans(const pairDbl& means) {
		Mean1 = means.first;
		Mean2 = means.second;
	}

	// Increases P values in consideration of means
	//	@x: first value to increase
	//	@y: second value to increase
	void IncreasePVars(double x, double y);

	// Increases S values with check-up for exceeding
	//	@x: first value to increase
	//	@y: second value to increase
	//	@cID: chrom's ID, needed for exception only
	//	@return: false if all right 
	bool IncreaseSVars(ULLONG x, ULLONG y, chrid cID);

	// True if Signal coefficient is undefined (exceeding digital limit)
	inline bool IsUndefS() { return static_cast<bool>(VarS[0] == Undef); }

	// Calculates anr returns coefficients
	CC GetR() const;
};

// Wrapper for ChromSizes<RegionsRange> object to implement the interface for the common functionality
// of Bedf and DefRegions.
// Allows pass through these objects via iterator in a common way.
class ShellGenRegions : public Chroms<RegionsRange>
{
public:
	// Returns an iterator pointing to the first Region of pointed chrom
	//	@it: chromosome's iterator
	const Regions::Iter ChromBegin(Chroms<RegionsRange>::cIter it) const { return Data(it).first; }

	// Returns an iterator referring to the past-the-end Region of pointed chrom
	//	@it: chromosome's iterator
	const Regions::Iter ChromEnd(Chroms<RegionsRange>::cIter it) const { return Data(it).second; }

	// Constructor by DefRegions
	ShellGenRegions(const DefRegions& gRgns) {
		for(DefRegions::cIter it=gRgns.cBegin(); it!=gRgns.cEnd(); it++)
			AddVal(CID(it), RegionsRange(gRgns.Data(it).Begin(), gRgns.Data(it).End()));
	}

	// Constructor by Features
	ShellGenRegions(const Features& fs) {
		for(Features::cIter it=fs.cBegin(); it!=fs.cEnd(); it++)
			AddVal(CID(it), RegionsRange(fs.FeaturesBegin(it), fs.FeaturesEnd(it)));
	}

	// Gets number of regions of pointed chrom
	//	@it: chromosome's iterator
	inline chrlen RegionsCount(Chroms<RegionsRange>::cIter it) const {
		return Data(it).second - Data(it).first; }

	// Gets summary length of regions of pointed chrom
	//	@it: chromosome's iterator
	chrlen RegionsLength(Chroms<RegionsRange>::cIter it) const { 
		// need singleton?
		chrlen len = 0;
		for(Regions::Iter iter=Data(it).first; iter<Data(it).second; iter++)
			len += iter->Length();
		return len;
	}
};

// 'Results' accumulates calculated coeficients for each chromosome and total,
// and provides theirs output
class Results : public Chroms<CC>
{
public:
	enum eCCPrint {
		cIND = 0x1,	// output local CC
		cTTL = 0x2,	// output total CC
	};

private:
	CC _total;
	eCCPrint _printCC;

public:
	// creates inctance with defined local/total inquires
	inline Results(UINT ccPrint) :
		_printCC(eCCPrint(ccPrint)) {}

	// Gets inquire about local CC: true if should be calculated
	inline bool GiveLocal() const { return (_printCC & cIND) != 0; }
	
	// Gets inquire about total CC: true if should be calculated
	inline bool GiveTotal() const { return (_printCC & cTTL) != 0; }

	// Sets total results
	inline void	SetTotal(CC total) { _total = total; }

	// Prints results for each chromosome and total or total only (if it was defined).
	//	@printTitles: if true, print titles, otherwise does not only in case of single result
	void Print(bool printTitles);
};

typedef Array<chrlen>	chrlens;

// 'ChromLengths' represents chrlen array with extended methods for calculating CC
class ChromLengths : public chrlens
{
private:
	mutable double _mean;

	// Gets a pair of means of subarrays for this instance and given array
	pairDbl GetMeans(bool notGotThis, const chrlens& arr, long begin, long end) const;

public:
	inline ChromLengths(chrlen cLen = 0) : _mean(0), chrlens(cLen) {}

	// Gets a pair of means of subarrays for this instance and given array,
	// and set _mean for each array in case of whole arrays
	pairDbl SetMeans(const ChromLengths& arr, long begin=0, long end=0) const;

	// Returns maximal value in region
	chrlen GetMaxVal(long begin=0, long end=0) const;

	// Multiplys value in region to ratio
	void MultiplyVal(float ratio, long begin=0, long end=0);

	// Multiplys all values in region to ratio synchronously for pair of arrays
	//static void SynchMultiplyVal(Array<T>& arr1, Array<T>& arr2, 
	//	float ratio1, float ratio2, long begin=0, long end=0);

	// Calculates correlation coefficients
	//	@cID: chrom's ID. Needed for exception only
	//	@ecc: identifier what coefficient to calculate
	//	@arr: second array to compare
	//	@begin: low boundary of calculated range
	//	@end: high boundary of calculated range
	//	@return: pair of coefficients
	CC const GetR (chrid cID, CCkey::eCC ecc, const ChromLengths& arr, long begin=0, long end=0);

	// Accumulates variances and covariances for calculating CC
	//	@cID: chrom's ID. Needed for exception only
	//	@ecc: identifier what coefficient to calculate
	//	@ccaggr: accumulative aggregate
	//	@arr: second array to compare
	//	@begin: low boundary of calculated range
	//	@end: high boundary of calculated range
	void AccumVars (chrid cID, CCkey::eCC ecc, CCaggr & ccaggr, const chrlens& arr,
		long begin=0, long end=0) const;
};

// 'ChromsMap' is a container of Chroms's arrays (one array for one chromosome).
// Provides methods to calculate R.
class ChromsMap : public Obj, public Chroms<ChromLengths>
{
protected:
	// keeps data pointers and temporary variables needed for reading constructor only
	class Pocket
	{
	protected:
		const ChromSizes& _cSizes;	// chrom sizes only to control chrom length exceeding
		ChromLengths*  _map;	// current chroms map (to avoid call indexator each time)
		chrlen		_cSize;		// size of current adding chromosome

	public:
		inline Pocket(const ChromSizes& cs) : _cSizes(cs) {}

		// sets current chrom's size and reserves chrom's map
		inline void Reserve(chrid cID, ChromLengths* map) {
			(_map = map)->Reserve((_cSize = _cSizes[cID])/_Space);
		}
	};

	static BYTE	_Space;

	double	_binWidth;	// width of bins of histogram; if 0, no histogram
	eRS		_printFRes;	// sign to print results for each feature from 'template' and how to sort it

private:
	// Gets a pair of total genome means between this and second ChromsMap
	//	@gRgn: genome defined regions
	//	@map: second ChromsMap
	//	return: a pair of means for each set of chroms
	pairDbl GetGenomeMean(const DefRegions& gRgn, const ChromsMap& map) const;

public:
	inline ChromsMap() : 
		_binWidth (Options::GetFVal(oBINWIDTH)),
		_printFRes(eRS(Options::GetIVal(oFRES))) {}

	inline const ChromLengths & operator[] (chrid cID) const { return At(cID).Data; }

	// Calculates r for each region and fills results
	//	@cc: type of correlation coefficient
	//	@wig: ChromsMap object to correlate with
	//	@shGRgns: shell of treated chrom's regions
	//	@results: object to fill results
	void CalcRegionsR(CCkey::eCC ecc, const ChromsMap& wig,
		const ShellGenRegions& shGRgns, Results& results);

	// Calculates r and fills results
	//	@cc: type of correlation coefficient
	//	@wig: ChromsMap object to correlate with
	//	@gRgns: real chrom's regions
	//	@templ: template with defined regions or NULL
	//	@results: object to fill results
	void CalcR(CCkey::eCC ecc, const ChromsMap& wig, const DefRegions& gRgns, 
		const Features* templ, Results & results);

#ifdef DEBUG
	// Prints regions.
	//	@rgnCnt: number of regions to print or all if not specified
	void Print(chrlen rgnCnt=0) const;
#endif
};

// 'WigMap' encapsulates wig variableStep format. It's a wig-specialized wrapper of ChromsMap.
class WigMap : public ChromsMap
{
private:
	static const BYTE	_CntFields = 2;	// number of field in tab file for data line

	class WigPocket : public ChromsMap::Pocket
	{
		p_ulong _recCnt;	// number of: first - registered, second - written WIG records

	public:
		inline WigPocket(const ChromSizes& cs) : _recCnt(0, 0), Pocket(cs) {}

		// Gets number of registered and written WIG records
		inline p_ulong RecCount() const { return _recCnt; }

		// Adds region for current chromosome and increases record counters.
		//	@start:	region's start position
		//	@size:	region's length
		//	@val:	region's value
		//	@spotter: spotter to collect statistics
		//	return: 1 if region is added; otherwise 0
		void AddRegion(chrlen start, chrlen size, chrlen val, Spotter& spotter);
	};

	// Gets an item's title
	//	@pl: true if plural form
	inline const string& ItemTitle(bool pl = false) const { return FT::ItemTitle(FT::WIG, pl); };

	// Initializes instance from wig file
	//	@spotter: spotter to control ambiguities
	//	@cSizes: chrom sizes to control chrom length exceedeing
	//	return: numbers of all and initialied items for given chrom
	p_ulong InitDerived	(Spotter& spotter, const ChromSizes& cSizes);

	// Adds chromosome and set it as current
	//	@cID: adding chromosome's ID
	//	@pocket: initializing temporary variables
	inline void AddChrom (chrid cID, Pocket& pocket) {
		pocket.Reserve(cID, &(AddEmptyElem(cID).Data));
	}

public:
	// Creates wigMap object.
	//	@fName: file name
	//	@cSizes: chrom sizes to control the chrom length exceedeng; if NULL, no control
	//	@info: type of feature ambiguties that should be printed
	//	@printName: true if file name should be printed unconditionally, otherwise deneds on info
	//	@primary: if true object is primary
	//	@gRgns: genome's regions
	inline WigMap(const char* fName, eInfo info, bool printName, bool primary, DefRegions& gRgns)
	{
		Spotter spotter(info, false, FT::WIG);
		Init(NULL, fName, spotter, gRgns.ChrSizes(), info > iLAC || printName, primary);
	}
};

// 'DensMap' encapsulates density map. It's a density-specialized wrapper of ChromsMap.
class DensMap : public ChromsMap
{
private:
	class DensPocket : public ChromsMap::Pocket
	{
	private:
		const Reads&  _reads;
		chrlen	_wCurrLen,		// length of uncsanning rest of current window
				_wi;			// current window's index
		Reads::cItemsIter	_currIt,
							_endIt;

		// Fills number or reads from window started with @start position
		void ScanWindow(chrlen start);

	public:
		inline DensPocket(const Reads& rs, const ChromSizes& cs) : 
			_reads(rs), Pocket(cs) {}

		// Fills number or reads from region @rgn
		//void ScanRegion(const Region&);
		void ScanChrom();

		// Adds chrom by Reads iter
		//	@it: Reads iter
		//	@map: pointer to the chroms map
		void AddChrom (Reads::cIter it, ChromLengths* map);
	};

	// Adds chromosome and set it as current
	//	@it: chromosome's bed iterator
	//	@pocket: initializing temporary variables
	inline void AddChrom (Reads::cIter it, DensPocket& pocket) {
		pocket.AddChrom(it, &(AddEmptyElem(CID(it)).Data));
	}

	// Gets an item's title
	//	@plural: true if plural form
	inline const string& ItemTitle(bool plural = false) const
	{ return strEmpty; }	// never called: to close Obj virtual method

	// Initializes empty instance to close Obj virtual method
	p_ulong InitDerived(Spotter&, const ChromSizes&)
	{ return make_pair(0, 0); }

public:
	DensMap(const Reads& reads, const ChromSizes& cs);

};

// 'ChromRanges' represented a set of chromosome's ranges.
struct ChromRanges : public ItemIndexes
/*
	To keep chromosome's number and first/last ranges indexes Features::ChromItemsInd is used;
	FirstInd/LastInd are keeping first/last ranges indexes,
	and lengths of all features from fs1, fs2 are keeping additionally for efficiency.
*/
{
	chrlen FeatrsLen1;		// length of all features of chromosome in fs1
	chrlen FeatrsLen2;		// length of all features of chromosome in fs2

	inline ChromRanges() : ItemIndexes(), FeatrsLen1(0), FeatrsLen2(0) {}

	inline ChromRanges(chrlen firstInd, chrlen lastInd, chrlen len1, chrlen len2) : 
		ItemIndexes(firstInd, lastInd), FeatrsLen1(len1), FeatrsLen2(len2) {}
};

// 'JointedBeds' represents two bed-files as a chromosomes collection and theirs joint features (ranges).
class JointedBeds : Chroms<ChromRanges>
/*
 * This is fast but complicated implementation of calculating algorithm.
 */
{
private:
	// 'PairR' pepresetns a pair of R classes accumulated variances & covariance and calculated r.
	class PairR
	/*
	 * First in pare is signal R, second - Pearson R.
	 */
	{
	private:
		// Class 'R' accumulates variances & covariance and calculated r
		class R {
		private: 
			double	_var1, _var2,	// variances 
					_cov;			// covariance: SUM( X-Xmean)^2 * (Y-Ymean)^2 )
	// Mean value for every range may accept only one of two values: 0-mean or 1-mean,
	// and they are static within chromosome.
	// So for efficiency we can keep theirs in 3 arrays:
	// two square's arrays: arrays of squared values from each bed;
	// there are only 2 combinations, but they are duplicated for the simplicity acces,
	// and one crossing array: array of all combinations
	// of multiplications different values from fs1 and fs2
			double	_sqMeans1[4], _sqMeans2[4],	// square's arrays for bed, fs2
					_crossMeans[4];				// crossing array
		public:
			inline R () : _var1(0), _var2(0), _cov(0) {}
			// Keeps mean values for fs1, fs2. 
			//	@clear: if true, clear instance for treatment of new chromosome
			void Init(double mean1, double mean2, bool clear);
			// Accumulates next length of range
			void Increment(chrlen len, char val);
			// Returns coefficient of correlation
			inline double Get() const {	return CCkey::CalcR(_var1, _var2, _cov); }
		};

		R	_s;	// signal R
		R	_p;	// Pearson R
		CCkey::eCC	_ecc;

	public:
		inline PairR (CCkey::eCC ecc) : _ecc(ecc) {}
		// Keeps mean values for fs1, fs2. 
		//	@clear: if true, clear instance for treatment of new chromosome
		void Init(double mean1, double mean2, bool clear) {
			if( CCkey::IsS(_ecc) )	_s.Init(0, 0, clear);
			if( CCkey::IsP(_ecc) )	_p.Init(mean1, mean2, clear);
		}
		// Accumulates next length of range
		void Increment(chrlen len, char val) {
			if( CCkey::IsS(_ecc) )	_s.Increment(len, val);
			if( CCkey::IsP(_ecc) )	_p.Increment(len, val);
		}
		// Returns coefficient of correlation
		CC Get() const {
			CC res;
			if( CCkey::IsS(_ecc) )	res.SetS(_s.Get());
			if( CCkey::IsP(_ecc) )	res.SetP(_p.Get());
			return res;
		}
	};

	// 'Range' represens a range.
	struct Range 
	// Range is a joint feature with start position and joint (combined) value.
	// All features from fs1, fs2 are splitted by contiguous ranges with predefined value:
	// VAL1 (only the first Features has a feature here), or
	// VAL2 (only the second Features has a feature here), or
	// VAL1 & VAL2 (both of the Beds have a feature here), or
	// 0 (no features for both of the Beds)
	{
		chrlen	Start;	// start position of range in standard chromosomal coordinates
		char	Val;	// value of range

		inline Range() : Start(0), Val(0) {}
		inline Range(chrlen pos, char val) : Start(pos), Val(val) {}
	};

	vector<Range>	_ranges;

public:
	// Fills ChromRanges & Range by given two beds.
	// Beds chromosomes should be checked as Treated.
	// Both fs1 & fs2 must be valid: no duplicated, crossed, adjacent, coverage features;
	// in other case R may be wrong
	JointedBeds(Features & fs1, Features & fs2);
	
	// Calculates r and fills results
	//	@cc: type of correlation coefficient
	//	@cSizes: chrom sizes
	//	@results: object to fill results
	void CalcR(CCkey::eCC ecc, const ChromSizes& cSizes, Results& results);
#ifdef DEBUG
	void	Print();
#endif
};

// 'CorrPair' represents pair of objects to compare, and methods for recognizing types and calculation CC
class CorrPair
{
private:
	struct FileType	// keeps pointers to the common methods
	{
		typedef Obj*	(CorrPair::*CreateObj)	(const char*, bool);
		typedef void	(*DeleteObj)			(Obj*);		// static function
		typedef bool	(CorrPair::*FillGenRegns)(DefRegions&);

		CreateObj		Create;			// constructor of type
		DeleteObj		Delete;			// destructor of type
		FillGenRegns	FillGenRgns;	// filler DefRegions by common chroms
	};
	
	static int _FileTypesCnt;		// count of treatment file's types
	static FileType _FileTypes[];

	Obj*	_firstObj;
	Obj*	_secondObj;
	DefRegions&	_gRgns;		// initial genome regions to correlate
	Features*	_templ;
	CCkey::eCC	_ecc;
	BYTE	_typeInd;		// type of corr. files: 0 - wig, 1 - bedF, 2 - bedR
	bool	_printWarn;		// true if print line warnings
	bool	_printName;		// print file names unconditionally
	BaseItems::eInfo	_info;

	// Returns true if info output is not laconic
	inline bool IsNotLac()	const { return _info > BaseItems::iLAC; }
	// Returns true if Reads are treating
	inline bool IsWig()		const { return !_typeInd; }
	// Returns true if Features are treating
	inline bool IsBedF()	const {	return _typeInd == 1; }
	// Returns true if BedE are treating
	inline bool IsAlign()	const { return _typeInd == 2; }

	// Creates features bed object.
	Obj* CreateBedF	(const char* fName, bool isPrimary);
	
	// Creates alignment object
	Obj* CreateBedR	(const char* fName, bool isPrimary);

	// Creates wigMap object.
	Obj* CreateWig	(const char* fName, bool isPrimary);

	// Checks first & second bedF for common chroms.
	//	@gRgns: doesn't used, states for common syntax only, to call via function pointer
	//	return: true if there are common chroms
	bool FillComnonChromsBedF(DefRegions& gRgns);

	// Checks first & second bedF for common chroms
	// and fill DefRegions by coverages or density common chroms
	//	@gRgns: DefRegions to fill
	//	return: true if there are common chroms
	bool FillComnonChromsMap(DefRegions& gRgns);

	// Calculates r for genome features.
	//	@firstBed: first Features to correlate
	//	@results: results holder
	inline void CalcCCBedF(Features& firstBed, Results& results) {
		// DefRegions with common chroms in it is not needed since common chroms 
		// are set automatically by JointedBeds
		JointedBeds(firstBed, *((Features*)_secondObj)).CalcR(_ecc, _gRgns.ChrSizes(), results);

		// DEBUG control check by BedMap
		//BedMap bMap(*((Features*)_secondObj), gRgns);
		//BedMap(firstBed, gRgns).CalcR(_ecc, bMap, results);
	}

	// Calculates r for genome coverages or density.
	//	@gRgns: DefRegions with common chroms in it
	//	@results: results holder
	inline void CalcCCMap(DefRegions& gRgns, Results& results) {
		((ChromsMap*)_firstObj)->CalcR(_ecc, *((ChromsMap*)_secondObj), gRgns, _templ, results);
	}

	// Checks file extisting and extention validity
	//	@fName: file's name
	//	@abortInvalid: if true throw extention if checking is false
	//	return: index in _FileTypes or _FileTypesCnt if checking is false
	BYTE CheckFileExt(const char * fName, bool abortInvalid);

public:
	// Creates an instance with checking primary object.
	//	@primefName: primary file's name
	//	@gRgns: genome regions
	//	@templ: name of template bed file, or NULL if undefined
	//	@multiFiles: true if more then one secondary files are placed
	CorrPair(const char* primefName, DefRegions& gRgns, const char* templ, bool multiFiles);

	~CorrPair();
	
	// Adds secondary object, calculates and prints CCkey.
	void CalcCC(const char * fName);
};

#ifdef DEBUG

class BedMap : public Chroms<ChromLengths>
/*
 * Class 'BedMap' represents bed-file as a list of byte's arrays (one array for one chromosome),
 * where each byte in array represents one nucleotide.
 * Bytes corresponding to the nucleotides included in features have value 1, not included - 0.
 * This is rather slow and simple implementation of calculating algorithm.
 */
{
public:
	BedMap	(const Features& bed, DefRegions& gRgns);
	
	// Calculates r and fills results
	//	@cc: type of correlation coefficient
	//	@bMap: BedMap object to correlate with
	//	@results: object to fill results
	void CalcR(CCkey::eCC ecc, BedMap& bMap, Results& results);

	//void  Write(string fileName);
};

class TestCC
{
	class Sample
	{
		static const int ArrLen = 100;
		ChromLengths _arr;
	
	public:
		Sample(const string& fname);
	
		inline void GetR(const Sample& sample)
		{ _arr.GetR(1, CCkey::ccP, sample._arr).Print(); }

		void Print()
		{ for(int i=0; i<ArrLen; i++)	cout << _arr[i] << EOL;	}
	};

public:
	TestCC();
};

#endif	// DEBUG
