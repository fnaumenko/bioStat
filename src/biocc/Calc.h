/**********************************************************
Calc.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 3.12.2020
-------------------------
Provides classes for calculating CC
***********************************************************/
#pragma once

#include "Data.h"
#include "bioCC.h"

typedef pair<double, double> pairDbl;
typedef pair<Regions::Iter, Regions::Iter> RegionsRange;

// 'eRS' defines types of printed template regions (features) CC
enum eRS {
	rsOFF = 0,		// not set; never pointed in command line
	rsR = 1,		// sorted by regions
	rsC = 2			// sorted by coefficients
};

// 'PSums' accumulates sums in single-pass range-based algorithm and performes CC calculation
class PSums 
{
	double _sumX = 0, _sumY = 0;		// sum of signal values
	double _sumXY = 0;					// sum of the products of the values of both signals
	double _sumSqrX = 0, _sumSqrY = 0;	// sum of squared signal values
	float	_valX = 0, _valY = 0;
	genlen	_len = 0;

public:
	void Clear() { _len=_sumX=_sumY=_sumXY=_sumSqrX=_sumSqrY=_valX=_valY=0; }

	// Adds range length and correlated range values
	void AddVal(chrlen len, float valX, float valY);

	// Returnes Pearson CC
	double PCC();

	// Returnes Signal CC
	double SCC() { return _sumXY / sqrt(_sumSqrX) / sqrt(_sumSqrY); }
};

// 'CC' represents a pair of corr. coeficients: Pearson(first) and signal(second), and provides theirs output.
struct CC : private pairDbl
{
public:
	enum eCC {	// defines types of CC
		ccP = 0x1,	// Pearson correlation coefficient: first bit
		ccS = 0x2,	// signal correlation coefficient: second bit
	};

private:
	static const int _Empty = -2;	// uninitialised coefficient
	static eCC Stated;				// setting for the current session

	static void Print(double val);

public:
	static const int Undef = 666;	// undefined coefficient

	static inline void Set(eCC stated) { Stated = stated; }

	// Returns true if parameter is signal coefficient
	static inline bool IsS() { return Stated & ccS; }

	// Returns true if parameter is Pearson coefficient
	static inline bool IsP() { return Stated & ccP; }

	// Returns true if parameter is both coefficients
	static inline bool IsBoth() { return Stated & ccS && Stated & ccP; }

	inline CC() { first = second = _Empty; }

	inline CC(double p, double s) { first = p; second = s; }

	inline bool operator < (const CC& ccres) const
	{ return first != _Empty ? first < ccres.first : second < ccres.second; }

	// Sets Pearson coefficient
	inline double GetP() const { return first; }

	// Sets signal coefficient; undefined value by default
	inline double GetS() const { return second; }

	// Sets Pearson coefficient
	inline void SetP(double val) { first = val; }

	// Sets signal coefficient; undefined value by default
	inline void SetS(double val = 0) { second = val; }

	// Sets CC by sums in one pass algorithm
	void Set(PSums& sums);

	// Returns true if even one value in pair is setting
	inline bool NotEmpty() const { return first != _Empty || second != _Empty; }

	// returns single value
	inline double GetSingleVal() const { return first != _Empty ? first : second; }

	// set single negative value to absolute value
	void SetSingleAbsVal();

	void Print() const;
};

// 'PrintMngr': print manager. Prints local & total results according to laconic verbose level
static class PrintMngr
{
public:
	enum ePrint {
		pLOC = 0x1,	// output local CC
		pTOT = 0x2,	// output total CC
	};

private:
	static ePrint PrintCC;
	static Obj::eInfo	Verbose;

public:
	static void Init(int ccPrint, Obj::eInfo verb);

	// Returns true if info output is not laconic
	inline static bool IsNotLac() { return Verbose > Obj::eInfo::LAC; }

	// Returns true if local CC should be calculated
	inline static bool IsPrintLocal() { return PrintCC & pLOC; }

	// Returns true if total CC should be calculated
	inline static bool IsPrintTotal() { return PrintCC & pTOT; }

	// Returns verbose level
	inline static Obj::eInfo Verb() { return Verbose; }

	static void PrintChrom(chrid cID, const CC& cc);

	static void PrintTotal(const CC& cc);
} pr;

// 'ValPos' incapsilates valued position
struct ValPos
{
	chrlen	Pos;
	float	Val;

	inline ValPos(chrlen pos, float val) : Pos(pos), Val(val) {}

	// Compares two Regions by start position. For sorting a container.
	// Basic abstract method implementation.
	static inline bool CompareByStartPos(const ValPos& v1, const ValPos& v2)
	{ return v1.Pos < v2.Pos;	}

#ifdef _DEBUG
	inline void Print() const { cout << Pos << TAB << Val << LF; }
#endif	// _DEBUG
};

// 'BaseCover' represents base class for Cover and ReadDens
class BaseCover : public Items<ValPos>
{
	// BaseCover is represented by a set of ALL valued ranges, including zero valued.
	// Thus regions span the entire chromosome space.
	// Zero valued ranges are not represented in input WIG file, they are added during initialization.

	double	_binWidth;	// width of bins of histogram; if 0, no histogram
	eRS		_printFRes;	// sign to print results for each feature from 'template' and how to sort it

	// Checks the element for the new potential start/end positions for all possible ambiguous.
	// Items<> abstract method implementation
	//  return: true if item should be accepted; otherwise false
	inline bool CheckPrevPos(const Region&, ItemsIter, Spotter&) { return true; }

	// Gets a copy of region by container iterator.
	// Items<> abstract method implementation.
	inline Region const Regn(cItemsIter it) const { return Region(); }

protected:
	BaseCover() :
		_binWidth(Options::GetFVal(oBINWIDTH)),
		_printFRes(Options::Assigned(oPRFCC) ? eRS(Options::GetIVal(oPRFCC)) : eRS::rsOFF)
	{}

	// Adds unzero valued position to the coverage
	//	@start: unzero region start position
	//	@prevEnd: previoues unzero region end position
	//	@val: unzero region value
	void AddPos(chrlen start, chrlen prevEnd, float val);

public:
	// Calculates and prints corr. coefficients, using single-pass range-based algorithm
	//	@cv: compared cover
	//	@gRgns: def regions (chrom sizes)
	//	@templ: template to define treated regions
	void CalcR(const BaseCover& cv, const DefRegions& gRgns, const Features* templ);
};

// 'Cover' represents chrom's cover
class Cover : public BaseCover
{
	// Adds frag/read to the container.
	// Abstract BaseItems<> method implementation.
	//	@rgn: Region with mandatory fields
	//	@spotter: temporary values & ambiguities
	//	return: true if frag/read was added successfully
	bool AddItem(const Region& rgn, Spotter& spotter);

	// Adds last zero range to close the coverage.
	//	Abstract BaseItems<> method implementation 
	//	@spotter: used do get last item end
	//	return: count of added items
	UINT FinishItems(const Spotter& spotter);

	// Gets an item's title
	//	Obj abstract method implementation.
	//	@pl: true if plural form
	inline const string& ItemTitle(bool pl = false) const
	{ return FT::ItemTitle(FT::eType::BGRAPH, pl); };

	// Items<> method redefinition 
	inline void SortIfNecessary(Spotter&, ULONG, bool) {}

	// Initializes instance from wig file
	//	@spotter: spotter to control ambiguities
	//	@cSizes: chrom sizes to control chrom length exceedeing
	//	return: numbers of all and initialied items for given chrom
	p_ulong InitDerived(Spotter& spotter, const ChromSizes& cSizes);

public:
	// Creates new instance by wig-file name
	// Invalid instance wil be completed by throwing exception.
	//	@title: title printed before file name or NULL
	//	@fName: file name
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@info: type of feature ambiguties that should be printed
	//	@printfName: true if file name should be printed unconditionally, otherwise deneds on info
	//	@abortInval: true if invalid instance should abort excecution
	Cover(const char* title, const string& fName, ChromSizes& cSizes, 
		eInfo info, bool printfName, bool abortInval) : BaseCover()
	{
		Spotter spotter(FT::eType::BGRAPH, info);
		Init(title, fName, spotter, cSizes, info > eInfo::LAC || printfName, abortInval, 4);
	}
};

// 'ReadDens' represents Read density as a cover
class ReadDens : BasicReads, public BaseCover
{
	map<chrlen, UINT> _map;

	// Adds read to the container.
	// Abstract BaseItems<> method implementation.
	//	@rgn: Region with mandatory fields
	//	@spotter: temporary values & ambiguities
	//	return: true if read was added successfully
	bool AddItem(const Region& rgn, Spotter& spotter);

	// Items<> method empty redefinition
	inline bool CheckLastPos(const Region&, Spotter&) {	return true; }

	// Abstract BaseItems<> method empty implementation 
	inline UINT FinishItems(const Spotter&)	{ return 0;	}

	// Gets an item's title.
	//	Obj abstract method implementation.
	//	@pl: true if plural form
	inline const string& ItemTitle(bool pl = false) const
	{ return FT::ItemTitle(FT::eType::ABED, pl); };

	// Fills items from intermediate container
	//	BaseItems<> abstract method implementation.
	void FillChromItems(const Spotter& spotter);

public:
	// Creates new instance by abed/bam-file name
	// Invalid instance wil be completed by throwing exception.
	//	@title: title printed before file name or NULL
	//	@fName: file name
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@info: type of feature ambiguties that should be printed
	//	@printfName: true if file name should be printed unconditionally, otherwise deneds on info
	//	@abortInval: true if invalid instance should abort excecution
	ReadDens(const char* title, const string& fName, ChromSizes& cSizes,
		eInfo info, bool printfName, bool abortInval) : BaseCover()
	{
		Spotter spotter(FT::GetType(fName.c_str(), true), info);
		Init(title, fName, spotter, cSizes, info > eInfo::LAC || printfName, abortInval);
	}
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
	// Class 'R' accumulates variances & covariance and calculated r
	class R {
	private:
		double	_var1, _var2;	// variances 
		double	_cov;			// covariance: SUM( X-Xmean)^2 * (Y-Ymean)^2 )
		/*
		Mean value for every range may accept only one of two values: 0-mean or 1-mean,
		and they are static within chromosome.
		So for efficiency we can keep theirs in 3 arrays:
		two square's arrays: arrays of squared values from each bed;
		there are only 2 combinations, but they are duplicated for the simplicity acces,
		and one crossing array: array of all combinations
		of multiplications different values from fs1 and fs2
		*/
		double	_sqMeans1[4], _sqMeans2[4];	// square's arrays for bed, fs2
		double	_crossMeans[4];				// crossing array

		double CalcR() const {
			return  !_var1 && !_var2 ? CC::Undef : _cov / (sqrt(_var1) * sqrt(_var2));
		}	// not sqrt(_var1*_var2) 

	public:
		inline R() : _var1(0), _var2(0), _cov(0) {}

		// Keeps mean values for fs1, fs2. 
		//	@clear: if true, clear instance for treatment of new chromosome
		void Init(double mean1, double mean2, bool clear);

		// Accumulates next length of range
		void Increment(chrlen len, char val);

		// Returns coefficient of correlation
		inline double Get() const { return CalcR(); }
	};

	// 'PairR' pepresetns a pair of R classes accumulated variances & covariance and calculated r.
	// First in pair is signal R, second - Pearson R.
	struct PairR : pair<R, R>
	{
		// Keeps mean values for fs1, fs2. 
		//	@clear: if true, clear instance for treatment of new chromosome
		void Init(double mean1, double mean2, bool clear) {
			if (CC::IsS())	first.Init(0, 0, clear);
			if (CC::IsP())	second.Init(mean1, mean2, clear);
		}

		// Accumulates next length of range
		void Increment(chrlen len, char val) {
			if (CC::IsS())	first.Increment(len, val);
			if (CC::IsP())	second.Increment(len, val);
		}

		// Returns coefficient of correlation
		CC Get() const {
			CC res;
			if (CC::IsS())	res.SetS(first.Get());
			if (CC::IsP())	res.SetP(second.Get());
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

		inline Range(chrlen pos = 0, char val = 0) : Start(pos), Val(val) {}
	};

	vector<Range>	_ranges;

public:
	// Fills ChromRanges & Range by given two beds.
	// Beds chromosomes should be checked as Treated.
	// Both fs1 & fs2 must be valid: no duplicated, crossed, adjacent, coverage features;
	// in other case R may be wrong
	JointedBeds(Features& fs1, Features& fs2);

	// Calculates r and fills results
	//	@cSizes: chrom sizes
	void CalcR(const ChromSizes& cSizes);

#ifdef _DEBUG
	void	Print();
#endif
};


// 'CorrPair' represents pair of objects to compare, and methods for recognizing types and calculation CC
class CorrPair
{
private:
	struct FileType	// keeps pointers to the common methods
	{
		typedef Obj* (CorrPair::* CreateObj)	(const char*, bool);
		typedef void (*DeleteObj)				(Obj*);		// static function

		CreateObj		Create;			// constructor of type
		DeleteObj		Delete;			// destructor of type
	};

	static int _FileTypesCnt;		// count of treatment file's types
	static FileType _FileTypes[];

	Obj* _firstObj;
	Obj* _secondObj;
	DefRegions& _gRgns;		// initial genome regions to correlate
	Features* _templ;
	BYTE	_typeInd;		// type of file types: 0 - wig, 1 - bedF, 2 - bedR
	bool	_printWarn;		// true if print line warnings
	bool	_printName;		// print file names unconditionally

	// Returns true if Reads are treating
	inline bool IsWig()		const { return _typeInd == 0; }

	// Returns true if Features are treating
	inline bool IsBedF()	const { return _typeInd == 1; }
	
	// Returns true if BedE are treating
	inline bool IsAlign()	const { return _typeInd == 2; }

	// Creates features bed object.
	Obj* CreateBedF(const char* fName, bool isPrimary);

	// Creates alignment object
	Obj* CreateBedR(const char* fName, bool isPrimary) {
		return new ReadDens(NULL, fName, _gRgns.ChrSizes(), PrintMngr::Verb(), _printName, isPrimary);
	}

	// Creates covering object.
	//	@fName: file name
	//	@primary: if true object is primary
	Obj* CreateWig(const char* fName, bool isPrimary) {
		return new Cover(NULL, fName, _gRgns.ChrSizes(), PrintMngr::Verb(), _printName, isPrimary);
	}

	// Calculates r for genome features.
	//	@firstBed: first Features to correlate
	inline void CalcCCBedF(Features& firstBed) {
		// common chroms are set automatically by JointedBeds
		JointedBeds(firstBed, *((Features*)_secondObj)).CalcR(_gRgns.ChrSizes());
	}

	// Calculates r for coverages and read densities
	void CalcCCWigCover(DefRegions& gRgns) {
		((Cover*)_firstObj)->CalcR( *((Cover*)_secondObj), gRgns, _templ );
	}

	// Checks file extisting and extention validity
	//	@fName: file's name
	//	@abortInvalid: if true throw extention if checking is false
	//	return: index in _FileTypes or _FileTypesCnt if checking is false
	BYTE CheckFileExt(const char* fName, bool abortInvalid);

public:
	// Creates an instance with checking primary object.
	//	@primefName: primary file's name
	//	@gRgns: genome regions
	//	@templ: name of template bed file, or NULL if undefined
	//	@multiFiles: true if more then one secondary files are placed
	CorrPair(const char* primefName, DefRegions& gRgns, const char* templ, bool multiFiles);

	~CorrPair();

	// Adds secondary object, calculates and prints CCkey.
	void CalcCC(const char* fName);
};
