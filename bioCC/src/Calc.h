/**********************************************************
Calc.h
Provides classes for calculating CC
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 04.24.2024
***********************************************************/
#pragma once

#include "bioCC.h"
#include "Options.h"
#include "DefRegions.h"
#include "Features.h"


typedef pair<double, double> pairDbl;
typedef pair<Regions::Iter, Regions::Iter> RegionsRange;

// 'eRS' defines types of printed template regions (features) CC
enum eRS {
	rsOFF = 0,		// not set; never pointed in command line
	rsR = 1,		// sorted by regions
	rsC = 2			// sorted by coefficients
};

// 'PrintMngr': print manager. Prints local & total results according to laconic verbose level
static class PrintMngr
{
public:
	enum ePrint {
		LOC = 0x1,	// output local CC
		TOT = 0x2,	// output total CC
	};

	// Initialisez by user setings
	static void Init(int ccPrint, eOInfo info, bool multiFiles);

	// Returns true if oinfo output is not laconic
	inline static bool IsNotLac() { return _OInfo > eOInfo::LAC; }

	// Returns true if name is printed
	inline static bool IsPrName () { return _PrName; }

	// Returns true if only name (without additional info) is printed
	inline static bool IsPrNameOnly() { return _PrName && _OInfo <= eOInfo::NM; }

	// Returns name ot null depending on the print name mode
	inline static const char* EchoName(const char* name) { return _PrName ? nullptr : name; }

	// Returns outputted info
	inline static eOInfo OutInfo() { return _OInfo; }

	// Throws an exception for an empty object
	static void CompleteEmpty(const char* fName, FT::eType type);

	// Returns true if local CC should be calculated
	inline static bool IsPrintLocal() { return _PrintCC & LOC; }

	// Returns true if total CC should be calculated
	inline static bool IsPrintTotal() { return _PrintCC & TOT; }

	// Prints PCC
	//	@cc: correlation coefficient
	//	@cID: calculated chrom pcc or total
	static void PrintCC(float cc, chrid cID = Chrom::UnID);

private:
	static ePrint	_PrintCC;
	static eOInfo	_OInfo;
	static bool		_PrName;
} pr;

using freqPair = pair<chrlen, UINT>;

// 'ValPos' incapsilates valued position
struct ValPos
{
	chrlen	Pos;
	float	Val;

	inline ValPos(chrlen pos = 0, float val = 0) : Pos(pos), Val(val) {}

	inline ValPos(const freqPair& vpos) : Pos(vpos.first), Val(float(vpos.second)) {}

	inline void Clear() { Pos = 0, Val = 0; }

	// Compares two Regions by start position. For sorting a container.
	// Basic abstract method implementation.
	static inline bool CompareByStartPos(const ValPos& v1, const ValPos& v2)
	{ return v1.Pos < v2.Pos;	}

#ifdef _DEBUG
	inline void Print() const { cout << Pos << TAB << Val << LF; }
#endif	// _DEBUG
};

// 'PlainCover' represents coverage as a collection of ValPos. Base class for Cover and ReadDens
class PlainCover : public Items<ValPos>
{
	// PlainCover is represented by a set of ALL valued ranges, including zero valued.
	// Thus regions span the entire chromosome space.
	// Zero valued ranges are not represented in input WIG file, they are added during initialization.

	const float _binWidth;		// width of bins of histogram; if 0, no histogram
	const eRS	 _printFRes;	// sign to print results for each feature from 'template' and how to sort it

protected:
	UniBedReader* _file = nullptr;	// for child constructor only
	size_t	_lastInd = 0;			// last index of recorded item; for child constructor only

	PlainCover() :
		_binWidth(Options::GetFVal(oBINWIDTH)),
		_printFRes(Options::Assigned(oPRFCC) ? eRS(Options::GetIVal(oPRFCC)) : eRS::rsOFF)
	{}

	// pass through file records
	template<typename T>
	void Pass(T* obj, UniBedReader& file) {
		_file = &file;
		file.Pass(*obj);
		_file = nullptr;
	}

	// Gets an item's title.
	//	@type: data type
	//	@pl: true if plural form
	//inline static const string& ItemTitle(FT::eType type, bool pl) { return FT::ItemTitle(type, pl); };

	// Adds unzero valued position to the coverage
	//	@vPos: unzero region start position and value
	//	@prevEnd: previoues unzero region end position
	//	return: true if position is added
	bool AddPos(const ValPos& vPos, chrlen prevEnd);

	// Adds chrom to the instance
	//	@cID: chrom
	//	@cLen: current chrom length
	//	@prevEnd: end of previous entry
	void AddChrom(chrlen cID, chrlen cLen, chrlen prevEnd);

public:
	// Calculates and prints corr. coefficients, using single-pass range-based algorithm
	//	@cv: compared cover
	//	@gRgns: def regions (chrom sizes)
	//	@templ: template to define treated regions
	//	return: true if calculation was actually done
	bool CalcR(const PlainCover& cv, const DefRegions& gRgns, const Features* templ);

	// Writes inner representation to BEDGRAPG file
	void Write(const string& fName) const;
};

// 'Cover' represents chrom's cover initialized by WIG data
class Cover : public PlainCover
{
	// Gets an item's title
	//	@pl: true if plural form
	//static const string& ItemTitle(bool pl = false) { return PlainCover::ItemTitle(FT::eType::BGRAPH, pl); }

	// Adds chrom to the instance
	//	@cID: chrom
	//	@cLen: current chrom length
	void AddChrom(chrlen cID, chrlen cLen) { PlainCover::AddChrom(cID, cLen, _file->PrevItemEnd()); }

	// Initializes instance from wig file
	//	return: numbers of all and initialied items for given chrom
	void InitWiggle(BedReader& file, const ChromSizes& cSizes);

public:
	// Creates new instance by wig-file name
	// Invalid instance wil be completed by throwing exception.
	//	@fName: file name
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@prfName: true if file name should be printed unconditionally, otherwise deneds on oinfo
	//	@abortInval: true if invalid instance should abort excecution
	Cover(const char* fName, ChromSizes& cSizes, eOInfo oinfo, bool prfName, bool abortInval);

	// Adds cover's item
	bool operator()()
	{ AddPos(ValPos(_file->ItemStart(), _file->ItemValue()), _file->PrevItemEnd()); return true; }

	// Closes current chrom, open next one
	//	@param cID: current chrom ID
	//	@param cLen: current chrom length
	//	@param cnt: current chrom items count
	//	@param nextcID: next chrom ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID) { AddChrom(cID, cLen); }

	// Closes last chrom
	//	@param cID: current chrom ID
	//	@param cLen: current chrom length
	//	@param cnt: current chrom items count
	//	@param nextcID: next chrom ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, size_t tCnt) { AddChrom(cID, cLen); }
};

class ReadDens : public PlainCover
{
	using rfreq = map<chrlen, UINT>;

	rfreq* _freq = nullptr;		// reads frequency; for constructor only

	// Gets an item's title.
	//	@pl: true if plural form
	//static const string& ItemTitle(bool pl = false) { return PlainCover::ItemTitle(FT::eType::ABED, pl); };

	// Adds chrom to the instance
	//	@cID: current chrom ID
	//	@cLen: current chrom length
	void AddChrom(chrid cID, chrlen cLen);

public:
	// Creates new instance by abed/bam-file name
	// Invalid instance wil be completed by throwing exception.
	//	@fName: file name
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@printfName: true if file name should be printed unconditionally, otherwise deneds on oinfo
	//	@abortInval: true if invalid instance should abort excecution
	ReadDens(const char* fName, ChromSizes& cSizes, eOInfo oinfo, bool printfName, bool abortInval);

	// Adds Read
	bool operator()() {
		(*_freq)[_file->ItemStrand() ? _file->ItemStart() : _file->ItemEnd()]++;
		return true;
	}

	// Closes current chrom, open next one
	//	@param cID: current chrom ID
	//	@param cLen: current chrom length
	//	@param cnt: current chrom items count
	//	@param nextcID: next chrom ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID) { if(cnt) AddChrom(cID, cLen); }

	// Closes last chrom
	//	@param cID: last chrom ID
	//	@param cLen: current chrom length
	//	@param cnt: last chrom items count
	//	@param tCnt: total items count
	void operator()(chrid cID, chrlen cLen, size_t cnt, size_t tCnt) { if (cnt) AddChrom(cID, cLen); }
};

// 'ChromRanges' represented a set of chromosome's ranges.
struct ChromRanges : public ItemIndices
	/*
		To keep chromosome's number and first/last ranges indexes Features::ChromItemsInd is used;
		FirstInd/LastInd are keeping first/last ranges indexes,
		and lengths of all features from fs1, fs2 are keeping additionally for efficiency.
	*/
{
	chrlen FeatrsLen1;		// length of all features of chromosome in fs1
	chrlen FeatrsLen2;		// length of all features of chromosome in fs2

	inline ChromRanges() : ItemIndices(), FeatrsLen1(0), FeatrsLen2(0) {}

	inline ChromRanges(chrlen firstInd, chrlen lastInd, chrlen len1, chrlen len2) :
		ItemIndices(firstInd, lastInd), FeatrsLen1(len1), FeatrsLen2(len2) {}
};

// 'JointedBeds' represents two bed-files as a chromosomes collection and theirs joint features (ranges).
//	Fast but a bit complicated implementation of calculating algorithm.
class JointedBeds : Chroms<ChromRanges>
{
private:
	// 'Range' represens a joint feature with start position and joint (combined) value.
	struct Range
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

	vector<Range> _ranges;

public:
	// Fills ChromRanges & Range by given two beds.
	// Beds chromosomes should be checked as Treated.
	// Both fs1 & fs2 must be valid: no duplicated, crossed, adjacent, coverage features;
	// in other case R may be wrong
	JointedBeds(const Features& fs1, const Features& fs2);

	// Calculates r and fills results
	//	@cSizes: chrom sizes
	//	return: true if calculation was actually done
	bool CalcR(const ChromSizes& cSizes);

#ifdef _DEBUG
	void	Print();
#endif
};


// 'CorrPair' represents pair of objects to compare, and methods for recognizing types and calculation CC
class CorrPair
{
	// Pointers to the common methods
	struct FileType {
		void* (CorrPair::* Create)(const char*, bool);	// type constructor
		void(*Delete)(void*);							// type destructor
	};

	static FileType _FileTypes[];

	void* _firstObj = nullptr;
	void* _secondObj = nullptr;
	Features* _templ = nullptr;
	DefRegions& _gRgns;		// initial genome regions to correlate
	FT::eType	_type;		// type of compared files
	BYTE	_typeInd;		// type of file types: 0 - wig, 1 - bedF, 2 - bedR

	// Returns true if Features are treating
	inline bool IsBedF()	const { return _typeInd == 1; }

	// Throws an exception if no items, or prints LF while needed
	template<typename T>
	void CheckItemsCount(T* obj, const char* fName)	{
		if (!obj->ItemsCount()) 
			PrintMngr::CompleteEmpty(fName, _type);
		else if (PrintMngr::IsPrNameOnly())
			dout << LF;
	}

	// Creates features bed object.
	//	@fName: file name
	//	@primary: if true object is primary
	void* CreateBedF(const char* fName, bool isPrimary);

	// Creates alignment object
	//	@fName: file name
	//	@primary: if true object is primary
	void* CreateBedR(const char* fName, bool isPrimary);

	// Creates covering object.
	//	@fName: file name
	//	@primary: if true object is primary
	void* CreateWig(const char* fName, bool isPrimary);

	// Calculates r for genome features.
	//	@firstBed: first Features to correlate
	bool CalcCCBedF(Features& first) {
		// common chroms are set automatically by JointedBeds
		return JointedBeds(first, *((Features*)_secondObj)).CalcR(_gRgns.ChrSizes());
	}

	// Calculates r for coverages and read densities
	bool CalcCCCover(DefRegions& gRgns) {
		return ((Cover*)_firstObj)->CalcR(*((Cover*)_secondObj), gRgns, _templ);
	}

	// Checks file extisting and extention validity
	//	@fName: file's name
	//	@abortInvalid: if true throw extention if checking is false
	//	return: index in _FileTypes or _FileTypesCnt if checking is false
	char CheckFileExt(const char* fName, bool abortInvalid);

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
