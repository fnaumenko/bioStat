/**********************************************************
Calc.ccp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 3.12.2020
-------------------------
Provides classes for calculating CC
***********************************************************/

#include "Calc.h"
#include "Data.h"

const char* AND = " and";
const string ForCorrelation = " for correlation";
const string sFormat = " format";

/********************  PSums *********************/

// Adds range length and correlated range values
void PSums::AddVal(chrlen len, float valX, float valY)
{
	const double valXLen = valX * len, valYLen = valY * len;

	_sumX += valXLen;	_sumSqrX += valX * valXLen;
	_sumY += valYLen;	_sumSqrY += valY * valYLen;
	_sumXY += valXLen * valY;
	_len += len;
	//cout << ":\tlen: " << len << TAB << valX << TAB << valY;// << LF;
	//cout << "\t_len: " << _len << LF;
	//cout << "\tsumX: " << _sumX << "\tsumY: " << _sumY << "\tsumSqrX: " << _sumSqrX << "\tsumSqrY: " << _sumSqrY << LF;
}

// Returnes Pearson CC
double PSums::PCC()
{
	double pcc = (_len * _sumXY - _sumX * _sumY) /
		sqrt(_len * _sumSqrX - _sumX * _sumX) /
		sqrt(_len * _sumSqrY - _sumY * _sumY);
	return pcc;

	//return (_len * _sumXY - _sumX * _sumY) /
	//	sqrt(_len * _sumSqrX - _sumX * _sumX) /
	//	sqrt(_len * _sumSqrY - _sumY * _sumY);
	// not sqrt((_len * _sumSqrX - _sumX * _sumX) * (_len * _sumSqrY - _sumY * _sumY)); because of possinle overwlov
}

/********************  end of PSums *********************/

/********************  ÑÑ *********************/

CC::eCC CC::Stated;		// setting for the current session

void CC::Print(double val)
{
	if (val == _Empty)	return;
	dout << left << setw(12) << setfill(' ');
	if (val == Undef || _isnan(val))	dout << "UNDEF";
	else	dout << val;
}

// Sets CC by sums in one pass algorithm
void CC::Set(PSums& sums)
{
	if (IsP())	SetP(sums.PCC());
	if (IsS())	SetS(sums.SCC());
}

// set single negative value to absolute value
void CC::SetSingleAbsVal()
{
	if (first != _Empty)	first = fabs(first);
	else					second = fabs(second);
}

void CC::Print() const
{
	Print(first);
	Print(second);
	dout << LF;
}

/********************  end of ÑÑ *********************/

/********************  PrintMngr *********************/

PrintMngr::ePrint	PrintMngr::PrintCC;
Obj::eInfo	PrintMngr::Verbose;

void PrintMngr::Init(int ccPrint, Obj::eInfo verb)
{
	PrintCC = ePrint(ccPrint);
	Verbose = verb;
}

void PrintMngr::PrintChrom(chrid cID, const CC& cc)
{
	dout << Chrom::AbbrName(cID) << TAB;
	cc.Print();
	fflush(stdout);		// when called from a package 
}

void PrintMngr::PrintTotal(const CC& cc)
{
	if (Verbose > Obj::eInfo::LAC || IsPrintLocal())
		dout << "total:\t";
	cc.Print();
	fflush(stdout);		// when called from a package
}

/********************  end of PrintMngr *********************/

/************************ FeatureR ************************/

// 'FeatureR' represetns pair <feature-ID><feature-CC>
struct FeatureR : pair<chrlen, CC>
{
	inline FeatureR(chrlen i, const CC& res) { first = i; second = res; }

	inline bool operator < (const FeatureR& rccr) const { return second < rccr.second; }

	// returns single value
	inline double GetSingleVal() const { return second.GetSingleVal(); }

	// set single negative value to absolute value
	inline void SetSingleAbsVal() { second.SetSingleAbsVal(); }

	inline void Print() const { dout << first << TAB; second.Print(); }
};

// 'FeatureRs' represetns FeatureR collection, including methods to print collection and CC histogram
class FeatureRs : vector<FeatureR>
{
private:
	// Replace negative values by positive
	//	return: true if even one value had been replaced
	bool SetAbsVals()
	{
		bool holdNegative = false;
		for (vector<FeatureR>::iterator it = begin(); it != end(); it++)
			if (it->GetSingleVal() < 0) {
				it->SetSingleAbsVal();
				holdNegative = true;
			}
		return holdNegative;
	}

public:
	inline FeatureRs(chrlen cnt) { reserve(cnt); }

	inline void AddVal(chrlen ind, CC val) { push_back(FeatureR(ind + 1, val)); }

	void Print(int printFRes)
	{
		if (printFRes == rsOFF)	return;
		if (printFRes == rsC)		// soretd by feature; are sorted initially
			sort(begin(), end());	// by increase
		dout << "\n#RGN\tCC\n";
		for (vector<FeatureR>::iterator it = begin(); it != end(); it++)
			it->Print();
	}

	// Creates and prints histogram
	void PrintHist(double binWidth)
	{
		if (!binWidth)		return;
		SetAbsVals();
		sort(begin(), end());	// by increase
		// define factor
		// factor is a divisor of binWidth: 0.1--0.9=>10, 0.01--0.09=>100 etc
		short F = 10;
		for (; binWidth * F < 1; F *= 10);
		vector<FeatureR>::iterator it = begin();
		// then float instead of double because of wrong consolidation by round double
		double minBin = float(int(F * it->GetSingleVal())) / F;
		//float minBin = F*it->GetSingleVal();
		//		minBin = float(int(minBin))/F;
		double maxBin = F * (end() - 1)->GetSingleVal();
		int	maxdecBin = int(maxBin);
		if (maxBin - maxdecBin)	maxdecBin++;	// round up
		if (maxdecBin % 2)			maxdecBin++;	// get even bin
		maxBin = double(maxdecBin) / F;
		Array<int> hist(int((maxBin - minBin) / binWidth) + 1);		// histogram
		// consolidation: fill histogram
		for (; it != end(); it++)
			hist[int((maxBin - it->GetSingleVal()) / binWidth)]++;
		// print histogram
		dout << "BIN UP\tCOUNT\n";
		for (BYTE k = 0; k < hist.Length(); k++)
			dout << (maxBin - k * binWidth) << TAB << hist[k] << LF;
		//dout << "finish\n";
	}
};

/************************ end of FeatureR ************************/

/************************ class BaseCover ************************/

// Adds unzero valued position to the coverage
//	@start: unzero region start position
//	@prevEnd: previoues unzero region end position
//	@val: unzero region value5
void BaseCover::AddPos(chrlen start, chrlen prevEnd, float val)
{
	// is it a gap between prev and current regions?
	if (start != prevEnd && val)				// check for the unzero val to avoid duplicate zero items
		_items.push_back(ValPos(prevEnd, 0));	// add zero space after previous unzero region
	if (_items.back().Val != val)				// check prev val for wiggle_0 with regular tiems with the same span (MACS output
		_items.push_back(ValPos(start, val));	// add start of current unzero region
}

// Sets new position and returns current len
chrlen SetPos(chrlen& pos, chrlen currPos)
{
	const chrlen len = currPos - pos;
	pos = currPos;
	return len;
}

// Calculates and prints corr. coefficients, using single-pass range-based algorithm
//	@cv: compared cover
//	@rgns: def regions (chrom sizes)
//	@templ: template to define treated regions
void BaseCover::CalcR(const BaseCover& cv, const DefRegions& rgns, const Features* templ)
{
	CC totalCC;
	PSums	totalSums;
	bool fillLocRes = templ && (_binWidth || _printFRes >= 0);
#ifdef _DEBUG
	//Print("first", 15);
	//cv.Print("second", 15);	cout << LF;
	//rgns.Print();
	//if (templ)	templ->Print();
#endif
	// rgns is already limited by chroms represented in template, if it's defined
	for (DefRegions::cIter rit = rgns.cBegin(); rit != rgns.cEnd(); rit++) {
		CC chrCC;
		chrlen fCnt = 0;	// count of features
		chrid cID = CID(rit);
		if (!FindChrom(cID) || !cv.FindChrom(cID))	continue;
		Items<Featr>::cItemsIter itF, itFend;		// template feature iterator
		if (templ) {
			auto itC = templ->GetIter(cID);
			itF = templ->ItemsBegin(itC);
			itFend = templ->ItemsEnd(itC);
			fCnt = templ->ItemsCount(itC);
		}

		// local results
		CC locCC;
		chrlen i = 0;
		FeatureRs locResults(_binWidth ? fCnt : 0);	// create histogram

		// X, Y denotes values corresponding to the first and second sequence
		PSums cSums, locSums;
		chrlen len;					// length of current united (with equal X, Y values) region
		chrlen posN, pos = 0;		// position keeps posX or posY, current position, 
		float valX = 0, valY = 0;	// X, Y current value
		bool insideF = !templ;		// true if current position is inside current template feature
		bool closeF = false;		// true if the feature has just ended
		bool inTempl = templ;		// true if current position did not go beyond the border of the last feature 
		cItemsIter itX = ItemsBegin(cID), itY = cv.ItemsBegin(cID);
		const cItemsIter itXend = ItemsEnd(cID), itYend = cv.ItemsEnd(cID);

		while (itX != itXend && itY != itYend) {
			const chrlen posX = itX->Pos,	posY = itY->Pos;
			const float prevValX = valX,	prevValY = valY;	// X, Y current value

			//== set valX, valY
			if (posX > posY)			// more
				posN = posY, valY = itY++->Val;
			else {						// equal or less
				posN = posX; valX = itX++->Val;
				if (posX == posY) 	// equal
					valY = itY++->Val;
			}
			//== set pos, len
			if (inTempl)
				if (posN > itF->End) {							// exit current feature
					closeF = len = SetPos(pos, itF->End+1);		// take len including the end of the feature
					inTempl = ++itF != itFend;					// next feature
				}
				else if (!insideF && posN >= itF->Start)		// entry current feature
					insideF = len = (pos = posN) - itF->Start;	// cutoff range before feature start
				else len = SetPos(pos, posN);
			else len = SetPos(pos, posN);

			//== accumulate sums; the last interval, ended by pos, is saved!!
			if (insideF && len) {		// len can be 0 at the boundary of the feature
				//cout << pos;
				if (fillLocRes)		locSums.AddVal(len, prevValX, prevValY);
				cSums.AddVal(len, prevValX, prevValY);		// previous combined region
				if (PrintMngr::IsPrintTotal())
					totalSums.AddVal(len, valX, valY);
			}

			//== close feature, save its CC
			if (closeF) {
				if (fillLocRes) {			// add CC for each feature
					locCC.Set(locSums);
					locResults.AddVal(i++, locCC);
					locSums.Clear();
				}
				closeF = insideF = false;
			}
		}
		//== print current result
		chrCC.Set(cSums);
		if (PrintMngr::IsPrintLocal()) {
			if (templ)
				locResults.Print(_printFRes),
				locResults.PrintHist(_binWidth);	// print histogram
			PrintMngr::PrintChrom(cID, chrCC);
		}
	}
	if (PrintMngr::IsPrintTotal())	
		totalCC.Set(totalSums),
		PrintMngr::PrintTotal(totalCC);
}

/************************ end of BaseCover ************************/

/************************ class Cover ************************/

// Adds frag/read to the container.
// Abstract BaseItems<> method implementation.
//	@rgn: Region with mandatory fields
//	@spotter: temporary values & ambiguities
//	return: true if frag/read was added successfully
bool Cover::AddItem(const Region& rgn, Spotter& spotter)
{
	AddPos(rgn.Start, spotter.lastEnd, spotter.File().ItemValue());
	spotter.lastEnd = rgn.End;
	return true;
}

// Adds last zero range to close the coverage.
//	Abstract BaseItems<> method implementation 
//	@spotter: used do get last item end
//	return: count of added items
UINT Cover::FinishItems(const Spotter& spotter)
{
	_items.push_back(ValPos(spotter.lastEnd, 0));	// add zero region after last unzero one
	_items.push_back(ValPos(spotter.chrLen - 1, 0));// add chrom boundary to finish iterating in CalcR correctly
	return 2;
}

// Returns a pointer to the substring defined by key.
//	@str: null-terminated string to search the key
//	@key: null-terminated string to search for
//	return: a pointer to the substring after key, or NULL if key does not appear in str
const char* KeyStr(const char* str, const char* key)
{
	const char* strKey = strstr(str, key);
	return strKey ? (strKey + strlen(key)) : NULL;
}

// Checks definition or declaration line for key
//	@str: null-terminated string to search the key
//	@key: null-terminated string to search for
//	@file: file to print error message with line number
//	return: point to substring followed after the key
const char* CheckSpec(const char* str, const char* key, const TabFile& file)
{
	const char* strKey = KeyStr(str, key);
	if (!strKey)
		Err(string("absent or wrong '") + key + "' key", file.LineNumbToStr().c_str()).Throw();
	return strKey;
}

// Returns required int value with check
//	@str: null-terminated string to search the key
//	@key: null-terminated string to search for
//	@file: file to print error message with line number
//	return: key value, or throws an exception if key does not appear in str
inline chrlen GetIntKey(const char* str, const char* key, const TabFile& file)
{
	return atoi(CheckSpec(str, key, file) + 1);
}

// Returns int value
//	@str: null-terminated string to search the key
//	@key: null-terminated string to search for
//	return: key value, or 0 if key does not appear in str
chrlen GetIntKey(const char* str, const char* key)
{
	const char* line = KeyStr(str, key);
	return line ? atoi(line + 1) : 0;
}

// Initializes instance from wig file
//	@spotter: spotter to control ambiguities
//	@cSizes: chrom sizes to control chrom length exceedeing
//	return: numbers of all and initialied items for given chrom
p_ulong Cover::InitDerived(Spotter& spotter, const ChromSizes& cSizes)
{
	const char* keyTrackType = "track type=";
	const char* typeBGraph = "bedGraph";
	const char* typeWiggle = "wiggle_0";
	const char* keyVarStep = "variableStep";
	const char* keyFixStep = "fixedStep";

	BedInFile& file = (BedInFile&)spotter.File();
	const char* line = file.GetNextLine(false);			// for WIG the first line should be a track definition line

	line = CheckSpec(line, keyTrackType, file);		// check track type key
	const chrlen len = strchr(line, SPACE) - line;	// the length of wiggle type in definition line
	if (!len)	file.ThrowExcept("track type is not specified");
	if (!strncmp(line, typeBGraph, len))			// BedGraph
		return InitBed(spotter, cSizes);
	else if (!strncmp(line, typeWiggle, len)) {		// fixed or variable step
		line = file.GetNextLine(false);
		if (KeyStr(line, keyFixStep))				// fixed step
			file.ResetWigType(FT::eType::WIG_FIX, 0, BYTE(strlen(keyFixStep)) + 1);
		else if (KeyStr(line, keyVarStep))			// variableStep
			file.ResetWigType(FT::eType::WIG_VAR, 1, BYTE(strlen(keyVarStep)) + 1);
		else file.ThrowExcept(string(line) + ": absent or unknown wiggle data format");
	}
	else file.ThrowExcept("type '" + string(line, len) + "' does not supported");

	ULONG estItemCnt = file.EstItemCount();
	if (!estItemCnt)		return make_pair(0, 0);

	const char* keyChrom = "chrom";
	const char* keyStart = "start";
	const char* keyStep = "step";
	const char*	keySpan = "span";

	chrid	cID = Chrom::UnID;		// current chrom ID
	bool	skipChrom = false;		// if true skip data lines for current chrom
	chrlen	pos = 0, 				// region start position, previous region start position, 
			span = 1, step = 0,		// span spec, step spec (for fixedStep),
			firstInd = 0;			// first index in item's container for current chrom
	const bool fixedStep = file.Type() == FT::eType::WIG_FIX;

	Items::ReserveItems(estItemCnt);
	//estItemCnt = 0;	// then used as counter	
	do
		if (isdigit(*line)) {	// *** data line
			//estItemCnt++;
			if (skipChrom)	continue;
			if (fixedStep)	pos += step;
			else			pos = file.IntField(0);
			AddPos(pos, spotter.lastEnd, file.ItemValue());
			spotter.lastEnd = pos + span;	// use spotter.lastEnd instead of local variable becouse of calling it in FinishItems()
		}
		else {					// *** declaration line
			//== set specifications
			//=== check keyChrom
			line = file.ChromMark() - strlen(Chrom::Abbr);			// level the initial value of chrom mark position
			const char* line1 = CheckSpec(line, keyChrom, file);	// keyChrom initial position
			BYTE spaceCnt = line1 - line + 1;			// number of extra spaces before keyChrom (+1 for '='
			if (fixedStep) {							// fix declarative parameters
				//=== check keyStart & keyStep
				pos = GetIntKey(line1, keyStart, file);	// initial position
				line1 += strlen(keyStart) + 1;			// shift to scan the rest of the line faster
				step = GetIntKey(line1, keyStep, file);	// initial step
				line1 += strlen(keyStep) + 1;			// shift to scan the rest of the line faster
				pos -= step;		// shift 'back' before the first pass, where pos+=step will be invoke
			}
			if (!(span = GetIntKey(line1, keySpan)))	// initial span: both for fixed- and variableStep
				span = 1;
			//== check chrom
			if (file.GetNextChrom(spaceCnt)) {
				chrid nextCID = file.GetChrom();
				if (skipChrom = nextCID == Chrom::UnID)
					continue;								// negligible next chrom
				if (Chrom::NoCustom()) {					// are all chroms specified?
					if (cID != Chrom::UnID)					// skip first pass, while cutt chrom is still undefined
						AddChrom(cID, firstInd, spotter);
				}
				else {										// single chrom is specified
					if (!fixedStep && pos)		// rigion is initialized: items for the specified chrom are existed and saved
						break;		// the chrom itself will be saved after loop
					if (skipChrom = nextCID != Chrom::CustomID())
						continue;
				}
				cID = nextCID;
				spotter.lastEnd = 0;
				firstInd = ItemsCount();
				spotter.SetTreatedChrom(cID);
				if (cSizes.IsFilled())	spotter.chrLen = cSizes[cID];
			}
		}
	while (line = file.GetNextLine(false));
	// save last chrom
	if (cID != Chrom::UnID)		// is last chrom valid?
		AddChrom(cID, firstInd, spotter);
	//cout << " est/fact: " << float(estItemCnt) / ItemsCount() << SPACE;
	return make_pair(ItemsCount(), ItemsCount());	// option i=CNT prints: <name>: XXXX intervals
}

/************************ end of class Cover ************************/

/************************ class ReadDens ************************/

// Adds read to the container.
// Abstract BaseItems<> method implementation.
//	@rgn: Region with mandatory fields
//	@spotter: temporary values & ambiguities
//	return: true if read was added successfully
bool ReadDens::AddItem(const Region& rgn, Spotter& spotter)
{
	CheckStrand(spotter);
	_map[spotter.File().ItemStrand() ? rgn.Start : rgn.End]++;
	return true;
}

// Fills items from intermediate container
//	BaseItems<> abstract method implementation.
void ReadDens::FillChromItems(const Spotter& spotter)
{
	if (!_map.size())	return;
	chrlen prevEnd = 0;
	for (auto it = _map.begin(); it != _map.end(); it++) {
		AddPos(it->first, prevEnd, it->second);
		prevEnd = it->first + 1;
	}
	_items.push_back(ValPos(prevEnd, 0));			// add zero region after last unzero one
	_items.push_back(ValPos(spotter.chrLen - 1, 0));// add chrom boundary to finish iterating in CalcR correctly
	_map.clear();
}

/************************ end of class ReadDens ************************/

/************************ class JointedBeds ************************/

// Keeps mean values for fs1, fs2. 
//	@clear: if true, clear instance for treatment of new chromosome
void JointedBeds::R::Init(double mean1, double mean2, bool clear)
{
	double	d1 = 1 - mean1,
		d2 = 1 - mean2;

	_sqMeans1[0] = _sqMeans1[2] = mean1 * mean1;
	_sqMeans1[1] = _sqMeans1[3] = d1 * d1;
	_sqMeans2[0] = _sqMeans2[1] = mean2 * mean2;
	_sqMeans2[2] = _sqMeans2[3] = d2 * d2;
	_crossMeans[0] = mean1 * mean2;
	_crossMeans[1] = -mean2 * d1;
	_crossMeans[2] = -mean1 * d2;
	_crossMeans[3] = d1 * d2;
	if (clear)
		_cov = _var1 = _var2 = 0;
}

// Accumulates next length of range
void JointedBeds::R::Increment(chrlen len, char val)
{
	_cov += len * _crossMeans[val];
	_var1 += len * _sqMeans1[val];
	_var2 += len * _sqMeans2[val];
}

// Fills ChromRanges & Range by given two beds.
// Beds chromosomes should be checked as Treated.
// Both fs1 & fs2 must be valid: no duplicated, crossed, adjacent, coverage features;
// in other case R may be wrong
JointedBeds::JointedBeds(Features& fs1, Features& fs2)
{
	const Region fEnd = Region(CHRLEN_UNDEF, CHRLEN_UNDEF - 1);	// last chromosome's joint feature
	const char VAL1 = 0x1;	// value represented first Features's feature
	const char VAL2 = 0x2;	// value represented second Features's feature
	
	char val;							// current joint range value
	chrlen	pos, pos2;					// current positions and in fs2. Initially are equal
	chrlen	fi1, fi2;					// features indexes in fs1, fs2
	chrlen	firstInd = 0, lastInd = 0;	// current first, last feature indexes in JointedBeds
	Region	f1, f2;						// dedicated feature used for detecting
	BaseItems::cIter cit1, cit2;
	const BaseItems::cIter cit1end = fs1.End(), cit2end = fs2.End();

	_ranges.reserve(2 * (fs1.Count() + fs2.Count()));	// ranges

	for (cit1 = fs1.Begin(); cit1 != cit1end; cit1++) {
		if (!fs1.IsTreated(cit1))	continue;
		fi1 = fi2 = val = 0;
		cit2 = fs2.GetIter(CID(cit1));
		if(cit2 == cit2end)			continue;		// no chrom
		const chrlen fCnt1 = fs1.ItemsCount(cit1);	// count of features in fs1, fs2
		const chrlen fCnt2 = fs2.ItemsCount(cit2);
		f1 = fs1.Feature(cit1);
		f2 = fs2.Feature(cit2);
		// loop through current chromosome's features 
		while (fi1 < fCnt1 || fi2 < fCnt2) {
			pos = val & VAL1 ? (f1.End + 1) : f1.Start;
			pos2 = val & VAL2 ? (f2.End + 1) : f2.Start;
			if (pos < pos2) {
				val ^= VAL1;		// flip val for fs1
				if (!(val & VAL1))	// true when fs1 feature is closed (every second range)
					f1 = ++fi1 < fCnt1 ? fs1.Feature(cit1, fi1) : fEnd;
			}
			else if (pos > pos2) {
				pos = pos2;
				val ^= VAL2;		// flip val for fs2
				if (!(val & VAL2))	// true when fs2 feature is closed (every second range)
					f2 = ++fi2 < fCnt2 ? fs2.Feature(cit2, fi2) : fEnd;
			}
			else {
				val ^= VAL1 ^ VAL2;	// flip val for both beds
				if (!(val & VAL1))	// true when fs1 feature is closed 
					f1 = ++fi1 < fCnt1 ? fs1.Feature(cit1, fi1) : fEnd;
				if (!(val & VAL2))	// true when fs2 feature is closed 
					f2 = ++fi2 < fCnt2 ? fs2.Feature(cit2, fi2) : fEnd;
			}
			_ranges.push_back(Range(pos, val));	// add new joint feature
			lastInd++;
		}
		AddVal(CID(cit1), ChromRanges(
			firstInd,
			lastInd,
			fs1.FeaturesLength(cit1),
			fs2.FeaturesLength(cit2)
		));
		firstInd = lastInd;
	}
}

// Calculates r and fills results
//	@cSizes: chrom sizes
//	@results: object to fill results
void JointedBeds::CalcR(const ChromSizes& cSizes)
{
	const bool isPrLocal = PrintMngr::IsPrintLocal();
	const bool isPrTotal = PrintMngr::IsPrintTotal();
	const genlen gSize = isPrTotal ? cSizes.GenSize() : 0;	// genome's size

	chrlen	cSize;			// chromosome's size
	chrlen	start, stop;	// range's boundary positions
	chrlen	ri, len;		// index of range, length of range
	char	val;			// value of current range
	double	featrsLen1, featrsLen2;
	PairR	pairLocR, pairTotR;
	ChromRanges cRanges;	// current ChromRanges

	for (cIter it = cBegin(); it != cEnd(); it++) {
		cRanges = it->second.Data;
		cSize = cSizes[CID(it)];
		featrsLen1 = double(cRanges.FeatrsLen1);
		featrsLen2 = double(cRanges.FeatrsLen2);

		if (isPrLocal)	pairLocR.Init(featrsLen1 / cSize, featrsLen2 / cSize, true);
		if (isPrTotal)	pairTotR.Init(featrsLen1 / gSize, featrsLen2 / gSize, false);

		start = val = 0;	// first range
		for (ri = cRanges.FirstInd; ri <= cRanges.LastInd; ri++) {
			stop = _ranges[ri].Start;
			len = stop - start;
			if (isPrLocal)	pairLocR.Increment(len, val);
			if (isPrTotal)	pairTotR.Increment(len, val);
			// next range
			val = _ranges[ri].Val;
			start = stop;
		}
		len = cSize - stop;		// last range

		//== print current result
		if (isPrLocal) {
			pairLocR.Increment(len, 0);
			PrintMngr::PrintChrom(CID(it), pairLocR.Get());
		}
		if (isPrTotal)
			pairTotR.Increment(len, 0);
	}
	if (isPrTotal)
		PrintMngr::PrintTotal(pairTotR.Get());
}

#ifdef _DEBUG
void	JointedBeds::Print()
{
	chrlen	ri;				// index of range
	cout << "JointedBeds:\n";
	for (cIter it = cBegin(); it != cEnd(); it++) {
		cout << Chrom::AbbrName(CID(it)) << ":";
		for (ri = Data(it).FirstInd; ri <= Data(it).LastInd; ri++)
			cout << '\t' << _ranges[ri].Start << '\t' << int(_ranges[ri].Val) << LF;
	}
}
#endif
/************************ end of class JointedBeds ************************/

/************************ class CorrPair ************************/

static inline void DelWig(Obj* obj) { if (obj) delete (Cover*)obj; }
static inline void DelBedF(Obj* obj) { if (obj) delete (Features*)obj; }
static inline void DelBedR(Obj* obj) { if (obj) delete (ReadDens*)obj; }

CorrPair::FileType CorrPair::_FileTypes[] = {
	{ &CorrPair::CreateWig,  DelWig},
	{ &CorrPair::CreateBedF, DelBedF},
	{ &CorrPair::CreateBedR, DelBedR}
};

int CorrPair::_FileTypesCnt = sizeof(CorrPair::_FileTypes) / sizeof(CorrPair::FileType);

//void IgnoreOption(int opt, chrlen len)
//{
//	const string optName = Options::OptionToStr(opt);
//	Err("in order for the " + optName +
//		" to be valid, the space should be at least 5 times less than the shortest extended feature length of "
//		+ NSTR(len)).
//		Warning(". " + optName + " is ingored.");
//	Options::ResetIntVal(opt);
//}

// Creates an instance with checking primary object.
//	@cID: chromosome's ID
//	@primefName: primary file's name
//	@rgns: genome regions
//	@templ: name of template bed file, or NULL if undefined
//	@multiFiles: true if more then one secondary files are placed
CorrPair::CorrPair(const char* primefName, DefRegions& rgns, const char* templ, bool multiFiles) :
	_firstObj(NULL),
	_secondObj(NULL),
	_templ(NULL),
	_gRgns(rgns),
	_typeInd(CheckFileExt(primefName, true)),
	_printWarn(Options::GetBVal(oWARN)),
	_printName(multiFiles)
{
	if (templ)
		if (IsBedF()) {
			if (PrintMngr::IsNotLac())
				Err("ignored", string(sTemplate) + sBLANK + templ).Throw(false);
		}
		else {
			_templ = new Features(sTemplate, FS::CheckedFileName(templ),
				rgns.ChrSizes(), PrintMngr::Verb(), _printName, true, _printWarn);
			chrlen minDistance = _templ->GetMinDistance();
			chrlen extLen = Options::GetIVal(oEXTLEN);
			if (extLen > minDistance / 2) {
				dout << "WARNING: extended length of " << extLen
					<< " exceeds half the distance between the nearest features. Reduced to ";
				extLen = minDistance / 20;		// len /= 10, len *= 10; to round up to 10
				dout << (extLen *= 10) << LF;
			}
			_templ->Extend(extLen, rgns.ChrSizes(), PrintMngr::Verb());
		}
	if (PrintMngr::IsNotLac()) {
		if (CC::IsP())		dout << "Pearson";
		if (CC::IsBoth())	dout << AND << SPACE;
		if (CC::IsS())		dout << "Signal";
		dout << " CC between\n";
	}
	_firstObj = (this->*_FileTypes[_typeInd].Create)(primefName, true);
	_gRgns.Init();
	if (_printName || PrintMngr::IsNotLac()) {
		dout << AND;
		if (multiFiles)	dout << "...";
		dout << LF;
	}
}

CorrPair::~CorrPair() {
	if (_templ)		delete _templ;
	_FileTypes[_typeInd].Delete(_firstObj);
	_FileTypes[_typeInd].Delete(_secondObj);
}

// Adds secondary object, calculates and prints CCkey.
void CorrPair::CalcCC(const char* fName)
{
	_FileTypes[_typeInd].Delete(_secondObj);
	_secondObj = NULL;

	//== check type
	BYTE typeInd = CheckFileExt(fName, false);
	if (typeInd == _FileTypesCnt)	return;
	if (typeInd != _typeInd)
		return Err("different" + sFormat, fName).Throw(false);

	//== create object
	_secondObj = (this->*_FileTypes[_typeInd].Create)(fName, false);
	if (_secondObj->IsBad())	return;

	//== calculate r
	if (IsBedF()) {
		int expStep = Options::GetIVal(oEXTSTEP);
		if (expStep) {	// calculation r by step increasing expanding length
			int expLen = Options::GetIVal(oEXTLEN);
			for (int i = 0; i <= expLen; i += expStep) {
				Features corrBedF(*((Features*)_firstObj));

				corrBedF.Extend(i, _gRgns.ChrSizes(), PrintMngr::Verb());
				CalcCCBedF(corrBedF);
			}
			return;
		}
		else			// primary fs is already extended
			CalcCCBedF(*((Features*)_firstObj));
	}
	else				
		CalcCCWigCover(_gRgns);
}

// Creates features bed object.
//	@fName: file name
//	@primary: if true object is primary
Obj* CorrPair::CreateBedF(const char* fName, bool primary)
{
	Features* fs = new Features(NULL, fName, _gRgns.ChrSizes(),
		PrintMngr::Verb(), _printName, primary, _printWarn);
	// if primary is bad, it throws exception and quit, so check at once
	if (!fs->IsBad() && fs->IsStrandPres()) {
		if (fs->EOLNeeded())	dout << LF;
		Err("looks like an alignment!", fName).Warning();
	}
	if (primary)
		// expand primary fs f.e. to get true correlation
		// between narrow TFBS (primary) and wide recovered peaks (secondary)
		fs->Extend(Options::GetIVal(oEXTLEN), _gRgns.ChrSizes(), PrintMngr::Verb());
	return fs;
}

// Checks file extisting and extention validity
//	@fName: file's name
//	@abortInvalid: if true throw extention if checking is false
//	return: index in _FileTypes or _FileTypesCnt if checking is false
BYTE CorrPair::CheckFileExt(const char* fName, bool abortInvalid)
{
	if (FS::CheckFileExist(fName, abortInvalid))	return _FileTypesCnt;
	switch (FT::GetType(fName, Options::GetBVal(oALIGN))) {
	case FT::eType::BGRAPH:	return 0;
	case FT::eType::BED:	return 1;
	case FT::eType::ABED:
	case FT::eType::BAM:	return 2;
	}
	Err("unpredictable" + sFormat, fName).Throw(abortInvalid);
	return _FileTypesCnt;
}

/************************ end of class CorrPair ************************/
