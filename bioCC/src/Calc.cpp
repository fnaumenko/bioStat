/**********************************************************
Calc.ccp
Provides classes for calculating CC
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 11.22.2023
***********************************************************/

#include "Calc.h"
#include "ChromData.h"
#include <algorithm>    // std::sort

const string sFormat = " format";
const char* sUNDEF = "UNDEF";
const int Undef = -3;	// undefined coefficient

// Prints PCC
void PrintR(float cc) 
{
	dout << left << setw(12) << setfill(SPACE);
	if(cc == Undef || isNaN(cc))	dout << sUNDEF;
	else	dout << cc;
	dout << LF;
}

// Stores the operation state of the - executed or not
class R
{
	mutable bool _done = false;
protected:
	// Returns PCC
	//	@cov: covariance
	//	@var1: variance 1
	//	@var2: variance 2
	float GetR(const double& cov, const double& var1, const double& var2) const {
		_done = true;
		return cov ? (var1 && var2 ? float(cov / sqrt(var1) / sqrt(var2)) : Undef) : 1;
	}

public:
	// Returns true if calculation was actually done
	bool IsDone() const { return _done; }
};

/********************  PrintMngr *********************/

PrintMngr::ePrint	PrintMngr::_PrintCC;
eOInfo				PrintMngr::_OInfo;
bool				PrintMngr::_PrName;

// Initialisez by user setings
void PrintMngr::Init(int ccPrint, eOInfo info, bool multiFiles)
{
	_PrintCC = ePrint (ccPrint);
	_OInfo = info;
	_PrName = IsNotLac() || multiFiles;
}

// Throws an exception for an empty object
void PrintMngr::CompleteEmpty(const char* fName, FT::eType type)
{
	ostringstream ss;
	if (!_PrName)	ss << fName;
	if (_OInfo <= eOInfo::NM)	ss << SepCl;
	ss << "no " << FT::ItemTitle(type);
	if (Chrom::UserCID() != Chrom::UnID)
		ss << " for stated " << Chrom::ShortName(Chrom::UserCID());
	Err(ss.str()).Throw();
}

// Prints PCC
//	@cc: correlation coefficient
//	@cID: calculated chrom pcc or total
void PrintMngr::PrintCC(float cc, chrid cID)
{
	if (cID == Chrom::UnID) {	// total CC
		if (_OInfo > eOInfo::LAC || IsPrintTotal())
			dout << "total:\t";
	}
	else						// chrom CC
		dout << Chrom::AbbrName(cID) << TAB;

	PrintR(cc);
	fflush(stdout);		// when called from a package 
}

/********************  end of PrintMngr *********************/

/************************ class PlainCover ************************/

bool PlainCover::AddPos(const ValPos& vPos, chrlen prevEnd)
{
	bool add;
	// is it a gap between prev and current regions?
	if (add = (vPos.Val && vPos.Pos != prevEnd))					// check for the unzero val to avoid duplicate zero items
		_items.emplace_back(prevEnd);	// add zero space after previous unzero region
	if (add = (!_items.size() || _items.back().Val != vPos.Val))	// check prev val for wiggle_0 with regular tiems with the same span (MACS output)
		_items.push_back(vPos);				// add start of current unzero region
	return add;
}

void PlainCover::AddChrom(chrlen cID, chrlen cLen, chrlen prevEnd)
{
	if (!prevEnd)	return;			// skip first call with cID == Chrom::UnID
	// add 2 last items
	_items.emplace_back(prevEnd);	// add zero region after last unzero one
	if(prevEnd < cLen)
		_items.emplace_back(cLen);	// add chrom boundary to finish iterating in CalcR correctly

	// add chrom
	const chrlen lastInd = chrlen(_items.size());
	AddVal(cID, ItemIndices(chrlen(_lastInd), lastInd));	// minus added 1 first and 2 last zero intervals
										// !!! replace AddVal by emplace
	_lastInd = lastInd;
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
//	@gRgns: def regions (chrom sizes)
//	@templ: template to define treated regions
//	return: true if calculation was actually done
bool PlainCover::CalcR(const PlainCover& cv, const DefRegions& rgns, const Features* templ)
{
	// 'spR' - single-pass Pearson coefficient (R) calculater; keeps accumulates sums & calculates PCC
	class spR : public R
	{
		double	_sumX = 0, _sumY = 0;		// sum of signal values
		double	_sumXY = 0;					// sum of the products of the values of both signals
		double	_sumSqrX = 0, _sumSqrY = 0;	// sum of squared signal values
		genlen	_len = 0;

	public:
		void Clear() { _sumX = _sumY = _sumXY = _sumSqrX = _sumSqrY = 0; _len = 0; }

		// Adds range length and correlated range values
		void AddVal(chrlen len, float valX, float valY) {
			const double valXLen = double(valX) * len;
			const double valYLen = double(valY) * len;

			_sumX += valXLen;	_sumSqrX += valX * valXLen;
			_sumY += valYLen;	_sumSqrY += valY * valYLen;
			_sumXY += valXLen * valY;
			_len += len;
		}

		// Returnes Pearson CC
		float PCC() {
			return GetR(_len * _sumXY - _sumX * _sumY,
				_len * _sumSqrX - _sumX * _sumX,
				_len * _sumSqrY - _sumY * _sumY);
		}
	};

	// 'FeatureR' represetns pair <feature-ID><feature-PCC>
	struct FeatureR : pair<chrlen, float>
	{
		inline FeatureR(chrlen id, float cc) { first = id; second = cc; }

		inline bool operator < (const FeatureR& rccr) const { return second < rccr.second; }

		void Print() const { dout << first << TAB; PrintR(second); }
	};

	// 'FeatureRs' represetns FeatureR collection, including methods to print collection and CC histogram
	class FeatureRs : vector<FeatureR>
	{
		// Creates and prints histogram
		void PrintHist(float binWidth)
		{
			// ** set abs values and sort
			for (auto& i : *this)
				if (i.second < 0)	i.second = -i.second;
			sort(begin(), end());	// by increase
			
			// ** define factor - a divisor of binWidth: 0.1--0.9=>10, 0.01--0.09=>100 etc
			short F = 10;
			for (; binWidth * F < 1; F *= 10);

			// ** define min scaled bin value
			auto it = begin();
			// then float instead of double because of wrong consolidation by round double
			float minBin = float(int(it->second * F)) / F;
			
			// ** define max scaled bin value
			auto itEnd = prev(end());	// pointed to the LAST item!
			size_t undefCnt = 0;		// count of regions with undefined CC
			while (itEnd->second > 1)	// eliminate UNDEF CC
				itEnd--, undefCnt++;
			float maxBin = float(F * itEnd->second);
			{	// round maxBin
				int	maxdecBin = int(maxBin);
				if (maxBin - maxdecBin)	maxdecBin++;	// round up
				if (maxdecBin % 2)		maxdecBin++;	// get even bin
				maxBin = float(maxdecBin) / F;
			}
			vector<int> hist(size_t((maxBin - minBin) / binWidth) + 1, 0);		// histogram

			// ** fill histogram by consolidated values
			while (it <= itEnd)
				hist[int((maxBin - it++->second) / binWidth)]++;

			// ** cut off low bins with zero value
			size_t lim = hist.size() - 1;
			for (size_t k = lim; k; k--)
				if (!hist[k]) lim--;
				else break;
			// ** print histogram
			dout << "BIN UP\tCOUNT\n";
			for (BYTE k = 0; k <= lim; k++)
				dout << (maxBin - k * binWidth) << TAB << hist[k] << LF;
			if (undefCnt)
				dout << sUNDEF << TAB << undefCnt << LF;
		}

	public:
		inline FeatureRs(chrlen cnt) { reserve(cnt); }

		inline void AddVal(chrlen ind, float val) { emplace_back(ind + 1, val); }

		// Prints result and histogram
		void Print(eRS printFRes, float binWidth)
		{
			if (printFRes != rsOFF) {
				if (printFRes == rsC)		// soretd by feature; are sorted initially
					sort(begin(), end());	// by increase
				dout << "\n#RGN\tCC\n";
				for (const auto& cc : *this)		cc.Print();
			}
			if (binWidth)	PrintHist(binWidth);
		}
	};

	const bool fillLocRes = templ && (_binWidth || _printFRes);
	spR chrR, totR;
#ifdef _DEBUG
	//Print("first", 0);
	//cv.Print("second", 0);	cout << LF;
	////rgns.Print();
	//if (templ)	templ->Print();
#endif
	// rgns is already limited by chroms represented in template, if it's defined
	for (auto rit = rgns.cBegin(); rit != rgns.cEnd(); rit++) {
		chrid cID = CID(rit);
		if (!FindChrom(cID) || !cv.FindChrom(cID))	continue;
		chrlen fCnt = 0;							// count of templ features
		Items<Featr>::cItemsIter itF, itFend;		// template feature iterator
		if (templ) {
			auto itC = templ->GetIter(cID);
			itF = templ->ItemsBegin(itC);
			itFend = templ->ItemsEnd(itC);
			fCnt = chrlen(templ->ItemsCount(itC));
		}

		// local results
		spR locR;
		chrlen ind = 0;				// item index
		FeatureRs locResults(_binWidth ? fCnt : 0);	// create histogram
		chrlen len;					// length of current united (with equal X, Y values) region
		chrlen posN, pos = 0;		// max(posX,posY), current position, 
		float valX = 0, valY = 0;	// first sequence, second sequence current value
		bool insideF = !templ;		// true if current position is inside current template feature
		bool inTempl = templ;		// true if current position did not go beyond the border of the last feature 
		bool closeF = false;		// true if the feature has just ended
		const auto itXend = ItemsEnd(cID), itYend = cv.ItemsEnd(cID);

		chrR.Clear();
		// loop through cover items (intervals)
		for (auto itX = ItemsBegin(cID), itY = cv.ItemsBegin(cID);
			itX != itXend && itY != itYend; ) {
			const chrlen posX = itX->Pos, posY = itY->Pos;
			const float prevValX = valX, prevValY = valY;	// X, Y current value

			//== set valX, valY
			if (posX > posY)			// more
				posN = posY, valY = itY++->Val;
			else {						// equal or less
				posN = posX; valX = itX++->Val;
				if (posX == posY) 		// equal
					valY = itY++->Val;
			}
			//== set pos, len
			if (inTempl)
				if (posN > itF->End) {						// exit current feature
					closeF = true;
					len = SetPos(pos, itF->End);			// take len including the end of the feature
					inTempl = ++itF != itFend;				// next feature
				}
				else if (!insideF && posN >= itF->Start)	// entry current feature
					insideF = true,
					len = (pos = posN) - itF->Start;		// cutoff range before feature start
				else 
					len = SetPos(pos, posN);
			else 
				len = SetPos(pos, posN);

			//== accumulate sums; the last interval, ended by pos, is saved!!
			if (insideF && len) {		// len can be 0 at the boundary of the feature
				if (fillLocRes)		
					locR.AddVal(len, prevValX, prevValY);
				chrR.AddVal(len, prevValX, prevValY);		// previous combined region
				if (PrintMngr::IsPrintTotal())
					totR.AddVal(len, prevValX, prevValY);
					//totR.AddVal(len, valX, valY);
			}

			//== close feature, save loc CC
			if (closeF) {
				if (fillLocRes) {			// add CC for each feature
					locResults.AddVal(ind++, locR.PCC());
					locR.Clear();
				}
				closeF = insideF = false;
			}
		}
		
		//== print current result
		if (PrintMngr::IsPrintLocal()) {
			if (templ)
				locResults.Print(_printFRes, _binWidth);
			PrintMngr::PrintCC(chrR.PCC(), cID);
		}
	}
	if (PrintMngr::IsPrintTotal())
		PrintMngr::PrintCC(totR.PCC());
	return chrR.IsDone() || totR.IsDone();
}

// Writes inner representation to BEDGRAPG file
void PlainCover::Write(const string& fName) const
{
	ofstream file;
	file.open(fName);
	file << "track type=bedGraph\n";
	for (const auto& c : Container()) {
		const string& chr = Chrom::AbbrName(c.first);
		const auto itEnd = --ItemsEnd(c.second.Data);
		for (auto it = ItemsBegin(c.second.Data); it != itEnd; it++)
			if(it->Val)
				file << chr << TAB << it->Pos << TAB << next(it)->Pos << TAB << it->Val << LF;
	}
	file.close();
}

/************************ end of PlainCover ************************/

/************************ class Cover ************************/

// Initializes instance from wig file
//	return: numbers of all and initialied items for given chrom
void Cover::InitWiggle(BedReader& file, const ChromSizes& cSizes)
{
	static const string keyChrom = "chrom";
	static const string keyStart = "start";
	static const string keyStep = "step";
	static const string keySpan = "span";

	bool	skipChrom = false;		// if true skip data lines for current chrom
	chrlen	prevEnd = 0,			// previous region end position, 
			span = 1, step = 0;		// span spec, step spec (for fixedStep),
	const char* line;
	const bool fixedStep = file.Type() == FT::eType::WIG_FIX;
	chrid	cID = Chrom::UnID, nextCID = cID;	// current, next chrom ID
	ULONG	cItemCnt = 0,	// count of total accepted intervals
			itemCnt = 0,	// count of accepted intervals of current chrom, total
			recCnt = 0;		// count of total records
	BYTE	offset = 0;		// offset to value delimiter (TAB or SPACE); for Variable Step, only increases
	ValPos	vPos;			// region start position & value
	char	firstC;
	Timer timer(UniBedReader::IsTimer);
	// pointer to initialise position & value function
	void (*setValPos)(ValPos & vPos, chrlen step, BYTE & offset, const char* s);

	if(fixedStep) {
		firstC = 'f';
		setValPos = [](ValPos& vPos, chrlen step, BYTE&, const char* s) {
			vPos.Pos += step;
			vPos.Val = float(atof(s));
		};
	}
	else {
		firstC = 'v';
		setValPos = [](ValPos& vPos, chrlen, BYTE& offset, const char* s) {
			vPos.Pos = atoi(s);
			for (s += offset; *s != SPACE && *s != TAB; ++s, ++offset);	// go to value
			vPos.Val = float(atof(++s));
		};
	}

	while (line = file.GetNextLine(false))
		if (*line == firstC) {		// *** declaration line
			// *** set specifications
			// ** check keyChrom
			line = file.ChromMark() - strlen(Chrom::Abbr);					// level the initial value of chrom mark position
			const char* line1 = line = file.CheckSpec(line, keyChrom) + 1;	// keyChrom initial position

			if (fixedStep) {							// fix declarative parameters
				// * check keyStart & keyStep
				vPos.Pos = file.GetIntKey(line1, keyStart);	// initial position
				line1 += keyStart.length() + 1;			// shift to scan the rest of the line faster
				step = file.GetIntKey(line1, keyStep);	// initial step
				vPos.Pos -= step;							// shift 'back' before the first pass, where pos+=step will be invoke
				line1 += keyStep.length() + 1;			// shift to scan the rest of the line faster
			}
			line1 = TabReader::KeyStr(line1, keySpan);
			span = line1 ? atoi(line1 + 1) : 1;			// initial span: both for fixed- and variableStep
			// * check chrom
			if (file.GetNextChrom(nextCID, line)) {
				if (Chrom::IsSetByUser()) {				// are all chroms specified?
					if (cID != Chrom::UnID)				// skip first pass, when curr chrom is still undefined
						PlainCover::AddChrom(cID, cSizes[cID], prevEnd),
						itemCnt += cItemCnt;
				}
				else {								// single chrom is specified
					if (!fixedStep && vPos.Pos)		// region is initialized: items for the specified chrom are existed and saved
						break;						// the chrom itself will be saved after loop
					if (skipChrom = nextCID != Chrom::UserCID())
						continue;
				}
				cID = nextCID;
				cItemCnt = prevEnd = offset = 0;
				vPos.Clear();
			}
		}
		else {				// *** data line
			if (skipChrom)	continue;
			setValPos(vPos, step, offset, file.GetLine());
			cItemCnt += AddPos(vPos, prevEnd);
			prevEnd = vPos.Pos + span;
			recCnt++;
		}
	// save last chrom
	if (cID != Chrom::UnID) {	// is last chrom valid?
		PlainCover::AddChrom(cID, cSizes[cID], prevEnd);
		itemCnt += cItemCnt;
	}
	// print stats
	if (PrintMngr::OutInfo() >= eOInfo::STD) {
		if (itemCnt == 1)	itemCnt = 0;	// single interval is equal to an empty coverage
		UniBedReader::PrintItemCount(itemCnt, FT::ItemTitle(FT::eType::WIG_FIX, itemCnt != 1));
		if (PrintMngr::OutInfo() == eOInfo::STAT)
			dout << " (" << recCnt << " data lines)";
		if (!(Timer::Enabled && UniBedReader::IsTimer))	dout << LF;
	}
	timer.Stop(1, true, PrintMngr::OutInfo() > eOInfo::NM);
}

// Creates new instance by wig-file name
// Invalid instance wil be completed by throwing exception.
//	@fName: file name
//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
//	@prfName: true if file name should be printed unconditionally, otherwise deneds on oinfo
//	@abortInval: true if invalid instance should abort excecution
Cover::Cover(const char* fName, ChromSizes& cSizes, eOInfo oinfo, bool prfName, bool abortInval)
	: PlainCover()
{
	UniBedReader file(fName, FT::eType::BGRAPH, &cSizes, 4, 0, oinfo, prfName, true, abortInval);

	ReserveItems(file.EstItemCount());	// EstItemCount() > 0 even for empty file, because of track line
	if (file.Type() == FT::eType::BGRAPH)
		Pass(this, file);
	else
		InitWiggle((BedReader&)file.BaseFile(), cSizes);
	
	if (Options::GetBVal(oWRITE)) {
		const string ext = FS::GetExt(fName);
		Write(FS::FileNameWithoutExt(fName) + "_out." + ext);
	}

	//PrintEst(file.EstItemCount());
}

/************************ end of class Cover ************************/

/************************ class ReadDens ************************/

// Adds chrom to the instance
//	@cID: current chrom ID
//	@cLen: current chrom length
void ReadDens::AddChrom(chrid cID, chrlen cLen)
{
	// fill items
	chrlen prevEnd = 0;
	for (const freqPair& i : *_freq) {
		AddPos(i, prevEnd);
		prevEnd = i.first + 1;
	}
	_freq->clear();
	PlainCover::AddChrom(cID, cLen, prevEnd);
}

// Creates new instance by abed/bam-file name
// Invalid instance wil be completed by throwing exception.
//	@fName: file name
//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
//	@printfName: true if file name should be printed unconditionally, otherwise deneds on oinfo
//	@abortInval: true if invalid instance should abort excecution
ReadDens::ReadDens(const char* fName, ChromSizes& cSizes, eOInfo oinfo, bool printfName, bool abortInval)
	: PlainCover()
{
	RBedReader file(fName, &cSizes, BYTE(Options::GetRDuplLevel(oDUPL)), oinfo, printfName, abortInval);
	rfreq freq;
	_freq = &freq;

	ReserveItems(file.EstItemCount());
	Pass(this, file);
	_freq = nullptr;

	//PrintEst(file.EstItemCount());
}

/************************ end of class ReadDens ************************/

// Fills ChromRanges & Range by given two beds.
// Beds chromosomes should be checked as Treated.
// Both fs1 & fs2 must be valid: no duplicated, crossed, adjacent, coverage features;
// in other case R may be wrong
JointedBeds::JointedBeds(const Features& fs1, const Features& fs2)
{
	const Region fEnd = Region(CHRLEN_UNDEF, CHRLEN_UNDEF - 1);	// last chromosome's joint feature
	const char VAL1 = 0x1;	// value represented first Features's feature
	const char VAL2 = 0x2;	// value represented second Features's feature
	const auto cit1end = fs1.cEnd(), cit2end = fs2.cEnd();

	chrlen	firstInd = 0, lastInd = 0;	// current first, last feature indexes in JointedBeds
	Region	r1, r2;						// dedicated feature used for detecting

	_ranges.reserve(2 * (fs1.Count() + fs2.Count()));	// ranges
	for (auto cit1 = fs1.cBegin(); cit1 != cit1end; cit1++) {
		if (!fs1.IsTreated(cit1))	continue;
		auto cit2 = fs2.GetIter(CID(cit1));
		if(cit2 == cit2end)			continue;		// no chrom
		const chrlen fCnt1 = chrlen(fs1.ItemsCount(cit1));	// count of features in fs1, fs2
		const chrlen fCnt2 = chrlen(fs2.ItemsCount(cit2));
		r1 = fs1.Regn(cit1);
		r2 = fs2.Regn(cit2);
		char val = 0;							// current joint range value
		// loop through current chromosome's features 
		for (chrlen fi1 = 0, fi2 = 0; fi1 < fCnt1 || fi2 < fCnt2;) {
			chrlen pos = val & VAL1 ? (r1.End + 1) : r1.Start;
			const chrlen pos2 = val & VAL2 ? (r2.End + 1) : r2.Start;
			if (pos < pos2) {
				val ^= VAL1;		// flip val for fs1
				if (!(val & VAL1))	// true when fs1 feature is closed (every second range)
					r1 = ++fi1 < fCnt1 ? fs1.Feature(cit1, fi1) : fEnd;
			}
			else if (pos > pos2) {
				pos = pos2;
				val ^= VAL2;		// flip val for fs2
				if (!(val & VAL2))	// true when fs2 feature is closed (every second range)
					r2 = ++fi2 < fCnt2 ? fs2.Feature(cit2, fi2) : fEnd;
			}
			else {
				val ^= VAL1 ^ VAL2;	// flip val for both beds
				if (!(val & VAL1))	// true when fs1 feature is closed 
					r1 = ++fi1 < fCnt1 ? fs1.Feature(cit1, fi1) : fEnd;
				if (!(val & VAL2))	// true when fs2 feature is closed 
					r2 = ++fi2 < fCnt2 ? fs2.Feature(cit2, fi2) : fEnd;
			}
			_ranges.emplace_back(pos, val);	// add new joint feature
			lastInd++;
		}
		AddVal(CID(cit1), ChromRanges(
			firstInd, lastInd,
			fs1.FeaturesLength(cit1),
			fs2.FeaturesLength(cit2)
		));
		firstInd = lastInd;
	}
}

// Calculates r and fills results
//	@cSizes: chrom sizes
//	@results: object to fill results
//	return: true if calculation was actually done
bool JointedBeds::CalcR(const ChromSizes& cSizes)
{
	// 'dsR' - discrete signal Pearson coefficient (R) calculater; keeps cumulative data & calculates PCC
	class dsR : public R	// DsR
	{
		double	_var1 = 0, _var2 = 0;	// variances 
		double	_cov = 0;				// covariance: SUM( X-Xmean)^2 * (Y-Ymean)^2 )
		/*
		Mean value for every range may accept only one of two values: 0-mean or 1-mean,
		and they are const within chromosome.
		So for efficiency we can keep theirs in 3 arrays:
		two square's arrays: arrays of squared values from each bed;
		there are only 2 combinations, but they are duplicated for the simplicity acces,
		and one crossing array: array of all combinations
		of multiplications different values from fs1 and fs2
		*/
		double _sqMeans1[4], _sqMeans2[4];	// square's arrays for fs1, fs2
		double _crossMeans[4];				// crossing array

	public:
		// Keeps mean values for both features. 
		//	@clear: if true, clear instance for treatment of new chromosome
		void Init(const double(&mean)[2], bool clear) {
			const double d1 = 1 - mean[0], d2 = 1 - mean[1];

			_sqMeans1[0] = _sqMeans1[2] = mean[0] * mean[0];
			_sqMeans1[1] = _sqMeans1[3] = d1 * d1;
			_sqMeans2[0] = _sqMeans2[1] = mean[1] * mean[1];
			_sqMeans2[2] = _sqMeans2[3] = d2 * d2;
			_crossMeans[0] = mean[0] * mean[1];
			_crossMeans[1] = -mean[1] * d1;
			_crossMeans[2] = -mean[0] * d2;
			_crossMeans[3] = d1 * d2;
			if (clear)
				_cov = _var1 = _var2 = 0;
		}

		// Accumulates next length of range
		void Increment(chrlen len, char val) {
			_cov += len * _crossMeans[val];
			_var1 += len * _sqMeans1[val];
			_var2 += len * _sqMeans2[val];
		}

		// Returns Pearson coefficient of correlation
		inline float PCC() const { return GetR(_cov, _var1, _var2); }
	};

	const bool isPrLocal = PrintMngr::IsPrintLocal();
	const bool isPrTotal = PrintMngr::IsPrintTotal();
	const genlen gSize = isPrTotal ? cSizes.GenSize() : 0;	// genome's size
	dsR chrR, totR;

	for (const auto& c : Container()) {						// loop through chroms
		const ChromRanges& cRanges = c.second.Data;			// current ChromRanges
		const chrlen cSize = cSizes[c.first];
		const double fsLen[2]{ 
			double(cRanges.FeatrsLen1) / cSize,
			double(cRanges.FeatrsLen2) / cSize
		};
		if (isPrLocal)	chrR.Init(fsLen, true);
		if (isPrTotal)	totR.Init(fsLen, false);

		chrlen start = 0, stop = 0, len;	// range's boundary positions, length of range
		char val = 0;						// value of current range
		for (auto ri = cRanges.FirstInd; ri <= cRanges.LastInd; ri++) {
			stop = _ranges[ri].Start;
			len = stop - start;
			if (isPrLocal)	chrR.Increment(len, val);
			if (isPrTotal)	totR.Increment(len, val);
			// next range
			val = _ranges[ri].Val;
			start = stop;
		}
		len = cSize - stop;		// last range

		//== print current result
		if (isPrLocal) {
			chrR.Increment(len, 0);
			PrintMngr::PrintCC(chrR.PCC(), c.first);
		}
		if (isPrTotal)
			totR.Increment(len, 0);
	}
	if (isPrTotal)
		PrintMngr::PrintCC(totR.PCC());

	return chrR.IsDone() || totR.IsDone();
}

#ifdef _DEBUG
void	JointedBeds::Print()
{
	chrlen	ri;				// index of range
	cout << "JointedBeds:\n";
	for (cIter it = cBegin(); it != cEnd(); it++) {
		cout << Chrom::AbbrName(CID(it)) << COLON;
		for (ri = Data(it).FirstInd; ri <= Data(it).LastInd; ri++)
			cout << TAB << _ranges[ri].Start << TAB << int(_ranges[ri].Val) << LF;
	}
}
#endif
/************************ end of class JointedBeds ************************/

/************************ class CorrPair ************************/

CorrPair::FileType CorrPair::_FileTypes[] = {
	{ &CorrPair::CreateWig,  [](void* obj) { delete (Cover*)obj; }},
	{ &CorrPair::CreateBedF, [](void* obj) { delete (Features*)obj; }},
	{ &CorrPair::CreateBedR, [](void* obj) { delete (ReadDens*)obj; }}
};

// Creates an instance with checking primary object.
//	@cID: chromosome's ID
//	@primefName: primary file's name
//	@rgns: genome regions
//	@tfName: name of template bed file, or NULL if undefined
//	@multiFiles: true if more then one secondary files are placed
CorrPair::CorrPair(const char* primefName, DefRegions& rgns, const char* tfName, bool multiFiles) :
	_gRgns(rgns),
	_typeInd(CheckFileExt(primefName, true))
{
	PrintMngr::Init(Options::GetIVal(oPRCC), eOInfo(Options::GetIVal(oVERB)), multiFiles);
	UniBedReader::IsTimer = PrintMngr::OutInfo() > eOInfo::LAC;
	if (tfName)
		if (IsBedF()) {
			if (PrintMngr::IsNotLac()) {
				ostringstream ss;
				ss << sTemplate << SPACE << tfName;
				Err("ignored", ss.str()).Throw(false);
			}
		}
		else {
			if (PrintMngr::IsPrName())	dout << sTemplate << SepCl;
			_templ = new Features(FS::CheckedFileName(tfName), _gRgns.ChrSizes(),
				Options::GetBVal(oOVERL), PrintMngr::OutInfo(), PrintMngr::IsPrName(), true);
			CheckItemsCount(_templ, tfName);

			chrlen extLen = Options::GetIVal(oEXTLEN);
			if (extLen) {
				chrlen minDistance = _templ->GetMinDistance();
				if (extLen > minDistance / 2) {
					extLen = minDistance / 2 - 1;		// len /= 10, len *= 10; to round up to 10
					ostringstream ss;
					ss	<< "extended length exceeds half the distance between the nearest features. Reduced to "
						<< extLen;
					Err(ss.str(), PrintMngr::EchoName(tfName)).Warning();
				}
				_templ->Extend(extLen, rgns.ChrSizes(), UniBedReader::eAction::ABORT);
			}
		}
	if (PrintMngr::IsNotLac()) 	dout << "Pearson CC between\n";

	_firstObj = (this->*_FileTypes[_typeInd].Create)(primefName, true);
	_gRgns.Init();
	if (PrintMngr::IsNotLac()) {
		dout << " and";
		if (multiFiles)	dout << "...";
		dout << LF;
	}
}

CorrPair::~CorrPair() {
	delete _templ;
	_FileTypes[_typeInd].Delete(_firstObj);
	_FileTypes[_typeInd].Delete(_secondObj);
}

// Adds secondary object, calculates and prints CCkey.
void CorrPair::CalcCC(const char* fName)
{
	{	//== check type
		char typeInd = CheckFileExt(fName, false);
		if (typeInd != _typeInd) {
			if (typeInd == vUNDEF)		return;		// already checked
			return Err("different" + sFormat, fName).Throw(false);
		}
	}
	_FileTypes[_typeInd].Delete(_secondObj);
	_secondObj = nullptr;

	//== create object
	try { _secondObj = (this->*_FileTypes[_typeInd].Create)(fName, false); }
	catch (const Err& e) { dout << e.what() << LF; return; }
	
	//== calculate r
	bool done;
	if (IsBedF()) {
		if (done = CalcCCBedF(*((Features*)_firstObj))) {		// 'zero extended'
			const int extStep = Options::GetIVal(oEXTSTEP);
			if (extStep) {		// calculation r by step increasing expanding length
				const int extLen = Options::GetIVal(oEXTLEN);
				if (extLen < extStep)
					Err("extending length is less then extending step. Extension has stopped here.",
						PrintMngr::EchoName(fName)).Warning();
				else {
					Features bedF(*((Features*)_firstObj));
					for (int i = extStep; i <= extLen; i += extStep) {
						dout << "primer extended by " << i << ":\n";
						if (!bedF.Extend(extStep, _gRgns.ChrSizes(), UniBedReader::eAction::ABORT))
							break;
						CalcCCBedF(bedF);
					}
				}
			}
		}
	}
	else				
		done = CalcCCCover(_gRgns);
	if (!done)
		Err("no " + FT::ItemTitle(_type) + " for common " + Chrom::Title(true)).
			Throw(false, true);
}

// Creates features bed object.
//	@fName: file name
//	@primary: if true object is primary
void* CorrPair::CreateBedF(const char* fName, bool primary)
{
	Features* obj = new Features(fName, _gRgns.ChrSizes(), Options::GetBVal(oOVERL),
		PrintMngr::OutInfo(), PrintMngr::IsPrName(), primary);
	CheckItemsCount(obj, fName);
	if (obj->NarrowLenDistr())
		Err("looks like an alignment but is handled as ordinary bed!", PrintMngr::EchoName(fName)).Warning();
	return obj;
}

// Creates alignment object
void* CorrPair::CreateBedR(const char* fName, bool isPrimary)
{
	ReadDens* obj = new ReadDens(fName, _gRgns.ChrSizes(),
		PrintMngr::OutInfo(), PrintMngr::IsPrName(), isPrimary);
	CheckItemsCount(obj, fName);
	return obj;
}

// Creates covering object.
//	@fName: file name
//	@primary: if true object is primary
void* CorrPair::CreateWig(const char* fName, bool isPrimary)
{
	Cover* obj = new Cover(fName, _gRgns.ChrSizes(),
		PrintMngr::OutInfo(), PrintMngr::IsPrName(), isPrimary);
	CheckItemsCount(obj, fName);
	return obj;
}

// Checks file extisting and extention validity
//	@fName: file's name
//	@abortInvalid: if true throw extention if checking is false
//	return: index in _FileTypes or _FileTypesCnt if checking is false
char CorrPair::CheckFileExt(const char* fName, bool abortInvalid)
{
	if (!FS::CheckFileExist(fName, abortInvalid)) {
		switch (_type = FT::GetType(fName, Options::GetBVal(oALIGN))) {
		case FT::eType::BGRAPH:	return 0;
		case FT::eType::BED:	return 1;
		case FT::eType::ABED:
		case FT::eType::BAM:	return 2;
		}
		Err("unpredictable" + sFormat, fName).Throw(abortInvalid);
	}
	return vUNDEF;
}

/************************ end of class CorrPair ************************/
