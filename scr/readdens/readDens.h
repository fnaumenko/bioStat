/**********************************************************
readDens.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 08.06.2019
-------------------------
Provides option emum, ChromReadDistrib, GenReadDistrib classes
***********************************************************/

#pragma once
#include "Data.h"

enum {
	oGEN,
	oCHROM,
	oFBED,
	oEXTLEN,
	oGAPLEN,
	oDUPL,
	//oDIFFSZ,
	oMINSCR,
	oSPACE,
	oSERV,
	oINFO,
	oALARM,
	oFREQ,
	oOUTFILE,
	oTIME,
	oVERSION,
	oSUMM,
	oHELP
};

// 'GenReadDistrib' generates and keeps chrom Read's distributions
class ChromReadDistrib
{
public:
	// 'WinFreq' represents windows frequencies collection
	class WinFreq : public map<short, UINT>
	{
	public:
		typedef map<short, UINT>::iterator Iter;			// iterator
		typedef map<short, UINT>::const_iterator cIter;		// constant iterator

		// Sets count of empty windows (zero frequency)
		//	wCnt: count of windows
		void SetZeroFreq(chrlen wCnt)
		{
			if(!wCnt)	return;
			int cnt = 0;	// number of nonzero windows
			for(cIter it=begin(); it!=end(); cnt += it++->second);
			if(wCnt != cnt)		(*this)[0] = wCnt - cnt;
		}

		//// Gets total count of windows
		//chrlen WinCnt()
		//{
		//	chrlen cnt = 0;
		//	for(cIter it = begin(); it != end(); it++)
		//		cnt += it->second;
		//	return cnt;
		//}
#ifdef DEBUG
		void Print()
		{
			for(cIter it=begin(); it!=end(); it++)
				cout << it->first << TAB << it->second << EOL;
		}
#endif
	};

	// 'Pocket' represents temporary members needed while scanning
	class Pocket
	{
	private:
		// 'ScanVals' represents temporary values needed while scanning
		struct ScanVals
		{
		private:
			ULONG	_rCnt;		// accumulative count of Reads
			UINT	_freq;		// window's frequency
			chrlen	_rgnLen;	// accumulative length of all regions
			chrlen	_gapLen;	// accumulative length of all gaps
			chrlen	_lastEnd;	// the end of last Region
			chrlen	_wPrevInd;	// index of last window (or number of previous windows minus 1)
	
		public:
			inline ScanVals() { _rCnt = _freq = _rgnLen = _gapLen = _lastEnd = _wPrevInd = 0; }

			// Increases gap length by given Region
			//	@rgn: region
			void IncrGapLen(const Region& rgn)
			{
				_gapLen += rgn.Start - _lastEnd;
				_lastEnd = rgn.End;
			}

			// Increases region length by given Region
			//	@rgn: region
			inline void IncrRgnLen(const Region& rgn) {	_rgnLen += rgn.Length(); }


			// Adds Read to the distribution
			//	@centre: Read's centre
			//	@wFreq: windows frequency distribution
			void AddRead(chrlen centre, WinFreq& wFreq)
			{
				chrlen wInd = (centre - _gapLen) / WinLen;	// curr Read window index
				if(wInd > _wPrevInd) {				// Read belongs to next window
					if(_freq)						// save freq for previous windows
						wFreq[_freq]++,				// increment number of windows with given freq
						_freq = 0;
					_wPrevInd = wInd;
				}
				_freq++;
				_rCnt++;
			}

			// Gets count of windows	
			inline chrlen WinCnt() const { return _rgnLen / WinLen + int(bool(_rgnLen % WinLen)); }

			// Gets Gets count Reads	
			inline ULONG ReadCnt() const { return _rCnt; }

			// Gets window's frequency
			inline UINT Freq()	const { return _freq; }
			
			//	Gets density
			inline float Dens()	const { return ReadDens(_rCnt, _rgnLen); }
		};

		ScanVals _vals[Gr::Cnt];	// FG/BG scan values
		Reads::cItemsIter _rit;		// current Reads iterator
		Reads::cItemsIter _ritEnd;	// iterator referring to the past-the-end Read for current chrom

	public:
		Pocket(const Reads &bedR, Reads::cIter cit) :
			_rit(bedR.ReadsBegin(cit)), _ritEnd(bedR.ReadsEnd(cit)) {}
		
		// Increases defined gap length
		//	@rgn: defined region
		inline void IncrDefGapLen(const Region& rgn) {	_vals[Gr::BG].IncrGapLen(rgn); }

		// Scans region: increases number of Reads, total length, frequency and other scan values
		//	@g: ground
		//	@rgn: scanned region
		//	@wFreq: windows frequencies collection
		void ScanRegion(Gr::Type g, const Region& rgn, WinFreq* wFreq)
		{
			chrlen centre;	// Read's centre

			_vals[g].IncrGapLen(rgn);
			_vals[g].IncrRgnLen(rgn);
			for(; _rit!=_ritEnd; _rit++) {
				centre = _rit->Centre();
				if(centre >= rgn.End)	return;
				_vals[g].AddRead(centre, wFreq[g]);
			}
		}

		// Gets window's frequency
		inline UINT Freq(Gr::Type g)		const { return _vals[g].Freq(); }

		//	Gets density
		inline float Dens(Gr::Type g)		const { return _vals[g].Dens(); }

		// Gets count of all windows	
		inline chrlen WinCnt(Gr::Type g)	const { return _vals[g].WinCnt(); }

		// Gets count Reads	
		inline ULONG ReadCnt(Gr::Type g)	const { return _vals[g].ReadCnt() ; }
	};

	WinFreq	_wFreq[Gr::Cnt];	//	Gr::FG/BF windows frequencies collection
	ULONG	_rCnt[Gr::Cnt];		//	Gr::FG/BF count of Reads
	float	_dens[Gr::Cnt];		//	Gr::FG/BF density

public:
	static chrlen WinLen;
	static bool	IsFG;		// true if foreground distribution is defined
	static bool	PrintFreq;	// true if windows frequency should be printed

	// Scans chrom to define Read density and windows frequencies
	//	@bedR: scanning Reads in BED
	//	@cit: scanning chrom iterator
	//	@rgns: defined regions
	//	@templ: Gr::FG features in BED
	void Scan(const Reads &bedR, Reads::cIter cit, const Regions &rgns, const Features *templ)
	{
		chrlen prevEnd;
		Regions tRgns;				// template regions
		Regions::Iter tit, tit_end;	// template regions iterator, end iterator
		Pocket pocket(bedR, cit);

		if(templ) {
			templ->FillRegions(CID(cit), tRgns);
			tit = tRgns.Begin();
		}
		for(Regions::Iter it = rgns.Begin(); it != rgns.End(); it++) {
			pocket.IncrDefGapLen(*it);
			if(templ) {
				prevEnd = it->Start;
				tit_end = tRgns.ExtEnd(tit, it->End);
				for(; tit != tit_end; tit++) {
					pocket.ScanRegion(Gr::BG, Region(prevEnd, tit->Start), _wFreq);
					prevEnd = tit->End;
					pocket.ScanRegion(Gr::FG, *tit, _wFreq);
				}
				pocket.ScanRegion(Gr::BG, Region(prevEnd, it->End), _wFreq);
			}
			else
				pocket.ScanRegion(Gr::BG, *it, _wFreq);
		}
		if(templ)	ScanComplete(Gr::FG, pocket);
		ScanComplete(Gr::BG, pocket);
	}

	// Finalize scanning
	//	@g: ground
	//	@pocket: scanning values
	void ScanComplete(Gr::Type g, Pocket& pocket)
	{
		_wFreq[g][pocket.Freq(g)]++;				// add last frequency
		_wFreq[g].SetZeroFreq(pocket.WinCnt(g));	// add zero frequency
		_rCnt[g] = pocket.ReadCnt(g);
		_dens[g] = pocket.Dens(g);
	}

	// Gets total number of Reads in distribution
	inline ULONG ReadCount(Gr::Type g) const { return _rCnt[g]; }

	//	Gets density
	//inline float Dens(Gr::Type g) const { return _dens[g]; }

	void PrintDens	()
	{
		dout << setprecision(3);
		if(IsFG)	dout << _dens[Gr::FG] << TAB;
		dout << _dens[Gr::BG];
		if(IsFG && _dens[Gr::BG])
			dout << TAB << _dens[Gr::FG]/_dens[Gr::BG];
		//else	dout << TAB << (ReadCount(Gr::FG) + ReadCount(Gr::BG));
		else	dout << TAB << ReadCount(Gr::BG);
		dout << EOL;
		if(PrintFreq) {
			dout << "reads:\t";
			if(IsFG)	PrintTwoDistr();
			else		PrintDistr();
			dout << EOL;
		}
	}

	void PrintTwoDistr()
	{
		WinFreq::cIter it		= _wFreq[Gr::FG].begin();	// FG iterator
		WinFreq::cIter it_end	= _wFreq[Gr::FG].end();
		WinFreq::cIter BGit		= _wFreq[Gr::BG].begin();	// BG iterator
		WinFreq::cIter BGit_end	= _wFreq[Gr::BG].end();
		bool equal;

		dout << _rCnt[Gr::FG] << TAB << _rCnt[Gr::BG] << EOL;
		const char* wins = " wins";
		dout << "r/win\t" << Gr::Title(Gr::FG) << wins << TAB << Gr::Title(Gr::BG) << wins << EOL;
		PrintHorLine(2*8 + DigitsCount(_rCnt[Gr::BG]));

		while(true)
			if(it->first < BGit->first) {
				dout << it->first << TAB << it->second << EOL;
				if(++it == it_end)	{ it = BGit; it_end = BGit_end;	break; }
			}
			else {
				dout << BGit->first << TAB;
				if(equal = it->first == BGit->first)
					dout << it->second;
				dout << TAB << BGit->second << EOL;
				if(equal && ++it == it_end)	break;
				if(++BGit == BGit_end)		break;
			}
		// print rest of distr
		for(; it != it_end; it++) {
			dout << it->first << TAB;
			if(it_end == BGit_end)	dout << TAB;
			dout << it->second << EOL;
		}
	}

	void PrintDistr()
	{
		dout << _rCnt[Gr::BG] << EOL;
		for(WinFreq::cIter it = _wFreq[Gr::BG].begin(); it!=_wFreq[Gr::BG].end(); it++)
			dout << it->first << TAB << it->second << EOL;
	}
};

// 'GenReadDistrib' generates and keeps a collection of chrom Read's distributions
class GenReadDistrib : public Chroms<ChromReadDistrib>
{
public:
	// Initializes In-|Out-Read's distribution container:
	// for each chromosome initializes default Regions and PairReadDistrib.
	// Only chromosomes marked as 'Treated' would be treated.
	//	@bedR: Reads wich are distribeted
	//	@gRegn: define genome Regions
	//	@templ: features determining peaks or NULL
	GenReadDistrib(const Reads& bedR, DefRegions& gRegn, const Features* templ)
	{
		ChromReadDistrib::IsFG = templ ? true : false;
		Reserve(bedR.ChromCount());
		for(Reads::cIter it = bedR.cBegin(); it != bedR.cEnd(); it++)
			if(bedR.IsTreated(it))
				AddEmptyElem(CID(it)).Data.Scan(bedR, it, gRegn[CID(it)], templ);
	}

	
	// Returns total number of Reads
	//	@g: ground
	ULLONG	ReadCount(Gr::Type g) const
	{
		ULLONG res = 0;
		for(cIter it = cBegin(); it != cEnd(); it++)
			res += it->second.Data.ReadCount(g);
		return res;
	}

	void PrintDensity() 
	{
		//chrid	cID;
		//float	fgDens, bgDens;
		//ULONG	rCnt[GrCnt], wCnt[GrCnt];	// total count of reads in locations, total count of windows in locations
		//memset(rCnt, 0, GrCnt*sizeof(ULONG));
		//memset(wCnt, 0, GrCnt*sizeof(ULONG));

		// header
		const char* mDend = "Mean density, ";
		int w;
		dout << "Mean density, " << UnitDens << EOL;
		if(ChromReadDistrib::IsFG) {
			const char* ratio = "ratio";
			dout << TAB << Gr::Title(Gr::FG) << TAB << Gr::Title(Gr::BG) << TAB << ratio << EOL;
			w = 3*8 + strlen(ratio);
		}
		else
			w = strlen(mDend) + strlen(UnitDens);
		PrintHorLine(w);

		Iter itPenult = End()--;	// penultimate, last valid iterator
		for(Iter it=Begin(); it!=End(); it++) {
			dout << Chrom::AbbrName(CID(it)) << TAB;
			it->second.Data.PrintDens();
			//if(ChromCount() > 1 && it != itPenult)	dout << EOL;
		}
		//if( ChromCount() > 1 ) {
			dout << Total << SepClTab;
			if(ChromReadDistrib::IsFG)	
				dout << "FG reads " << ReadCount(Gr::FG) << "\tBG reads ";
			dout << ReadCount(Gr::BG) << EOL;
			//if(_isFG)	dout << (inDens = GetTotalDensity(0, rCnt, wCnt)) << TAB << TAB;
			//dout << (outDens = GetTotalDensity(1, rCnt, wCnt));
			//if(_isFG)	dout << TAB << TAB << inDens/outDens;
			//dout << EOL;
		//}
	}
};
