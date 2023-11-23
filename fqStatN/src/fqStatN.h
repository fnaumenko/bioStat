#pragma once
#include "TxtFile.h"

enum optValue {
	oOUTFILE,
	oTIME,
	oSUMM,
	oVERSION,
	oHELP
};

static class StatN
{
private:
	struct TemplN {
		ULONG Count;		// 'N' count 
		ULONG CountRead;	// Reads count
		vector<char> Pos;	// template with marked N positions; Array to ensure destructor
		//Array<char> Pos;	// template with marked N positions; Array to ensure destructor

		inline TemplN() : Count(0), CountRead(0) {}

		//inline TemplN(ULONG cnt, const Array<char>& pos)
		inline TemplN(ULONG cnt, const vector<char>& pos)
			: Count(cnt), CountRead(1) { Pos = pos;	}

		TemplN(const TemplN& t) { 
			Count = t.Count;
			CountRead = t.CountRead;
			Pos = t.Pos;
		}
	};

	static inline bool Compare (TemplN& el1, TemplN& el2)
	{ return el1.CountRead > el2.CountRead; }

public:
	static void Scan (FqReader& fqFile);

} statN;
