#pragma once
#include "FqReader.h"

enum optValue {
	oDOUT_FILE,
	oTIME,
	oSUMM,
	oVERSION,
	oHELP
};

static class StatN
{
private:
	struct TemplN
	{
		size_t Count = 0;		// 'N' count 
		size_t CountRead = 0;	// Reads count
		vector<char> Pos;	// template with marked N positions

		TemplN(size_t cnt, const vector<char>& pos) : Count(cnt), CountRead(1), Pos(pos) {}
	};

	static bool Compare (TemplN& el1, TemplN& el2) { return el1.CountRead > el2.CountRead; }

public:
	static void Scan (FqReader& fqFile);

} statN;
