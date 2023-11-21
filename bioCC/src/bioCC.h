/**********************************************************
bioCC.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 13.12.2021
-------------------------
Provides option emum and FileList class
***********************************************************/
#pragma once

enum optValue {
	oALIGN,
	oGENOM,
	oCHROM,
	oGAPLEN,
	oDUPL,
	oOVERL,
	oFILE_LIST,
	oFBED,
	oEXTLEN,
	oEXTSTEP,
	oPRCC,
	oBINWIDTH,
	oPRFCC,
	oVERB,
	oWRITE,
	oOUTFILE,
	oTIME,
	oVERSION,
	oSUMM,
	oHELP
};


// 'FileList' represents file's names incoming from argument list or from input list-file.
class FileList
/*
 * Under Windows should be translated with Character Set as not Unicode
 * (only 'Use Multi-Byte Character Set' or 'Not Set').
 */
{
private:
	char **_files = nullptr;	// file names
	short _count = 0;			// count of file names

public:
	// Constructor for list from input file.
	// Lines begining with '#" are commetns and would be skipped.
	FileList(const char* fileName);

	~FileList();
	
	// Gets count of file's names.
	inline short Count() const { return _count; }
	
	inline char** Files() const { return _files; }
	
	inline const char* operator[](int i) const { return _files[i]; }

#ifdef _DEBUG
	void Print() const {
		if( _files )
			for(short i=0; i<_count; i++)
				cout << _files[i] << LF;
		else
			cout << "Empty\n";
	}
#endif
};
