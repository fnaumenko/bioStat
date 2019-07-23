/**********************************************************
bioCC.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 08.06.2019
-------------------------
Provides option emum and FileList class
***********************************************************/
#pragma once

enum optValue {
	oALIGN,
	oGEN,
	oGAPLEN,
	oDUPL,
	//oDIFFSZ,
	oFILE_LIST,
	oCHROM,
	oCC,
	oSPACE,
	oFBED,
	oEXTLEN,
	oEXTSTEP,
	oPRCC,
	oBINWIDTH,
	oFRES,
	oFNORM,
	//oCROSS,
	//oCIRC,
	//oSTEP,
	oINFO,
	oWARN,
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
	char **_files;		// file names
	short _count;		// count of file names
	bool _memRelease;	// true if memory should be free in destructor

public:
	// Constructor for argument list.
	FileList	(char* files[], short cntFiles);
	
	// Constructor for list from input file.
	// Lines begining with '#" are commetns and would be skipped.
	FileList	(const char* fileName);

	~FileList()	{
		if( _files && _memRelease ) {
			for(short i=0; i<_count; i++)
				free(_files[i]);
			delete [] _files;
		}
	}
	
	// Gets count of file's names.
	inline short Count() const { return _count; }
	
	inline char** Files() const { return _files; }
	
	inline const char* operator[](int i) const { return _files[i]; }

#ifdef DEBUG
	void Print() const {
		if( _files )
			for(short i=0; i<_count; i++)
				cout << _files[i] << EOL;
		else
			cout << "Empty\n";
	}
#endif
};

