/**********************************************************
bioCC.h
Provides option emum and FileList class
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 07/25/2024
***********************************************************/
#pragma once

enum optValue {
	oALIGN,
	oGENOM,
	oCHROM,
	oGAP_LEN,
	oDUPL,
	oOVERL,
	oFILE_LIST,
	oFBED,
	oEXT_LEN,
	oEXT_STEP,
	oPR_CC,
	oBIN_WIDTH,
	oPR_FCC,
	oVERB,
	oWRITE,
	oDOUT_FILE,
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

#ifdef MY_DEBUG
	void Print() const;
#endif
};
