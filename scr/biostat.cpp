/************************************************************************************
biostat: biostatistical package for NGS data
This is a command shell for calling statistical programs.
	
Copyright (C) 2019 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 21.07.2019
-------------------------
This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the	GNU General Public License for more details.
************************************************************************************/

#include <iostream>	
#include <iomanip>      // std::setw
#ifdef _WIN32
	#include <windows.h>
#else
	#include <string.h>
	#include <stdio.h>
	#include <sstream>
	#include <sys/stat.h>	// struct stat
	#include <unistd.h>		// getcwd() & realink
	#include <limits.h>		// PATH_MAX
	typedef unsigned char BYTE;
#endif

#define HPH		'-'
#define BLANK	' '
#define EOL		'\n'	// LF, 10
#define TAB		'\t'

using namespace std;

struct pairCommand { const char* first; const char* second; };

const pairCommand commands[] = {
	{ "cc",			"bioCC" },
	{ "fragdist",	"fragDist" },
	{ "readdens",	"readDens" },
	{ "valign",		"vAlign" },
	{ "fqstatn",	"fqStatN" },
};
const BYTE	commCnt = sizeof(commands)/sizeof(pairCommand);

const string appName = "biostat";
const char* optSumm = "--summ";

int PrintUsage(bool prTitle);
void CallApp(BYTE ind, const char* params, size_t paramsLen);


int main(int argc, char* argv[])
{
	if(argc < 2)	return PrintUsage(true);
	if(*argv[1] == HPH) {
		if(argv[1][1]=='h' || (*(argv[1]+1)==HPH && !strcmp((const char*)(argv[1]+2), "help")))
			return PrintUsage(true);
		cout << "wrong option; use -h|--help for help\n";
		return 1;
	}
	// define command
	char ind = -1;		// index
	for(char i=0; i<commCnt; i++)
		if(!strcmp(argv[1], commands[i].first))	{ ind = i; break; }
	if(ind == -1) {
		cout << "unrecognized command: " << argv[1] << EOL;
		return PrintUsage(false);
	}
	
	size_t len = 0;
	for(int i=2; i<argc; i++)	len += strlen(argv[i]) + 1;
	if(len)	CallApp(ind, argv[2], len);
	else	CallApp(ind, optSumm, strlen(optSumm));

	return 0;
}

int PrintUsage(bool prTitle)
{
	const char* ver = "1.0";
	if(prTitle)
		cout << appName.c_str() << ": statistical tools for NGS data\nVersion: "
			 << ver << "\n\nUsage:\t"
			 << appName.c_str() << " <command> [options]\n";
	cout << "\nCommands:\n";
	for(BYTE i=0; i<commCnt; i++) {
		cout << setw(10) << commands[i].first << setw(2) << BLANK;
		//cout << commands[i].second << EOL;
		CallApp(i, optSumm, strlen(optSumm));
	}
	cout << EOL;
	return !prTitle;
}

#ifndef _WIN32
// Returns true if file exists in called folder
bool IsFileExist(const char* fname)
{
	// find real folder from which the main app was called
	char result[PATH_MAX];
	ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
	string path = string(result, (count > 0) ? count : 0);
	// replace main app name by util name
	path = path.substr(0, path.rfind(appName)) + string(fname);
	// check if util is exists in this folder
	struct stat st;
	return (!stat(path.c_str(), &st) && st.st_mode & S_IFREG);
}
#endif

//void CallApp(const char* app, const char* params, size_t paramsLen)
void CallApp(BYTE ind, const char* params, size_t paramsLen)
{
	const char* missUtil = " is missing\n";
	const char* app = commands[ind].second;

#ifndef _WIN32
	if(!IsFileExist(app)) {	cerr << app << missUtil; return; }
#endif

	size_t appLen = strlen(app) + 1;	// strlen with closing 0!
	char* comm = new char[appLen + paramsLen];
	memcpy(comm, app, appLen);			// as we know appLen it's faster than strcpy
	memcpy(comm + appLen, params, paramsLen);
	{	// replace 0 with blank
		size_t i = --appLen;
		appLen += paramsLen;
		for(; i<appLen; i++)
			if(!comm[i])	comm[i] = BLANK;
	}
	memset(comm + ++appLen, '\0', 1);		// set final closing 0

#ifdef _WIN32
	STARTUPINFO si;     
	PROCESS_INFORMATION pi;

	ZeroMemory( &si, sizeof(si) );		// set the size of the structures
	si.cb = sizeof(si);
	ZeroMemory( &pi, sizeof(pi) );
	// start the program up
	if( CreateProcess(NULL, (LPSTR)comm, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi) )
	{
		// Wait until child process exits.
		WaitForSingleObject( pi.hProcess, INFINITE );
		
		// Close process and thread handles. 
		CloseHandle( pi.hProcess );
		CloseHandle( pi.hThread );
	}
	else {
		unsigned int err = GetLastError();
		if(err == 2)	cerr << app << missUtil;
		else			cerr << "error " << err << EOL;
	}
#else
	FILE *fp = popen(comm, "r");
	if(fp) {
		char buff[256];
		while(fgets(buff, sizeof(buff), fp)) cout << buff;
		pclose(fp);
	}
	else 	cerr << app << " open error\n";
#endif

	delete [] comm;
}
