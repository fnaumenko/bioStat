/************************************************************************************
biostat: biostatistical package for NGS data
This is a command shell for calling statistical programs.

Copyright (C) 2019 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 09.01.2022
-------------------------
This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the	GNU General Public License for more details.
************************************************************************************/

#include <iostream>	
#include <iomanip>      // setw
#ifdef _WIN32
#include <windows.h>
//#include <cstring>
#include <algorithm>
#else
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>	// struct stat
#include <unistd.h>		// getcwd() & realink
#include <limits.h>		// PATH_MAX
//typedef unsigned char BYTE;
#endif

typedef unsigned char BYTE;

#define HPH		'-'
#define SPACE	' '
#define EOL		'\n'	// LF, 10

using namespace std;

using pairCommand = pair<const char*, const char*>;

constexpr pairCommand commands[] = {
	{ "cc",			"bioCC" },
	{ "calldist",	"callDist" },
	//{ "readdens",	"readDens" },
	{ "valign",		"vAlign" },
	{ "fqstatn",	"fqStatN" },
};
constexpr BYTE maxCommLen = sizeof(commands[1]);
constexpr BYTE commCnt = sizeof(commands) / sizeof(commands[0]);

const string appName = "biostat";
const char* optSumm = "--summ";

int PrintUsage(bool prTitle);
int CallApp(BYTE ind, const char* argv[] = NULL, int paramsCnt = 3);

// Incapsulates command line to launch utility
class CommLine
{
	char* _comm;	// C-string contained utility name and parameters

public:
	// Constructor
	//	@app: utility name
	//	@argv: program arguments
	//	@argc: number of program arguments
	CommLine(const char* app, const char* argv[], int argc) {
		size_t	appLen = strlen(app);
		size_t	shift = 0, parsLen;
		size_t	commLen = appLen + 1;	// command line length: +1 for closing 0

		// calculate commLen
		for (int i = 2; i < argc; commLen += strlen(argv[i++]) + 1);	// +1 for blank

		// in Debug mode the utilities are launched from their original (development) folders
#ifdef _DEBUG
		const char* pathHead = "D:\\Documents\\source\\reposCPP\\Release\\";
		const char* pathTail = "\\Release\\";
		commLen += strlen(pathHead) + strlen(pathTail) + appLen;
#endif // _DEBUG
		_comm = new char[commLen + 1];	// +1 for 0, 2 for ".\"
		memset(_comm, SPACE, commLen);
#ifdef _DEBUG
		memcpy(_comm, pathHead, shift = strlen(pathHead));
		memcpy(_comm + shift, app, appLen);
		memcpy(_comm + (shift += appLen), pathTail, strlen(pathTail));
		memcpy(_comm + (shift += strlen(pathTail)), app, appLen);
		shift += appLen;
#else
		memcpy(_comm, app, shift = appLen);
#endif // _DEBUG

		for (int i = 2; i < argc; i++) {
			memcpy(_comm + ++shift, argv[i], parsLen = strlen(argv[i]));
			shift += parsLen;
		}
		_comm[commLen] = '\0';		// end string
	}

	~CommLine() { delete[] _comm; }

	const char* Get() const { return _comm; }

	size_t Size() const { return strlen(_comm); }
};

#ifndef _WIN32
// Returns true if file exists in called folder
bool IsFileExist(const char* fname)
{
	// find real folder from which the main app was called
	char result[PATH_MAX];
	ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
	string path = string(result, (count > 0) ? count : 0);
	// replace main app name by util name
	path = path.substr(0, path.rfind(appName)) + string(fname);
	// check if util is exists in this folder
	struct stat st;
	return (!stat(path.c_str(), &st) && st.st_mode & S_IFREG);
}
#endif

int main(int argc, char* argv[])
{
	if (argc < 2)	return PrintUsage(true);
	if (*argv[1] == HPH) {
		char secLit = *(argv[1] + 1);
		if (secLit == 'h' || (secLit == HPH && !strcmp((const char*)(argv[1] + 2), "help")))
			return PrintUsage(true);
		cout << "wrong option; use -h|--help for help\n";
		return 1;
	}
	// define command
	char ind = -1;		// index
	for (char i = 0; i < commCnt; i++)
		if (!strcmp(argv[1], commands[i].first)) { ind = i; break; }
	if (ind == -1) {
		cout << "unrecognized command: " << argv[1] << EOL;
		return PrintUsage(false);
	}
	return CallApp(ind, (const char**)argv, argc);
}

int PrintUsage(bool prTitle)
{
	const char* ver = "1.0";
	if (prTitle)
		cout << appName.c_str() << ": statistical tools for NGS data\nVersion: "
		<< ver << "\n\nUsage:\t"
		<< appName.c_str() << " <command> [options]\n";
	cout << "\nCommands:\n";
	for (BYTE i = 0; i < commCnt; i++) {
		cout << setw(10) << commands[i].first << setw(2) << SPACE;
		CallApp(i);
	}
	cout << EOL;
	return !prTitle;
}

// Calls package
//	@ind: index of package in commands[]
//	@argv: program arguments
//	@argc: number of program arguments
//	return: exit code
int CallApp(BYTE ind, const char* argv[], int argc)
{
	const char* missUtil = ":\t this utility is missing\n";
	const char* app = commands[ind].second;
	const char* helpParams[] = { NULL, NULL, optSumm };
	int ret = 0;

	CommLine cm(app, argv ? argv : helpParams, argc);
	//cout << cm.Get() << EOL;	return 0;	// control output
#ifdef _WIN32
	STARTUPINFO si;
	PROCESS_INFORMATION pi;

	ZeroMemory(&si, sizeof(si));		// set the size of the structures
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	// start the program up
	WCHAR target[maxCommLen];
	MultiByteToWideChar(CP_ACP, 0, cm.Get(), -1, target, maxCommLen);
	if (CreateProcess(NULL, LPWSTR(target), NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi)) {
		WaitForSingleObject(pi.hProcess, INFINITE);	// Wait until child process exits.
		// Close process and thread handles. 
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);
	}
	else {
		ret = GetLastError();
		if (ret == 2)	cerr << app << missUtil;
		else			cerr << "error " << ret << EOL;
	}
#else
	if (IsFileExist(app)) {
		FILE* fp = popen(cm.Get(), "r");
		if (fp) {
			char buff[256];
			while (fgets(buff, sizeof(buff), fp)) cout << buff;
			pclose(fp);
		}
		else { cerr << app << " open error\n";	ret = 1; }
	}
	else { cerr << app << missUtil; ret = 2; }
#endif
	return ret;
}
