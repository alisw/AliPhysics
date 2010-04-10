// $Id$

/************************************************************************
**
**
** This file is property of and copyright by the Department of Physics
** Institute for Physic and Technology, University of Bergen,
** Bergen, Norway, 2006
** This file has been written by Sebastian Bablok,
** sebastian.bablok@ift.uib.no
**
** Important: This file is provided without any warranty, including
** fitness for any particular purpose.
**
**
*************************************************************************/

//  @file   AliHLTPendolinoListReader.cxx
//  @author Sepastian Bablok
//  @date   
//  @brief  Helper class for pendolino list handling
//  @note   maintained by matthias.richter@cern.ch

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include <TObjString.h>
#include <TList.h>

#include "AliHLTPendolinoListReader.h"

ClassImp(AliHLTPendolinoListReader)

using namespace std;
//using namespace alice::hlt::pendolino;

//typedef pair<string, TMap> Key_Val_Pair; 
//typedef TMap<TString, TMap> AliasMap;


int AliHLTPendolinoListReader::kMAX_LINE_LENGTH = 256;


AliHLTPendolinoListReader::AliHLTPendolinoListReader()
  : TObject()
  , fValid(false)
  , fCalibObjList()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}


AliHLTPendolinoListReader::AliHLTPendolinoListReader(const char* filename) 
  : TObject()
  , fValid(false)
  , fCalibObjList()
{
  // see header file for class documentation
	if (ReadListFromFile(filename)) {
		fValid = true;
	} else {
		// log that filename does ot matches a file
	}
}


AliHLTPendolinoListReader::~AliHLTPendolinoListReader() 
{
  // see header file for class documentation
	if (fValid) {
		fCalibObjList.Delete();
	}
}


bool AliHLTPendolinoListReader::ReadListFromFile(const std::string filename)
{
  // see header file for class documentation
	return ReadListFromFile(filename.c_str());
}


bool AliHLTPendolinoListReader::ReadListFromFile(const char* filename)
{
  // see header file for class documentation
	ifstream pfile;
	unsigned int count = 0;

	// open property file
	if (!filename || filename[0]==0) {
		cout << "Filename for calibration object list is empty" << endl;
//		ProxyLogger::getLogger()->createLogMessage(MSG_ERROR, LOG_SOURCE_PROXY,
//				"PropertyReader received empty filename to read properties.");
		return false;
	}

	pfile.open(filename, ios_base::in);
	if (!pfile) {
		cout << "Error while opening list file or file does not exist."
				<< endl;
//		ProxyLogger::getLogger()->createLogMessage(MSG_ERROR, LOG_SOURCE_PROXY,
//				"PropertyReader is unable to open property file '" +
//				ProxyLogger::catos(filename) + "'.");
		return false;
	}
//	ProxyLogger::getLogger()->createLogMessage(MSG_INFO, LOG_SOURCE_PROXY,
//			"PropertyReader uses property file '" +	ProxyLogger::catos(filename)
//			+ "'.");

	while ((!pfile.eof()) && (!pfile.bad())) {
		char line[kMAX_LINE_LENGTH];
		char* ptrDelimiter;

		pfile.getline(line, kMAX_LINE_LENGTH);

		// skip empty lines
		if (strlen(line) == 0) {
			continue;
		}
		// skip comment lines (beginning with '#') and lines with blank chars
		// at beginnig
		if ((line[0] == '#') || (line[0] == ' ')) {
			continue;
		}

		// cut off comments at end
		ptrDelimiter = strstr(line, "#");
		if (ptrDelimiter != 0) {
			*ptrDelimiter = '\0';
		}
		// cut off end of line
		ptrDelimiter = strstr(line, " ");
		if (ptrDelimiter != 0) {
			*ptrDelimiter = '\0';
		}
		
        // seperate alias names by detectors
		ptrDelimiter = strstr(line, "=");
		if (ptrDelimiter != 0) {
			*ptrDelimiter = '\0';  	// replace '=' with '\0'
			ptrDelimiter++;			// points to values (Aliasname)
			if ((ptrDelimiter != 0) && (ptrDelimiter[0] != ' ') && (line[0] != ' ')) {
//				cout <<"getting detector"<< endl;
//				TString detector = new TString(line);
//				string value(ptrDelimiter);
				TList* pos = (TList*) fCalibObjList.GetValue(line);
				if (pos == 0) {
					// new detector
//					cout<<"inserting"<<endl;
					TList *aliasList = new TList();
			 		aliasList->Add(new TObjString(ptrDelimiter));
//					cout<<"including"<<endl;
					fCalibObjList.Add(new TObjString(line), (TObject*) aliasList); 
//							new TMap(new TObjString(ptrDelimiter), 0));
//					cout << "inserted: " << line << " with " << ptrDelimiter << endl;
				} else {
					if (pos->FindObject(ptrDelimiter) == 0) {
//						cout<<"adding"<<endl;
//						mProperties.insert(String_Pair(name, value));
						pos->Add(new TObjString(ptrDelimiter));
						// use value [Alias name] as key for map
//						cout << "added " << ptrDelimiter << " for " << line << endl;
					} else {
						cout << "   ~~~ ommiting duplicated Alias names (" <<
								ptrDelimiter << ")." << endl;
						continue;
					}
				}
				++count;
			} else {
//				ProxyLogger::getLogger()->createLogMessage(MSG_WARNING,
//						LOG_SOURCE_PROXY,
//						"PropertyReader encountered empty property for '" +
//						ProxyLogger::catos(line) + "'.");
				cout << "Missing property value for " << line << endl;
			}
		}

/*
	// file calibration object name in local vector
	string name(line);
	fCalibObjList.push_back(name);
	++count;
*/
	}
//	ProxyLogger::getLogger()->createLogMessage(MSG_DEBUG, LOG_SOURCE_PROXY,
//			"PropertyReader has read " + ProxyLogger::itos(count) +
//			" properties from file '" +	ProxyLogger::catos(filename) + "'.");

	// maybe check also for .bad() for returning a false
	pfile.close();
	fValid = true;
	return true;
}

void AliHLTPendolinoListReader::Print() const {
	TIter iter(&fCalibObjList);
	TObjString* detect;

	cout << "List Print:" << endl;
	while ((detect = (TObjString*) (iter.Next()))) {
		cout << " " << detect->String() << ":" << endl;
		TMap* value = (TMap*) fCalibObjList.GetValue(detect->String().Data());
		TIter valIt(value);
		TObjString* alias;
		while ((alias = (TObjString*) (valIt.Next()))) {
			cout << "   " << alias->String() << endl;
		}
		cout << endl;
	}
	cout << "---" << endl;
}

int AliHLTPendolinoListReader::RetrieveLastRunNumber(char* path) {
	string strPath = path;
	return AliHLTPendolinoListReader::RetrieveLastRunNumber(strPath);
}

int AliHLTPendolinoListReader::RetrieveLastRunNumber(string path) {
	int runNumber = 0;
	ifstream pfile;
	string filename;
	char line[25];

	filename = path + "/lastRunNumber";
	cout << "   --- PATH to lastRunNumber is " << filename << endl;
	
	pfile.open(filename.c_str(), ios_base::in);
	if (!pfile) {
		cout << "Error while opening last RunNumber file or file does not exist."
				<< endl;
		return runNumber;
	}
	if ((!pfile.eof()) && (!pfile.bad())) {
		pfile.getline(line, 25);
		runNumber = atoi(line);
	}

	return runNumber;
}


