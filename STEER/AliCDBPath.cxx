/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBPath						   //
//  Path string identifying the object:  			   //
//  "level0/level1/level2" 					   //
//  (was: "Detector/DBType/DetSpecType") 		           //
//  (example: "ZDC/Calib/Pedestals") 		         	   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBPath.h"

#include <TObjArray.h>
#include <TObjString.h>
#include <TRegexp.h>

#include "AliLog.h"

ClassImp(AliCDBPath)

//_____________________________________________________________________________
AliCDBPath::AliCDBPath() :
  TObject(),
  fPath(""),
  fLevel0(""),
  fLevel1(""),
  fLevel2(""),
  fIsValid(kTRUE),
  fIsWildcard(kFALSE)
{
// default constructor

}

//_____________________________________________________________________________
AliCDBPath::AliCDBPath(const AliCDBPath& other):
  TObject(other),
  fPath(other.fPath),
  fLevel0(""),
  fLevel1(""),
  fLevel2(""),
  fIsValid(other.fIsValid),
  fIsWildcard(other.fIsWildcard)
{
// constructor
	Init();
	InitPath();

}

//_____________________________________________________________________________
AliCDBPath::AliCDBPath(const char* level0, const char* level1,
	const char* level2):
  TObject(),
  fPath(""),
  fLevel0(level0), 
  fLevel1(level1), 
  fLevel2(level2),
  fIsValid(kTRUE),
  fIsWildcard(kFALSE)
{
// constructor

	fPath += level0;
	fPath += '/';
	fPath += level1;
	fPath += '/';
	fPath += level2;

	if ((IsWord(fLevel0) || fLevel0 == "*")
		&& (IsWord(fLevel1) || fLevel1 == "*")
		&& (IsWord(fLevel2) || fLevel2 == "*")) {

		fIsValid = kTRUE;
	} else {
		fIsValid = kFALSE;
		AliError(Form("Invalid AliCDBPath <%s/%s/%s>!", 
			level0, level1, level2));
	}

	Init();
}

//_____________________________________________________________________________
AliCDBPath::AliCDBPath(const char* path):
  TObject(),
  fPath(path),
  fLevel0(""),
  fLevel1(""),
  fLevel2(""),
  fIsValid(kTRUE),
  fIsWildcard(kFALSE)
{
// constructor

	Init();
	InitPath();	
}

//_____________________________________________________________________________
AliCDBPath::AliCDBPath(const TString& path):
  TObject(),
  fPath(path),
  fLevel0(""),
  fLevel1(""),
  fLevel2(""),
  fIsValid(kTRUE),
  fIsWildcard(kFALSE)
{
	Init();
	InitPath();
}

//_____________________________________________________________________________
void AliCDBPath::InitPath() {
// sets fLevel0, fLevel1, fLevel2, validity flagss from fPath

	TSubString strippedString = fPath.Strip(TString::kBoth);
	TString aString(strippedString);
	strippedString = aString.Strip(TString::kBoth, '/');

	TObjArray* anArray = TString(strippedString).Tokenize("/");
	Int_t paramCount = anArray->GetEntriesFast();
	
	if (paramCount == 1) {
		if (fPath == "*") {
			fLevel0 = "*";
			fLevel1 = "*";
			fLevel2 = "*";
			
			fIsValid = kTRUE;
		} else {
			fIsValid = kFALSE;
		}

	} else if (paramCount == 2) {
		fLevel0 = ((TObjString*) anArray->At(0))->GetString();
		TString aString =  ((TObjString*) anArray->At(1))->GetString();

		if (IsWord(fLevel0) && aString == "*") {
			fLevel1 = "*";
			fLevel2 = "*";
		
			fIsValid = kTRUE;			

		} else {
			fIsValid = kFALSE;
		}
		
	} else if (paramCount == 3) {
		fLevel0 = ((TObjString*) anArray->At(0))->GetString();
                fLevel1 = ((TObjString*) anArray->At(1))->GetString();
                fLevel2 = ((TObjString*) anArray->At(2))->GetString();

		if ((IsWord(fLevel0) || fLevel0 == "*")
	                && (IsWord(fLevel1) || fLevel1 == "*")
        	        && (IsWord(fLevel2) || fLevel2 == "*")) {
			
			fIsValid = kTRUE;
		} else {
			fIsValid = kFALSE;
		}

	} else {
		fIsValid = kFALSE;

	}
	
	if (!fIsValid) {
		AliInfo(Form("Invalid AliCDBPath <%s>!", fPath.Data()));
	} else {	
		fPath = Form("%s/%s/%s", fLevel0.Data(), fLevel1.Data(), fLevel2.Data());
	}
	
	delete anArray;
	
	Init();
}

//_____________________________________________________________________________
AliCDBPath::~AliCDBPath() {
// destructor

}

//_____________________________________________________________________________
Bool_t AliCDBPath::IsWord(const TString& str) {
// check if string is a word

	TRegexp pattern("^[a-zA-Z0-9_.]+$");

	return str.Contains(pattern);	
}

//_____________________________________________________________________________
void AliCDBPath::Init() {
// set fIsWildcard flag

	fIsWildcard = fPath.MaybeWildcard();	
}

//_____________________________________________________________________________
Bool_t AliCDBPath::Level0Comprises(const TString& str) const {
// check if Level0 is wildcard or is equal to str
	
	if (fLevel0 == "*") {
		return kTRUE;
	}

	return fLevel0 == str;
}

//_____________________________________________________________________________
Bool_t AliCDBPath::Level1Comprises(const TString& str) const {
// check if Level1 is wildcard or is equal to str

	if (fLevel1 == "*") {
		return kTRUE;
	}

	return fLevel1 == str;
}

//_____________________________________________________________________________
Bool_t AliCDBPath::Level2Comprises(const TString& str) const {
// check if Level2 is wildcard or is equal to str
	
	if (fLevel2 == "*") {
		return kTRUE;
	}

	return fLevel2 == str;
}

//_____________________________________________________________________________
Bool_t AliCDBPath::Comprises(const AliCDBPath& other) const {
// check if path is wildcard and comprises other

	return Level0Comprises(other.fLevel0) 
		&& Level1Comprises(other.fLevel1)
		&& Level2Comprises(other.fLevel2);
}
