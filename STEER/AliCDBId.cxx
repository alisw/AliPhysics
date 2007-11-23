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
//  class AliCDBId						   //
//  Identity of an object stored into a database:  		   //
//  path, run validity range, version, subVersion 		   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBId.h"
#include <Riostream.h>
#include <TObjArray.h>
#include <TObjString.h>

ClassImp(AliCDBId)

//_____________________________________________________________________________
AliCDBId::AliCDBId():
fPath(), 
fRunRange(-1,-1), 
fVersion(-1), 
fSubVersion(-1),
fLastStorage("new")
{
// constructor

}

//_____________________________________________________________________________
AliCDBId::AliCDBId(const AliCDBId& other):
TObject(),
fPath(other.fPath), 
fRunRange(other.fRunRange),
fVersion(other.fVersion), 
fSubVersion(other.fSubVersion),
fLastStorage(other.fLastStorage)
{
// constructor

}

//_____________________________________________________________________________
AliCDBId::AliCDBId(const AliCDBPath& path, Int_t firstRun, Int_t lastRun, 
	Int_t version, Int_t subVersion):
fPath(path), 
fRunRange(firstRun, lastRun), 
fVersion(version), 
fSubVersion(subVersion),
fLastStorage("new")
{
// constructor

} 

//_____________________________________________________________________________
AliCDBId::AliCDBId(const AliCDBPath& path, const AliCDBRunRange& runRange, 
	Int_t version, Int_t subVersion):
fPath(path), 
fRunRange(runRange), 
fVersion(version),
fSubVersion(subVersion),
fLastStorage("new")
{
// constructor

}

//_____________________________________________________________________________
AliCDBId* AliCDBId::MakeFromString(const TString& idString)
{
// constructor from string 
// string has the format as the output of AliCDBId::ToString:
// path: "TRD/Calib/PIDLQ"; run range: [0,999999999]; version: v0_s0

	AliCDBId* id = new AliCDBId("a/b/c",-1,-1,-1,-1);
	
	TObjArray* arr1 = idString.Tokenize(';');
	TIter iter1(arr1);
	TObjString *objStr1 = 0;
	while((objStr1 = dynamic_cast<TObjString*>(iter1.Next()))) {
		TString buff(objStr1->GetName());
		
		if(buff.Contains("path:")) {
			TString path(buff(buff.First('\"')+1, buff.Length()-buff.First('\"')-2));
			id->SetPath(path.Data());
		
		} else if (buff.Contains("run range:")) {
			TString firstRunStr(buff(buff.Index('[')+1, buff.Index(',')-buff.Index('[')-1));
			TString lastRunStr(buff(buff.Index(',')+1, buff.Index(']')-buff.Index(',')-1));
			id->SetRunRange(firstRunStr.Atoi(), lastRunStr.Atoi());
		
		} else if (buff.Contains("version:")) {	
			if (buff.Contains("_s")) {
				TString versStr(buff(buff.Last('v')+1, buff.Index('_')-buff.Last('v')-1));
				TString subVersStr(buff(buff.Last('s')+1, buff.Length()-buff.Last('s')-1));
				id->SetVersion(versStr.Atoi());
				id->SetSubVersion(subVersStr.Atoi());
			} else {
				TString versStr(buff(buff.Last('v')+1, buff.Length()-buff.Last('v')-1));
				id->SetVersion(versStr.Atoi());
			}
		}
	
	}
	
	delete arr1;
	
	return id;

}

//_____________________________________________________________________________
AliCDBId::~AliCDBId() {
//destructor

}

//_____________________________________________________________________________
Bool_t AliCDBId::IsValid() const {
// validity check

	if (!(fPath.IsValid() && fRunRange.IsValid())) {
		return kFALSE;
	}
	
	// FALSE if doesn't have version but has subVersion
	return !(!HasVersion() && HasSubVersion());
}

//___________________________________________________________________________
Bool_t AliCDBId::IsEqual(const TObject* obj) const {
// check if this id is equal to other id (compares path, run range, versions)

        if (this == obj) {
                return kTRUE;
        }

        if (AliCDBId::Class() != obj->IsA()) {
                return kFALSE;
        }
        AliCDBId* other = (AliCDBId*) obj;
	return fPath.GetPath() == other->GetPath() && fRunRange.IsEqual(&other->GetAliCDBRunRange()) &&
		fVersion == other->GetVersion() && fSubVersion == other->GetSubVersion();
}

//_____________________________________________________________________________
TString AliCDBId::ToString() const {
// returns a string of Id data

	TString result = Form("path: \"%s\"; run range: [%d,%d]",
				GetPath().Data(), GetFirstRun(), GetLastRun());

	if(GetVersion() >= 0) result += Form("; version: v%d", GetVersion());
	if(GetSubVersion() >= 0) result += Form("_s%d", GetSubVersion());
	return result;
}

//_____________________________________________________________________________
void AliCDBId::Print(Option_t* /*option*/) const {
// Prints ToString()

	cout << ToString().Data() << endl;

}
