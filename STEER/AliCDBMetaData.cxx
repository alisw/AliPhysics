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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Object meta data: full description of a run dependent database object     // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TRegexp.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>

#include "AliCDBMetaData.h"
#include "AliLog.h"


ClassImp(AliCDBMetaData)


//_____________________________________________________________________________
AliCDBMetaData::AliCDBMetaData() :
  TObject(),
  fName(""),
  fFirstRun(-1),
  fLastRun(-1),
  fVersion(-1),
  fPeriod(-1),
  fFormat(""),
  fResponsible("Duck, Donald"),
  fExtraInfo("")
{
// default constructor
// the default values mean no selection
  DecodeName();
}

//_____________________________________________________________________________
AliCDBMetaData::AliCDBMetaData
         (const char* name, Int_t firstRun, Int_t lastRun, Int_t period, 
	  const char* objFormat, const char* responsible, 
	  const char* extraInfo):
  TObject(),
  fName(name),
  fFirstRun(firstRun),
  fLastRun(lastRun),
  fVersion(-1),
  fPeriod(period),
  fFormat(objFormat),
  fResponsible(responsible),
  fExtraInfo(extraInfo)
{
// constructor
  DecodeName();
}

//_____________________________________________________________________________
AliCDBMetaData::AliCDBMetaData(const AliCDBMetaData& entry) :
  TObject(entry),
  fName(entry.fName),
  fFirstRun(entry.fFirstRun),
  fLastRun(entry.fLastRun),
  fVersion(entry.fVersion),
  fPeriod(entry.fPeriod),
  fFormat(entry.fFormat),
  fResponsible(entry.fResponsible),
  fExtraInfo(entry.fExtraInfo)
{
// copy constructor
  DecodeName();
}

//_____________________________________________________________________________
AliCDBMetaData& AliCDBMetaData::operator = (const AliCDBMetaData& entry)
{
// assignment operator
  fName = entry.fName;
  fFirstRun = entry.fFirstRun;
  fLastRun = entry.fLastRun;
  fVersion = entry.fVersion;
  fPeriod=entry.fPeriod;
  fFormat=entry.fFormat;
  fResponsible=entry.fResponsible;
  fExtraInfo=entry.fExtraInfo;
  DecodeName();
  return *this;
}

//_____________________________________________________________________________
void AliCDBMetaData::EncodeName(){
// Encode name from single elements ("Detector", "DBType", "DetSpecType" -> "Detector/DBType/DetSpecType")   
   fName = fDetector+'/'+fDBType+'/'+fDetSpecType;
   if(fDBType == "*" && fDetSpecType == "*") fName = fDetector+'/'+'*';
   if(fDetector == "*" && fDBType == "*" && fDetSpecType == "*") fName = "*";

}

//_____________________________________________________________________________
void AliCDBMetaData::DecodeName(){
// Decode name into single elements ("Detector/DBType/DetSpecType" -> "Detector", "DBType", "DetSpecType")   

 if(fName==""){fDetector=""; fDBType=""; fDetSpecType=""; return;}

 while(fName.EndsWith("/")) fName.Remove(fName.Last('/'));
 while(fName.BeginsWith("/")) fName.Remove(fName.First('/'),1);
 
 // fName= "fDetector/fDBType/fDetSpecType
 int nslashes=fName.CountChar('/');
 
 if(nslashes>2){AliError("Wrong format!\n");fDetector=""; fDBType=""; fDetSpecType="";}

 if(nslashes == 0){
   if(fName == "*"){fDetector="*"; fDBType="*"; fDetSpecType="*";}
   else{AliError("Wrong format!\n"); fDetector=""; fDBType=""; fDetSpecType="";}
 }
 if(nslashes == 1){
   if(fName.EndsWith("*"))
     {fDetector=fName(0, fName.Index('/')); fDBType="*"; fDetSpecType="*";}
   else {AliError("Wrong format!\n"); fDetector=""; fDBType=""; fDetSpecType="";}
 }

 if(nslashes == 2){
   int firstsl=fName.First('/'), lastsl=fName.Last('/'), lgth=fName.Length();
   fDetector=fName(0, firstsl); 
   fDBType=fName(firstsl+1, lastsl-(firstsl+1));
   fDetSpecType=fName(lastsl+1, lgth-(lastsl+1)); 
 }
 EncodeName();
}

//_____________________________________________________________________________
Bool_t AliCDBMetaData::IsStrictlyValid(Int_t runNumber, AliCDBMetaData* metaData) const
{
// check if the object is valid for runNumber. TRUE if metaData version is equal to this's version 

  if ((fFirstRun >= 0) && (runNumber < fFirstRun)) return kFALSE;
  if ((fLastRun >= 0) && (runNumber > fLastRun)) return kFALSE;
  if (metaData) {
    if ((metaData->fVersion >= 0) && (metaData->fVersion != fVersion)) 
      return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliCDBMetaData::IsValid(Int_t runNumber, AliCDBMetaData* metaData) const
{
// check if the object is valid for runNumber. TRUE if metaData version less or equal wrt to this's

  if ((fFirstRun >= 0) && (runNumber < fFirstRun)) return kFALSE;
  if ((fLastRun >= 0) && (runNumber > fLastRun)) return kFALSE;
  if (metaData) {
    if ((metaData->fVersion >= 0) && (metaData->fVersion < fVersion)) 
      return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliCDBMetaData::Compare(const TObject* object) const
{
// check whether this is preferred to object

  if (!object || !object->InheritsFrom(AliCDBMetaData::Class())) return 1;
  if (fVersion < ((AliCDBMetaData*)object)->GetVersion()) return -1;
  if (fVersion > ((AliCDBMetaData*)object)->GetVersion()) return 1;
  return 0;
}

//_____________________________________________________________________________
Bool_t AliCDBMetaData::Matches(const char* name, Int_t runNumber) const
{
// check whether name and run number match with this meta data

  if ((fFirstRun >= 0) && (runNumber < fFirstRun)) return kFALSE;
  if ((fLastRun >= 0) && (runNumber > fLastRun)) return kFALSE;
  if (!TString(name).Contains(TRegexp(fName))) return kFALSE;
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t operator == (const AliCDBMetaData& entry1, const AliCDBMetaData& entry2)
{
// compare two DB entries

  if (strcmp(entry1.GetName(), entry2.GetName()) != 0) return kFALSE;
  if (entry1.GetFirstRun() != entry2.GetFirstRun()) return kFALSE;
  if (entry1.GetLastRun() != entry2.GetLastRun()) return kFALSE;
  if (entry1.GetVersion() != entry2.GetVersion()) return kFALSE;
  return kTRUE;
}


