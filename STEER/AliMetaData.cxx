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
// base class of the metadata of run dependent objects                       //
// Derived classes: AliObjectMetaData, AliSelectionMetaData		     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TRegexp.h>

#include "AliMetaData.h"
#include "AliLog.h"


ClassImp(AliMetaData)


//_____________________________________________________________________________
AliMetaData::AliMetaData() :
  TObject(),
  fName(""),
  fFirstRun(-1),
  fLastRun(-1),
  fVersion(-1)
{
// default constructor
// the default values mean no selection
  DecodeName();
}

//_____________________________________________________________________________
AliMetaData::AliMetaData(const char* name, Int_t firstRun, Int_t lastRun, 
			 Int_t version) :
  TObject(),
  fName(name),
  fFirstRun(firstRun),
  fLastRun(lastRun),
  fVersion(version)
{
// constructor
  DecodeName();
}


//_____________________________________________________________________________
AliMetaData::AliMetaData(const AliMetaData& entry) :
  TObject(entry),
  fName(entry.fName),
  fFirstRun(entry.fFirstRun),
  fLastRun(entry.fLastRun),
  fVersion(entry.fVersion)
{
// copy constructor
  DecodeName();
}

//_____________________________________________________________________________
AliMetaData& AliMetaData::operator = (const AliMetaData& entry)
{
// assignment operator

  fName = entry.fName;
  fFirstRun = entry.fFirstRun;
  fLastRun = entry.fLastRun;
  fVersion = entry.fVersion;
  DecodeName();
  return *this;
}



//_____________________________________________________________________________
const char* AliMetaData::GetName() const
{
// get the name ("Detector/DBType/DetSpecType", example: "ZDC/Calib/Pedestals")

  return fName.Data();
}

//_____________________________________________________________________________
const char* AliMetaData::GetDetector() const
{
// get the detector's name (ZDC,ITS ...)

  return fDetector.Data();
}

//_____________________________________________________________________________
const char* AliMetaData::GetDBType() const
{
// get the database type (Calib, Align ...)

  return fDBType.Data();
}

//_____________________________________________________________________________
const char* AliMetaData::GetDetSpecType() const
{
// get the detector's specific type name (Pedestals, GainConst, DeadChannelMaps...)

  return fDetSpecType.Data();
}

//_____________________________________________________________________________
void AliMetaData::EncodeName(){
// Encode name from single elements ("Detector", "DBType", "DetSpecType" -> "Detector/DBType/DetSpecType")   
   fName = fDetector+'/'+fDBType+'/'+fDetSpecType;
   if(fDBType == "*" && fDetSpecType == "*") fName = fDetector+'/'+'*';
   if(fDetector == "*" && fDBType == "*" && fDetSpecType == "*") fName = "*";

}

//_____________________________________________________________________________
void AliMetaData::DecodeName(){
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
Bool_t AliMetaData::IsStrictlyValid(Int_t runNumber, AliMetaData* metaData) const
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
Bool_t AliMetaData::IsValid(Int_t runNumber, AliMetaData* metaData) const
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
Int_t AliMetaData::Compare(const TObject* object) const
{
// check whether this is preferred to object

  if (!object || !object->InheritsFrom(AliMetaData::Class())) return 1;
  if (fVersion < ((AliMetaData*)object)->GetVersion()) return -1;
  if (fVersion > ((AliMetaData*)object)->GetVersion()) return 1;
  return 0;
}

//_____________________________________________________________________________
Bool_t AliMetaData::Matches(const char* name, Int_t runNumber) const
{
// check whether name and run number match with this meta data

  if ((fFirstRun >= 0) && (runNumber < fFirstRun)) return kFALSE;
  if ((fLastRun >= 0) && (runNumber > fLastRun)) return kFALSE;
  if (!TString(name).Contains(TRegexp(fName))) return kFALSE;
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t operator == (const AliMetaData& entry1, const AliMetaData& entry2)
{
// compare two DB entries

  if (strcmp(entry1.GetName(), entry2.GetName()) != 0) return kFALSE;
  if (entry1.GetFirstRun() != entry2.GetFirstRun()) return kFALSE;
  if (entry1.GetLastRun() != entry2.GetLastRun()) return kFALSE;
  if (entry1.GetVersion() != entry2.GetVersion()) return kFALSE;
  return kTRUE;
}

