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
// meta data of run dependent objects                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TRegexp.h>

#include "AliMetaData.h"


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

}

//_____________________________________________________________________________
AliMetaData& AliMetaData::operator = (const AliMetaData& entry)
{
// assignment operator

  fName = entry.fName;
  fFirstRun = entry.fFirstRun;
  fLastRun = entry.fLastRun;
  fVersion = entry.fVersion;
  return *this;
}



//_____________________________________________________________________________
const char* AliMetaData::GetName() const
{
// get the name

  return fName.Data();
}


//_____________________________________________________________________________
Bool_t AliMetaData::IsValid(Int_t runNumber, AliMetaData* metaData) const
{
// check the validity of the object

  if ((fFirstRun >= 0) && (runNumber < fFirstRun)) return kFALSE;
  if ((fLastRun >= 0) && (runNumber > fLastRun)) return kFALSE;
  if (metaData) {
    if ((metaData->fVersion >= 0) && (metaData->fVersion != fVersion)) 
      return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliMetaData::Compare(const TObject* object) const
{
// check whether this is prefered to object

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
  if (TString(name).Contains(TRegexp(fName))) return kTRUE;
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
