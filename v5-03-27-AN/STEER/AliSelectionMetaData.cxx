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
// subsample of the object metadata, used to retrieve                        //
// a stored database object:                                                 // 
// name="Detector/DBType/DetSpecType", run validity, version                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TSystem.h>

#include "AliSelectionMetaData.h"
#include "AliMetaData.h"
#include "AliLog.h"


ClassImp(AliSelectionMetaData)


//_____________________________________________________________________________
AliSelectionMetaData::AliSelectionMetaData() :
  AliMetaData()
{
// default constructor
// the default values mean no selection
}

//_____________________________________________________________________________
AliSelectionMetaData::AliSelectionMetaData(const char* name, Int_t firstRun, Int_t lastRun, Int_t Version) :
  AliMetaData(name, firstRun, lastRun, Version)
{
// constructor
}

//_____________________________________________________________________________
AliSelectionMetaData::AliSelectionMetaData(const AliSelectionMetaData& entry) :
  AliMetaData(entry)
{
// copy constructor
}

//_____________________________________________________________________________
AliSelectionMetaData::AliSelectionMetaData(const AliMetaData& entry) :
  AliMetaData(entry)
{
// constructor of AliSelectionMetaData from AliMetaData
}


//_____________________________________________________________________________
AliSelectionMetaData& AliSelectionMetaData::operator = (const AliSelectionMetaData& entry)
{
// assignment operator
  fName = entry.fName,
  fFirstRun = entry.fFirstRun;
  fLastRun = entry.fLastRun;
  fVersion = entry.fVersion;
  return *this;
}

