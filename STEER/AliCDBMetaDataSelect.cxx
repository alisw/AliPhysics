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
// subsample of the CDB metadata, used to retrieve                           //
// a stored database object:                                                 // 
// name="Detector/DBType/DetSpecType", run validity, version                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TSystem.h>

#include "AliCDBMetaDataSelect.h"
#include "AliCDBMetaData.h"
#include "AliLog.h"


ClassImp(AliCDBMetaDataSelect)


//_____________________________________________________________________________
AliCDBMetaDataSelect::AliCDBMetaDataSelect() :
  AliCDBMetaData()
{
// default constructor
// the default values mean no selection
}

//_____________________________________________________________________________
AliCDBMetaDataSelect::AliCDBMetaDataSelect(const char* name, Int_t firstRun, Int_t lastRun, Int_t Version) :
  AliCDBMetaData(name, firstRun, lastRun)
{
// constructor
  fVersion=Version;
}

//_____________________________________________________________________________
AliCDBMetaDataSelect::AliCDBMetaDataSelect(const AliCDBMetaDataSelect& entry) :
  AliCDBMetaData(entry)
{
// copy constructor
}

//_____________________________________________________________________________
AliCDBMetaDataSelect::AliCDBMetaDataSelect(const AliCDBMetaData& entry) :
  AliCDBMetaData(entry)
{
// constructor of AliCDBMetaDataSelect from AliCDBMetaData
 fPeriod=-1;
 fFormat="";
 fResponsible="Duck, Donald";
 fExtraInfo="";
}


//_____________________________________________________________________________
AliCDBMetaDataSelect& AliCDBMetaDataSelect::operator = (const AliCDBMetaDataSelect& entry)
{
// assignment operator
  fName = entry.fName,
  fFirstRun = entry.fFirstRun;
  fLastRun = entry.fLastRun;
  fVersion = entry.fVersion;
  return *this;
}

