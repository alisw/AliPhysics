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

#include "AliObjectMetaData.h"
#include "AliMetaData.h"
#include "AliLog.h"


ClassImp(AliObjectMetaData)


//_____________________________________________________________________________
AliObjectMetaData::AliObjectMetaData() :
  AliMetaData(),
  fPeriod(-1),
  fFormat(""),
  fResponsible("Duck, Donald"),
  fExtraInfo("")
{
// default constructor
// the default values mean no selection
}

//_____________________________________________________________________________
AliObjectMetaData::AliObjectMetaData
         (const char* name, Int_t firstRun, Int_t lastRun, Int_t period, 
	  const char* objFormat, const char* responsible, 
	  const char* extraInfo):
  AliMetaData(name, firstRun, lastRun),
  fPeriod(period),
  fFormat(objFormat),
  fResponsible(responsible),
  fExtraInfo(extraInfo)
{
// constructor
}

//_____________________________________________________________________________
AliObjectMetaData::AliObjectMetaData(const AliObjectMetaData& entry) :
  AliMetaData(entry),
  fPeriod(entry.fPeriod),
  fFormat(entry.fFormat),
  fResponsible(entry.fResponsible),
  fExtraInfo(entry.fExtraInfo)
{
// copy constructor
}

//_____________________________________________________________________________
AliObjectMetaData& AliObjectMetaData::operator = (const AliObjectMetaData& entry)
{
// assignment operator
  fName = entry.fName;
  fFirstRun = entry.fFirstRun;
  fLastRun = entry.fLastRun;
  fPeriod=entry.fPeriod;
  fFormat=entry.fFormat;
  fResponsible=entry.fResponsible;
  fExtraInfo=entry.fExtraInfo;
  DecodeName();
  return *this;
}

//_____________________________________________________________________________
const int AliObjectMetaData::GetPeriod() const
{
// get the beam period

  return fPeriod;
}

//_____________________________________________________________________________
const char* AliObjectMetaData::GetFormat() const
{
// get the object's format

  return fFormat.Data();
}

//_____________________________________________________________________________
const char* AliObjectMetaData::GetResponsible() const
{
// get the object's responsible (the person who made it)

  return fResponsible.Data();
}

//_____________________________________________________________________________
const char* AliObjectMetaData::GetExtraInfo() const
{
// get the object's extra info

  return fExtraInfo.Data();
}

