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
// class that contains an object from the data base and knows about its      //
// validity range (meta data)                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliCDBEntry.h"

ClassImp(AliCDBEntry)


//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry() :
  TObject(),
  fObject(NULL),
  fMetaData()
{
// default constructor

}

//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry(const TObject* object, const AliCDBMetaData& metaData) :
  TObject(),
  fObject(object->Clone()),
  fMetaData(metaData)
{
// constructor

}

//_____________________________________________________________________________
AliCDBEntry::~AliCDBEntry()
{
// destructor

  delete fObject;
}


//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry(const AliCDBEntry& entry) :
  TObject(entry),
  fMetaData(entry.fMetaData)
{
// copy constructor

}

//_____________________________________________________________________________
AliCDBEntry& AliCDBEntry::operator = (const AliCDBEntry& entry)
{
// assignment operator

  delete fObject;
  fObject = entry.fObject->Clone();
  fMetaData = entry.fMetaData;
  return *this;
}



//_____________________________________________________________________________
const char* AliCDBEntry::GetName() const
{
// get the name

  return fMetaData.GetName();
}


//_____________________________________________________________________________
Int_t AliCDBEntry::Compare(const TObject* object) const
{
// check whether this is preferred to object

  if (!object || !object->InheritsFrom(AliCDBEntry::Class())) return 1;
  return fMetaData.Compare(&((AliCDBEntry*)object)->GetCDBMetaData());
}

