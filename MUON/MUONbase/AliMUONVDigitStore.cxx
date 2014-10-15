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

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONVDigitStore
///
/// Interface for a digit (or sdigit) container
///
/// It offers methods to Add, Find and Remove single elements, and
/// can create iterators to loop over (part of) the elements.
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONVDigitStore.h"

#include "AliLog.h"
#include "AliMUONVDigit.h"
#include <TClass.h>
#include <TString.h>
#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONVDigitStore)
/// \endcond

//_____________________________________________________________________________
AliMUONVDigitStore::AliMUONVDigitStore()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVDigitStore::~AliMUONVDigitStore()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t 
AliMUONVDigitStore::Add(TObject* object)
{
  /// Add an object, if it is of type AliMUONVDigit
  if (object)
  {
    AliMUONVDigit* digit = dynamic_cast<AliMUONVDigit*>(object);
    if (digit)
    {
      AliMUONVDigit* added = Add(*digit,AliMUONVDigitStore::kIgnore);
      if (!added)
      {
        AliError("Could not add digit through Add(TObject*) method");
      }
      else
      {
        return kTRUE;
      }
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONVDigitStore::Add(Int_t detElemId, 
                       Int_t manuId,
                       Int_t manuChannel,
                       Int_t cathode,
                       EReplacePolicy replace)
{
  /// Add a digit and return it
  AliMUONVDigit* digit = CreateDigit(detElemId,manuId,manuChannel,cathode);
  if (digit)
  {
    AliMUONVDigit* d = Add(*digit,replace);
    delete digit;
    return d;
  }
  return 0x0;
}

//____________________________________________________________________________
AliMUONVDigitStore* 
AliMUONVDigitStore::Create(const char* digitstoreclassname)
{
  /// Create a concrete digitStore, given its classname
  
  TClass* classPtr = TClass::GetClass(digitstoreclassname);
  if (!classPtr || !classPtr->InheritsFrom("AliMUONVDigitStore"))
  {
    return 0x0;
  }
  
  AliMUONVDigitStore* digitStore = 
    reinterpret_cast<AliMUONVDigitStore*>(classPtr->New());
  
  return digitStore;
}

//_____________________________________________________________________________
AliMUONVDigitStore*
AliMUONVDigitStore::Create(TTree& tree)
{
  /// Create store from the given tree (if possible).
  TString dataType = ( strcmp(tree.GetName(),"TreeD") == 0 ? "Digit" : 
                       (strcmp(tree.GetName(),"TreeS")== 9 ? "SDigit" : "")
                       );
  return static_cast<AliMUONVDigitStore*>(AliMUONVStore::Create(tree,dataType.Data()));
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONVDigitStore::FindObject(const TObject* object) const
{
  /// Find an object, if of AliMUONVDigit type.
  const AliMUONVDigit* digit = dynamic_cast<const AliMUONVDigit*>(object);
  if (digit)
  {
    return FindObject(digit->GetUniqueID());
  }
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONVDigitStore::FindObject(UInt_t uniqueID) const
{
  /// Find digit by its uniqueID
  
  return FindObject(AliMUONVDigit::DetElemId(uniqueID),
                    AliMUONVDigit::ManuId(uniqueID),
                    AliMUONVDigit::ManuChannel(uniqueID),
                    AliMUONVDigit::Cathode(uniqueID));
}

//_____________________________________________________________________________
Int_t 
AliMUONVDigitStore::GetSize(Int_t detElemId, Int_t cathode) const
{
  /// Return the number of digits we have for a given detection element
  TIter next(CreateIterator(detElemId,detElemId,cathode));
  Int_t n(0);
  while ( ( next() ) )
  {
    ++n;
  }
  return n;
}

