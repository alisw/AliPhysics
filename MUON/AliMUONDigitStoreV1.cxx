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

#include "AliMUONDigitStoreV1.h"

//-----------------------------------------------------------------------------
/// \class AliMUONDigitStoreV1
///
/// (Legacy) Implementation of AliMUONVDigitStore. 
/// Called legacy as the internal structure corresponds to what we
/// used to write as MUON.(S)Digits.root files, before the switch 
/// to data stores.
///
// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliLog.h"
#include "AliMUONDigit.h"
#include "AliMUONDigitStoreV1Iterator.h"
#include "AliMUONTOTCAStoreIterator.h"
#include "AliMUONTreeManager.h"
#include "AliMpConstants.h"
#include "AliMpDEManager.h"
#include "AliMpDEManager.h"
#include "AliMpStationType.h"
#include <Riostream.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONDigitStoreV1)
/// \endcond

namespace
{
  TString BaseName(const TString& name)
  {
    if ( name == "TreeS" ) return "MUONSDigit";
    if ( name == "TreeD" ) return "MUONDigit";
    return "";
  }
}

//_____________________________________________________________________________
AliMUONDigitStoreV1::AliMUONDigitStoreV1()
: AliMUONVDigitStore(), 
fDigits(new TObjArray(AliMpConstants::NofChambers())),
fChamberDigits(0x0)
{
  /// ctor
  fDigits->SetOwner(kTRUE);
  for ( Int_t i = 0; i < fDigits->GetSize(); ++i )
  {
    TClonesArray* tca = new TClonesArray("AliMUONDigit",100);
    tca->SetOwner(kTRUE);
    fDigits->AddAt(tca,i);
  }
  Clear();
  AliDebug(1,"");
}

//_____________________________________________________________________________
AliMUONDigitStoreV1::AliMUONDigitStoreV1(const AliMUONDigitStoreV1&)
: AliMUONVDigitStore(), 
fDigits(0x0),
fChamberDigits(0x0)
{
  /// copy ctor
  AliError("Please implement me");
}

//_____________________________________________________________________________
AliMUONDigitStoreV1& 
AliMUONDigitStoreV1::operator=(const AliMUONDigitStoreV1&)
{
  /// assignement operator
  AliError("Please implement me");
  return *this;
}

//_____________________________________________________________________________
AliMUONDigitStoreV1::~AliMUONDigitStoreV1()
{
  /// dtor
  delete fDigits;
}


//_____________________________________________________________________________
void 
AliMUONDigitStoreV1::Clear(Option_t* /*opt*/)
{
  /// Clear the tclonesarray, but keep the tobjarray's size constant.

  for ( Int_t i = 0; i <= fDigits->GetLast(); ++i ) 
  {
    ChamberDigits(i)->Clear("C");
  }  
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreV1::Add(const AliMUONVDigit& vdigit, EReplacePolicy replace)
{
  /// Try to add a digit to the store. Return whether the try was successfull
  /// or not.
  /// 
  /// If the digit is already there, the action taken depends on "replace"
  /// kAllow -> replacement will occur (i.e. return kTRUE)
  /// kDeny -> replacement will *not* occur (and returned value is kFALSE)
  /// kMerge -> both digits will be merged into one (return kTRUE)
  ///
  
  const AliMUONDigit* digit = dynamic_cast<const AliMUONDigit*>(&vdigit);
  
  if (!digit)
  {
    AliError(Form("Digit is not of the expected type (%s vs AliMUONdigit)",
                  vdigit.ClassName()));
    return 0x0;
  }
  
  Int_t index(-1);
  
  if ( replace != kIgnore ) 
  {
    AliMUONVDigit* alreadyThere = Find(*digit,index);
  
    if ( alreadyThere ) 
    {
      if ( replace == kDeny ) 
      {
        return 0x0;
      }
      if ( replace == kMerge ) 
      {
        alreadyThere->MergeWith(*digit);
        return alreadyThere;
      }
    }
  }
  
  Int_t detElemId = digit->DetElemId();
  Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
  TClonesArray* array = ChamberDigits(iChamber);
  if ( index < 0 )
  {
    index = array->GetLast()+1;
  }
  return (new((*array)[index]) AliMUONDigit(*digit));
}

//_____________________________________________________________________________
Bool_t
AliMUONDigitStoreV1::Connect(TTree& tree, Bool_t alone) const
{
  /// Connect this to the tree.
  
  AliMUONTreeManager tman;
  Bool_t ok(kTRUE);
  
  TString baseName(BaseName(tree.GetName()));
                     
  // Search for branch MUON(S)Digits1 to know if we need to set branch addresses
  // or to make them.
  TBranch* branch = tree.GetBranch(Form("%ss%d",baseName.Data(),1));
  
  Bool_t isMaking = (branch==0);
  
  if ( isMaking ) 
  {
    for ( Int_t i = 0; i < AliMpConstants::NofChambers(); ++i ) 
    {
      TString branchName(Form("%ss%d",baseName.Data(),i+1));
      ok = ok && tman.MakeBranch(tree,ClassName(),"TClonesArray",
                                 branchName.Data(),ChamberDigitsPtr(i));
    }
  }
  else
  {
    if ( alone && baseName != "MUONSDigit" ) 
    {
      // TreeS only has digits, so there's not need to play the branch status
      // game
      tman.UpdateBranchStatuses(tree,baseName.Data());
    }
    for ( Int_t i = 0; i < AliMpConstants::NofChambers(); ++i ) 
    {
      TString branchName(Form("%ss%d",baseName.Data(),i+1));
      ok = ok && tman.SetAddress(tree,branchName.Data(),ChamberDigitsPtr(i));
    }
  }
  
  return ok;
}

//_____________________________________________________________________________
TObject**
AliMUONDigitStoreV1::ChamberDigitsPtr(Int_t chamberId) const
{
  /// Get the address of the TClonesArray storing digits for chamberId.

  return fDigits->GetObjectRef(fDigits->UncheckedAt(chamberId));

  TObject* object = fDigits->At(chamberId);

  if (!object)
  {
    AliError(Form("Cannot get digits for chamberId=%d",chamberId));
    return 0x0;
  }
  else
  {
    return fDigits->GetObjectRef(object);
  }
}

//_____________________________________________________________________________
TClonesArray*
AliMUONDigitStoreV1::ChamberDigits(Int_t chamberId)
{
  /// Returns the tclonesarray storing digits for chamberId
  return static_cast<TClonesArray*>(fDigits->At(chamberId));
}

//_____________________________________________________________________________
const TClonesArray*
AliMUONDigitStoreV1::ChamberDigits(Int_t chamberId) const
{
  /// Returns the tclonesarray storing digits for chamberId
  return static_cast<TClonesArray*>(fDigits->At(chamberId));
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreV1::CreateDigit(Int_t detElemId, Int_t manuId,
                                 Int_t manuChannel, Int_t cathode) const
{
  return new AliMUONDigit(detElemId,manuId,manuChannel,cathode);
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreV1::Remove(AliMUONVDigit& digit)
{
  /// Remove one digit, and returns it, thus returning 0x0 if digit
  /// is not present.
  
  Int_t index;
  AliMUONVDigit* d(0x0);
  if ( ( d = FindIndex(digit.DetElemId(),digit.ManuId(),
                       digit.ManuChannel(),digit.Cathode(),index) ) )
  {
    Int_t iChamber = AliMpDEManager::GetChamberId(digit.DetElemId());
    TClonesArray* array = ChamberDigits(iChamber);
    array->RemoveAt(index);
    array->Compress();
  }
  return d;
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreV1::Find(const AliMUONVDigit& digit, Int_t& index) const
{
  /// Find a digit, and return its index.
  return FindIndex(digit.DetElemId(),digit.ManuId(),digit.ManuChannel(),digit.Cathode(),index);
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreV1::FindObject(Int_t detElemId, Int_t manuId, 
                                Int_t manuChannel, Int_t cathode) const
{
  /// Find a (trigger) digit
  Int_t index;
  return FindIndex(detElemId,manuId,manuChannel,cathode,index);
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreV1::FindIndex(Int_t detElemId, Int_t manuId, 
                               Int_t manuChannel, Int_t cathode, Int_t& index) const
{
  /// Find and return the index of a digit
  
  Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
  const TClonesArray* array = ChamberDigits(iChamber);
  if (!array) return 0x0;
  TIter next(array);
  AliMUONVDigit* digit;
  index=0;
  while ( ( digit = static_cast<AliMUONVDigit*>(next()) ) )
  {
    if ( digit->DetElemId() == detElemId &&
         digit->ManuId() == manuId &&
         digit->ManuChannel() == manuChannel &&
         digit->Cathode() == cathode ) 
    {
      return digit;
    }
    ++index;
  }
  return 0x0;
  
}

//_____________________________________________________________________________
TIterator* 
AliMUONDigitStoreV1::CreateIterator() const
{
  /// Return an iterator on the full store
  return new AliMUONTOTCAStoreIterator(fDigits,0,13);
}

//_____________________________________________________________________________
TIterator* 
AliMUONDigitStoreV1::CreateTrackerIterator() const
{
  /// Return an iterator on the tracker part of the store
  return new AliMUONTOTCAStoreIterator(fDigits,0,9);
}

//_____________________________________________________________________________
TIterator* 
AliMUONDigitStoreV1::CreateTriggerIterator() const
{
  /// Return an iterator on the trigger part of the store
  return new AliMUONTOTCAStoreIterator(fDigits,10,13);
}

//_____________________________________________________________________________
TIterator* 
AliMUONDigitStoreV1::CreateIterator(Int_t firstDetElemId, Int_t lastDetElemId,
                                    Int_t cathode) const
{
  /// Return an iterator on part of the store
  return new AliMUONDigitStoreV1Iterator(fDigits,firstDetElemId,lastDetElemId,cathode);
}

//_____________________________________________________________________________
Int_t
AliMUONDigitStoreV1::GetSize() const
{
  /// Return the number of digits we store
  Int_t n(0);
  
  for ( Int_t i = 0; i <= fDigits->GetLast(); ++i ) 
  {
    n += ChamberDigits(i)->GetEntries();
  }  
  return n;
}

//_____________________________________________________________________________
Bool_t 
AliMUONDigitStoreV1::HasMCInformation() const
{
  /// As this class is legacy, don't care about looping and loosing a bit of
  /// time...
  TIter next(CreateIterator());
  AliMUONVDigit* digit;
  while ( ( digit = static_cast<AliMUONVDigit*>(next()) ) )
  {
    if ( digit->HasMCInformation() ) return kTRUE;
  }
  return kFALSE;
}

