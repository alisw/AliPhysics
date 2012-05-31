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
/// \class AliMUONDigitStoreVImpl
///
/// Base implementation of VDigitStore, where digits are simply put
/// within a single TClonesArray (to get compact output on disk) and
/// fast search is achieved by using a separate index (a 2Dmap)
///
/// Note that this class is a base implementation only, because to add
/// a digit to a TClonesArray, you need to give the concrete class
/// name (made in subclass by overriding :
///
/// AliMUONVDigit* AddConcreteDigit(TClonesArray&,const AliMUONVDigit&, Int_t)
///
/// and
///
/// AliMUONVDigit* CreateDigit((Int_t detElemId, Int_t manuId,
/// Int_t manuChannel, Int_t cathode) const
/// 
/// methods.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONDigitStoreVImpl.h"

#include "AliLog.h"
#include "AliMUONDigitStoreVImplIterator.h"
#include "AliMUON2DMap.h"
#include "AliMUONTreeManager.h"
#include "AliMUONCalibParamNI.h"
#include "AliMpConstants.h"
#include <TClonesArray.h>
#include <TTree.h>
#include <TString.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONDigitStoreVImpl)
/// \endcond

namespace
{
  TString BaseName(const TString& name)
  {
    /// Name of the branch, depending on the tree name
    if ( name == "TreeS" ) return "MUONSDigit";
    if ( name == "TreeD" ) return "MUONDigit";
    return "";
  }
  
  Int_t InternalManuId(Int_t cathode, Int_t manuId)
  {
    /// Very local convention to insure that trigger digits are not
    /// mixed (as by default a same manuId can be in cath 0 or 1, which
    /// is never the case for tracker)
    ///
    /// WARNING : the resulting manuId must still be contained within 16 bits !
    ///
    return manuId | ( cathode << 15 );
  }
  
}

//_____________________________________________________________________________
AliMUONDigitStoreVImpl::AliMUONDigitStoreVImpl(const char* concreteClassName)
: AliMUONVDigitStore(),
  fDigits(new TClonesArray(concreteClassName,100)),
  fMap(0x0),
  fIndexed(kFALSE)
{
    /// ctor
}

//_____________________________________________________________________________
AliMUONDigitStoreVImpl::AliMUONDigitStoreVImpl(const AliMUONDigitStoreVImpl&)
: AliMUONVDigitStore(),
fDigits(0x0),
fMap(0x0),
fIndexed(kFALSE)
{
  /// copy ctor
  AliError("Please implement me");
}

//_____________________________________________________________________________
AliMUONDigitStoreVImpl& 
AliMUONDigitStoreVImpl::operator=(const AliMUONDigitStoreVImpl&)
{
  /// assignement operator
  AliError("Please implement me");
  return *this;
}


//_____________________________________________________________________________
AliMUONDigitStoreVImpl::~AliMUONDigitStoreVImpl()
{
  /// dtor
  delete fDigits;
  delete fMap;
}

//_____________________________________________________________________________
Bool_t 
AliMUONDigitStoreVImpl::Connect(TTree& tree, Bool_t alone) const
{
  /// Connect this to the tree.
  
  TString branchName(BaseName(tree.GetName()));
  branchName += ".";
  AliMUONTreeManager tman;
  Bool_t ok;
  
  if (tree.GetBranch(branchName.Data()))
  {
    if ( alone ) tman.UpdateBranchStatuses(tree,BaseName(tree.GetName()));
    ok = tman.SetAddress(tree,branchName.Data(),
                         const_cast<TClonesArray**>(&fDigits));
  }
  else
  {
    ok = tman.MakeBranch(tree,ClassName(),"TClonesArray",branchName.Data(),
                         const_cast<TClonesArray**>(&fDigits));
  }
  
  return ok;
}

//_____________________________________________________________________________
void 
AliMUONDigitStoreVImpl::Clear(Option_t*)
{
  /// Clear the internal digit array AND the index
  fDigits->Clear("C");
  ClearIndex();
}

//_____________________________________________________________________________
void
AliMUONDigitStoreVImpl::ClearIndex()
{
  /// Clear our internal index
  if ( fMap ) 
  {
    fMap->Clear();
  }
  else
  {
    fMap = new AliMUON2DMap(true);
  }
  fIndexed = kFALSE;
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreVImpl::Add(const AliMUONVDigit& vdigit, EReplacePolicy replace)
{
  /// Try to add a digit to the store. Return whether the try was successfull
  /// or not.
  /// 
  /// If the digit is already there, the action taken depends on "replace"
  /// kAllow -> replacement will occur (i.e. return kTRUE)
  /// kDeny -> replacement will *not* occur (and returned value is kFALSE)
  /// kMerge -> both digits will be merged into one (return kTRUE)
  ///
  
  if ( replace != kIgnore )
  {
    AliMUONVDigit* alreadyThere = Find(vdigit);
    if ( alreadyThere ) 
    {
      if ( replace == kDeny ) return 0x0;
      if ( replace == kMerge )
      {
        alreadyThere->MergeWith(vdigit);
        return alreadyThere;
      }
    }
  }
  
  
  Int_t n = fDigits->GetLast()+1;

  AliMUONVDigit* d = AddConcreteDigit(*fDigits,vdigit,n);
  
  if ( d )
  {
    UpdateIndex(*d,n);
  }
  
  return d;
}

//_____________________________________________________________________________
TIterator* 
AliMUONDigitStoreVImpl::CreateIterator() const
{
  /// Create an iterator over the full store
  return fDigits->MakeIterator();
}

//_____________________________________________________________________________
TIterator* 
AliMUONDigitStoreVImpl::CreateIterator(Int_t firstDetElemId, 
                                    Int_t lastDetElemId,
                                    Int_t cathode) const
{
  /// Create an iterator on a given part of the store
  (const_cast<AliMUONDigitStoreVImpl*>(this))->ReIndex();

  return new AliMUONDigitStoreVImplIterator(this,firstDetElemId,lastDetElemId,cathode);
}

//_____________________________________________________________________________
TIterator*
AliMUONDigitStoreVImpl::CreateTrackerIterator() const
{
  /// Create an iterator to loop over tracker digits only
  
  (const_cast<AliMUONDigitStoreVImpl*>(this))->ReIndex();

  return new AliMUONDigitStoreVImplIterator(this,100,1025);
}

//_____________________________________________________________________________
TIterator* 
AliMUONDigitStoreVImpl::CreateTriggerIterator() const
{
  /// Create an iterator to loop over trigger digits only
  (const_cast<AliMUONDigitStoreVImpl*>(this))->ReIndex();

  return new AliMUONDigitStoreVImplIterator(this,1100,1417);
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreVImpl::Find(const AliMUONVDigit& digit) const
{
  /// Find a given digit
  /// Note that we only check for the id of the digit (i.e. de,manu,...)
  /// not for the actual content (charge, ...) to decide whether
  /// it's the same digit or not
  
  return FindObject(digit.DetElemId(),digit.ManuId(),digit.ManuChannel(),
                    digit.Cathode());
}

//_____________________________________________________________________________
void
AliMUONDigitStoreVImpl::ReIndex()
{
  /// Recompute the fMap, which map (de,manu,ch) to an index within
  /// the fDigits array
  
  if ( fIndexed ) return;
  
  ClearIndex();
  
  TIter next(fDigits);
  AliMUONVDigit* d;
  Int_t digitIndex(0);
  
  while ( ( d = static_cast<AliMUONVDigit*>(next()) ) )
  {
    UpdateIndex(*d,digitIndex++);
  }
  
  fIndexed = kTRUE;
}

//_____________________________________________________________________________
void
AliMUONDigitStoreVImpl::UpdateIndex(const AliMUONVDigit& digit, Int_t index)
{
  /// Update the internal index given this new digit
  if (!fMap) fMap = new AliMUON2DMap(true);
  
  Int_t manuId = InternalManuId(digit.Cathode(),digit.ManuId());
                        
  AliMUONVCalibParam* param = 
  static_cast<AliMUONVCalibParam*>
  (fMap->FindObject(digit.DetElemId(),manuId));

  if (!param)
  {
    param = new AliMUONCalibParamNI(1,64,digit.DetElemId(),manuId,-1);
    fMap->Add(param);
  }
  param->SetValueAsInt(digit.ManuChannel(),0,index);
  fIndexed = kTRUE;
}

//_____________________________________________________________________________
Int_t
AliMUONDigitStoreVImpl::FindIndex(Int_t detElemId, Int_t internalManuId,
                                  Int_t manuChannel) const
{
  /// Find the index of a given (de,internalManu,ch) triplet
  
  AliMUONVCalibParam* param = 
  static_cast<AliMUONVCalibParam*>
  (fMap->FindObject(detElemId,internalManuId));
  
  if (param)
  {
    return param->ValueAsInt(manuChannel);
  }
  
  return -1;
}

//_____________________________________________________________________________
Int_t
AliMUONDigitStoreVImpl::FindIndex(const AliMUONVDigit& digit) const
{
  /// Find the index of a given digit
  return FindIndex(digit.DetElemId(),
                   InternalManuId(digit.Cathode(),digit.ManuId()),
                   digit.ManuChannel());
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreVImpl::FindObject(UInt_t uniqueID) const
{
  /// Find digit by its uniqueID
  
  return FindObject(AliMUONVDigit::DetElemId(uniqueID),
                    AliMUONVDigit::ManuId(uniqueID),
                    AliMUONVDigit::ManuChannel(uniqueID),
                    AliMUONVDigit::Cathode(uniqueID));
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreVImpl::FindObject(Int_t detElemId, Int_t manuId, Int_t manuChannel, 
                                Int_t cathode) const
{
  /// Find a digit

  (const_cast<AliMUONDigitStoreVImpl*>(this))->ReIndex();

  Int_t index = FindIndex(detElemId,InternalManuId(cathode,manuId),manuChannel);
  
  if (index>=0 ) 
  {
    return static_cast<AliMUONVDigit*>(fDigits->UncheckedAt(index));
  }
  
  return 0x0;
}

//_____________________________________________________________________________
Int_t 
AliMUONDigitStoreVImpl::GetSize() const
{
  /// Return the number of digits we hold
  return fDigits->GetLast()+1;
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreVImpl::Remove(AliMUONVDigit& digit)
{
  /// Remove a digit
  AliMUONVDigit* d = static_cast<AliMUONVDigit*>(fDigits->Remove(&digit));
  if (d) 
  {
    UpdateIndex(*d,-1);
  }
  return d;
}

