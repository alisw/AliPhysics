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

#include "AliMUON2DStoreValidator.h"

#include "AliLog.h"
#include "AliMUONCalibParam1I.h"
#include "AliMpIntPair.h"
#include "AliMpManuList.h"
#include "AliMUONV2DStore.h"
#include "TList.h"
#include "TObjArray.h"
#include "AliMUONCheckItem.h"
#include "AliMUONCheckItemIterator.h"
#include "AliMpDEManager.h"
#include "TObjString.h"
#include "AliMUONConstants.h"
#include "Riostream.h"

/// \class AliMUON2DStoreValidator
///
/// Determine which channels, manus, DEs, stations are missing
/// from a 2DStore. This is mainly to be used during (shuttle) preprocessing
/// to insure that what we'll put in the CDB is as complete as possible,
/// and to detect possible problem.
///
/// We made an effort to present the result of the validation in the most
/// concise way (i.e. if all channels of a DE are missing, we do not list 
/// them, we just write "DE dead" ;-) )
/// 
/// The list of missing things is kept in a structure of objects defined as :
/// 
/// fMissing = TObjArray[0..N tracking chambers]
/// fMissing[iChamber] = AliMUONCheckItem which contains n AliMUONCheckItem, 
/// where n is the number of DE for that chamber
/// fMissing[iChamber]->GetItem(de) = AliMUONCheckItem which contains m
/// AliMUONCheckItem where m is the number of Manu for that DE
/// fMissing[iChamber]->GetItem(de)->GetItem(manu) = AliMUONCheckItem which 
/// contains k TObjString = Form("%d",manuChannel)
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUON2DStoreValidator)
/// \endcond

//_____________________________________________________________________________
AliMUON2DStoreValidator::AliMUON2DStoreValidator() 
: TObject(),
  fManuList(0x0),
  fChambers(0x0)
{
    /// ctor
}

//_____________________________________________________________________________
AliMUON2DStoreValidator::~AliMUON2DStoreValidator()
{
  /// dtor
  delete fManuList;
  delete fChambers;
}

//_____________________________________________________________________________
AliMUONCheckItem* 
AliMUON2DStoreValidator::GetChamber(Int_t chamberID)
{
  /// Return (and create if not present) the given chamber
  /// chamberID in 0..NCh()
  
  if ( chamberID < 0 || chamberID >= AliMUONConstants::NCh() )
  {
    AliFatal(Form("Invalid chamber number %d",chamberID));
    return 0x0;
  }
  
  if (!fChambers) 
  {
    fChambers = new TObjArray(AliMUONConstants::NCh());
  }
    
  AliMUONCheckItem* chamber = 
    static_cast<AliMUONCheckItem*>(fChambers->At(chamberID));
  
  if (!chamber)
  {
    chamber = new AliMUONCheckItem(chamberID,
                                   AliMpDEManager::GetNofDEInChamber(chamberID),
                                   "Chamber");
    fChambers->AddAt(chamber,chamberID);
  }
  return chamber;
}

//_____________________________________________________________________________
AliMUONCheckItem* 
AliMUON2DStoreValidator::GetDE(Int_t detElemId)
{
  /// Return (and create if not present) a given detection element
  
  Int_t chamberID = AliMpDEManager::GetChamberId(detElemId);
  AliMUONCheckItem* chamber = GetChamber(chamberID);  
  AliMUONCheckItem* de = 
    static_cast<AliMUONCheckItem*>(chamber->GetItem(detElemId));
  if (!de)
  {
    AliDebug(3,Form("Did not find DE %4d into chamber %d, will create it",
                    detElemId,chamberID));
    de = new AliMUONCheckItem(detElemId,
                              AliMpManuList::NumberOfManus(detElemId),
                              "Detection Element");
    Bool_t ok = chamber->AddItem(detElemId,de);
    if (!ok)
    {
      AliError(Form("Could not add DE %4d into chamber %2d",detElemId,chamberID));
    }
  }
  return de;
}

//_____________________________________________________________________________
AliMUONCheckItem* 
AliMUON2DStoreValidator::GetManu(Int_t detElemId, Int_t manuId)
{
  /// Return (and create) a given manu
  
  AliMUONCheckItem* de = GetDE(detElemId);
  AliMUONCheckItem* manu = static_cast<AliMUONCheckItem*>(de->GetItem(manuId));
  if (!manu)
  {
    manu = new AliMUONCheckItem(manuId,AliMpManuList::NumberOfChannels(detElemId,manuId),"Manu");
    Bool_t ok = de->AddItem(manuId,manu);
    if (!ok)
    {
      AliError(Form("Could not add manu %4d into DE %4d",manuId,detElemId));
    }
    
  }
  return manu;
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::AddMissingChannel(Int_t detElemId, 
                                           Int_t manuId, Int_t manuChannel)
{
  /// Add one missing channel to the list of missing things
  
  AliDebug(3,Form("DE %4d Manu %4d Channel %2d is missing",
                  detElemId,manuId,manuChannel));

  AliMUONCheckItem* manu = GetManu(detElemId,manuId);
  Bool_t ok = manu->AddItem(manuChannel,new TObjString(Form("%2d",manuChannel)));
  if (!ok)
  {
    AliError(Form("Could not add channel %2d to manuId %4d in DE %4d",
                    manuChannel,manuId,detElemId));
  }
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::AddMissingManu(Int_t detElemId, Int_t manuId)
{
  /// Add one missing manu to the list of missing things
  
  AliDebug(3,Form("DE %4d Manu %4d is completely missing",
                  detElemId,manuId));

  Int_t n(AliMpManuList::NumberOfChannels(detElemId,manuId));

  for ( Int_t i = 0; i < n; ++i )
  {
    AddMissingChannel(detElemId,manuId,i);
  }
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::ReportManu(AliMUONCheckItem& manu)
{  
  /// Report list of missing channels from this manu
  
  TObjString* channel(0x0);
  AliMUONCheckItemIterator it(manu);
  
  it.First();
  
  while ( ( channel = static_cast<TObjString*>(it.Next()) ) )
  {
    cout << Form("\t\t\tChannel %s is missing",channel->GetString().Data()) << endl;
  }
  
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::ReportDE(AliMUONCheckItem& de)
{  
  /// Report list of missing manus from this de
  AliMUONCheckItem* manu(0x0);
  AliMUONCheckItemIterator it(de);
  
  cout << Form("\tDE %4d",de.GetID()) << endl;
  
  it.First();
  
  while ( ( manu = static_cast<AliMUONCheckItem*>(it.Next()) ) )
  {
    if ( manu->IsDead() )
    {
      cout << Form("\t\tManu %4d is missing",manu->GetID()) << endl;
    }
    else
    {
      ReportManu(*manu);
    }
  }
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::ReportChamber(AliMUONCheckItem& chamber)
{  
  /// Report list of missing de from this chamber
  
  AliMUONCheckItem* de(0x0);
  AliMUONCheckItemIterator it(chamber);
  
  it.First();
  
  while ( ( de = static_cast<AliMUONCheckItem*>(it.Next()) ) )
  {
    if ( de->IsDead() )
    {
      cout << Form("\tDE %4d is missing",de->GetID()) << endl;
    }
    else
    {
      ReportDE(*de);
    }
  }
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::Report(const TObjArray& chambers)
{
  /// Reports what is missing, trying to be as concise as possible.
  
  for ( Int_t iChamber = 0; iChamber <= chambers.GetLast(); ++iChamber )
  {
    AliMUONCheckItem* chamber = static_cast<AliMUONCheckItem*>(chambers.At(iChamber));
    if ( chamber )
    {
      if ( chamber->IsDead() )
      {
        cout << Form("Chamber %2d is missing",iChamber) << endl;
      }
      else
      {
        ReportChamber(*chamber);
      }
    }
  }
}

//_____________________________________________________________________________
TObjArray* 
AliMUON2DStoreValidator::Validate(const AliMUONV2DStore& store,
                                  Float_t invalidFloatValue)
{
  /// Validate the store. 
  /// The invalidFloatValue is used to decide if a store content value
  /// is valid or not.
  
  delete fChambers;
  fChambers = 0x0;
  
  if (!fManuList) fManuList = AliMpManuList::ManuList();

  // Now checks if some full manus are missing
  TIter next(fManuList);
  AliMpIntPair* p;

  while ( ( p = (AliMpIntPair*)next() ) )
  {
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
    AliMUONVCalibParam* test = 
      static_cast<AliMUONVCalibParam*>(store.Get(detElemId,manuId));
    if (!test)
    {
      // completely missing manu
      AddMissingManu(detElemId,manuId);
    }
    else
    {
      // manu is there, check all its channels
      for ( Int_t manuChannel = 0 ; manuChannel < test->Size(); ++manuChannel )
      {
        if ( AliMpManuList::DoesChannelExist(detElemId,manuId,manuChannel) &&
             ( test->ValueAsFloat(manuChannel,0) == invalidFloatValue ||
               test->ValueAsFloat(manuChannel,1) == invalidFloatValue ) )             
        {
          AddMissingChannel(detElemId,manuId,manuChannel);
        }
      }
    }
  }
  return fChambers;
}


