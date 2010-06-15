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
#include "AliMUONCheckItem.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include <Riostream.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>

//-----------------------------------------------------------------------------
/// \class AliMUON2DStoreValidator
///
/// Determine which channels, manus, DEs, stations are missing
/// from a VStore, which must be 2D, and the 2 dimensions must be
/// (detElemId,manuId).
/// This is mainly to be used during (shuttle) preprocessing
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
///
/// fMissing[iChamber] = AliMUONCheckItem which contains n AliMUONCheckItem, 
/// where n is the number of DE for that chamber
///
/// fMissing[iChamber]->GetItem(de) = AliMUONCheckItem which contains m
/// AliMUONCheckItem where m is the number of Manu for that DE
///
/// fMissing[iChamber]->GetItem(de)->GetItem(manu) = AliMUONCheckItem which 
/// contains k TObjString = Form("%d",manuChannel)
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUON2DStoreValidator)
/// \endcond

//_____________________________________________________________________________
AliMUON2DStoreValidator::AliMUON2DStoreValidator() 
: TObject(),
  fChambers(0x0),
  fStatus(0x0)
{
    /// ctor
}

//_____________________________________________________________________________
AliMUON2DStoreValidator::~AliMUON2DStoreValidator()
{
  /// dtor
  delete fChambers;
  delete fStatus;
}

//_____________________________________________________________________________
AliMUONCheckItem* 
AliMUON2DStoreValidator::GetChamber(Int_t chamberID)
{
  /// Return (and create if not present) the given chamber
  /// chamberID in 0..NCh()
  
  if ( chamberID < 0 || chamberID >= AliMpConstants::NofTrackingChambers() )
  {
    AliFatal(Form("Invalid chamber number %d",chamberID));
    return 0x0;
  }
  
  if (!fChambers) 
  {
    fChambers = new TObjArray(AliMpConstants::NofTrackingChambers());
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
                              AliMpDDLStore::Instance()->GetDetElement(detElemId)->NofManus(),
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
    manu = new AliMUONCheckItem(manuId,AliMpDDLStore::Instance()->GetDetElement(detElemId)->NofChannelsInManu(manuId),"Manu");
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

  Int_t n(AliMpDDLStore::Instance()->GetDetElement(detElemId)->NofChannelsInManu(manuId));

  for ( Int_t i = 0; i < n; ++i )
  {
    AddMissingChannel(detElemId,manuId,i);
  }
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::ReportManu(TList& lines, const AliMUONCheckItem& manu)
{  
  /// Report list of missing channels from this manu
  
  TObjString* channel(0x0);
  TIter next(manu.CreateIterator());
  
  while ( ( channel = static_cast<TObjString*>(next()) ) )
  {
    lines.Add(new TObjString(Form("\t\t\tChannel %s is missing or dead",
                                  channel->GetString().Data())));
  }
  
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::ReportDE(TList& lines, const AliMUONCheckItem& de)
{  
  /// Report list of missing manus from this de
  AliMUONCheckItem* manu(0x0);
  
  TIter next(de.CreateIterator());
  
  lines.Add(new TObjString(Form("DE %5d",de.GetID())));
  
  
  while ( ( manu = static_cast<AliMUONCheckItem*>(next()) ) )
  {
    if ( manu->IsDead() )
    {
      lines.Add(new TObjString(Form("\t\tManu %4d is missing or dead",manu->GetID())));
    }
    else
    {
      ReportManu(lines,*manu);
    }
  }
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::ReportChamber(TList& lines, const AliMUONCheckItem& chamber)
{  
  /// Report list of missing de from this chamber
  
  AliMUONCheckItem* de(0x0);
  TIter next(chamber.CreateIterator());
  
  while ( ( de = static_cast<AliMUONCheckItem*>(next()) ) )
  {
    if ( de->IsDead() )
    {
      lines.Add(new TObjString(Form("\tDE %4d is missing or dead",de->GetID())));
    }
    else
    {
      ReportDE(lines,*de);
    }
  }
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::Report(TList& lines) const
{ 
  /// 
  if (fChambers) 
  {
    Report(lines,*fChambers); 
  }
}

//_____________________________________________________________________________
void
AliMUON2DStoreValidator::Report(TList& lines, const TObjArray& chambers)
{
  /// Reports what is missing, trying to be as concise as possible.
  
  for ( Int_t iChamber = 0; iChamber <= chambers.GetLast(); ++iChamber )
  {
    AliMUONCheckItem* chamber = static_cast<AliMUONCheckItem*>(chambers.At(iChamber));
    if ( chamber )
    {
      if ( chamber->IsDead() )
      {
        lines.Add(new TObjString(Form("Chamber %2d is missing or dead",iChamber)));
      }
      else
      {
        ReportChamber(lines,*chamber);
      }
    }
  }
}

//_____________________________________________________________________________
TObjArray* 
AliMUON2DStoreValidator::Validate(const AliMUONVStore& store,
                                  AliMUONVStore* config)
{                                  
  /// Validate the store. Check only the presence of all manus (i.e.
  /// check nothing about the values themselves). 
  /// Absence of manus which are not in the config is considered as normal.
  
  Bool_t (*kCheck)(const AliMUONVCalibParam&,Int_t) = 0x0;
  return Validate(store,kCheck,config);
}

//_____________________________________________________________________________
TObjArray* 
AliMUON2DStoreValidator::Validate(const AliMUONVStore& store,
                                  Bool_t (*check)(const AliMUONVCalibParam&,Int_t),
                                  AliMUONVStore* config)
{
  /// Validate the store. 
  /// The check method is used to decide if a store content value
  /// is valid or not.
  
  delete fChambers;
  fChambers = 0x0;
  
  // Now checks if some full manus are missing

  AliMpManuIterator it;

  Int_t detElemId;
  Int_t manuId;
  
  while ( it.Next(detElemId,manuId) )
  {
    AliMUONVCalibParam* test = 
      static_cast<AliMUONVCalibParam*>(store.FindObject(detElemId,manuId));
    if (!test)
    {
      // completely missing manu
      if ( !config || ( config && config->FindObject(detElemId,manuId ) ) )
      {
        // manu is in the config but not in the store : that's an error
        AddMissingManu(detElemId,manuId);
      }
    }
    else
    {
      if (!check) continue;
      
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
      
      // manu is there, check all its channels
      for ( Int_t manuChannel = 0 ; manuChannel < test->Size(); ++manuChannel )
      {
        if ( de->IsConnectedChannel(manuId,manuChannel) &&
             !check(*test,manuChannel) )             
        {
          AddMissingChannel(detElemId,manuId,manuChannel);
        }
      }
    }
  }
  return fChambers;
  
}


//_____________________________________________________________________________
TObjArray* 
AliMUON2DStoreValidator::Validate(const AliMUONVStore& store,
                                  Float_t invalidFloatValue,
                                  AliMUONVStore* config)
{
  /// Validate the store. 
  /// The invalidFloatValue is used to decide if a store content value
  /// is valid or not.
  
  delete fChambers;
  fChambers = 0x0;
  
  // Now checks if some full manus are missing

  AliMpManuIterator it;
  Int_t detElemId;
  Int_t manuId;
  
  while ( it.Next(detElemId,manuId) )
  {
    AliMUONVCalibParam* test = 
      static_cast<AliMUONVCalibParam*>(store.FindObject(detElemId,manuId));
    if (!test)
    {
      if ( !config || ( config && config->FindObject(detElemId,manuId ) ) )
      {
        // completely missing manu
        AddMissingManu(detElemId,manuId);
      }
    }
    else
    {
      // manu is there, check all its channels
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
      
      for ( Int_t manuChannel = 0 ; manuChannel < test->Size(); ++manuChannel )
      {
        if ( de->IsConnectedChannel(manuId,manuChannel) &&
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


