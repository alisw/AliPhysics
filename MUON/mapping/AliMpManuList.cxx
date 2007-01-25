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

#include "AliMpManuList.h"

#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpCathodType.h"
#include "AliMpVSegmentation.h"
#include "TArrayI.h"
#include "TList.h"

///
/// \class AliMpManuList
///
/// A sort of cache for mapping information we use often (or that are
/// time consuming to recompute).
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMpManuList)
/// \endcond

//_____________________________________________________________________________
AliMpManuList::AliMpManuList()
{
  /// ctor
}

//_____________________________________________________________________________
AliMpManuList::~AliMpManuList()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t 
AliMpManuList::DoesChannelExist(Int_t detElemId, Int_t manuID, Int_t manuChannel)
{
  /// Whether a given (detElemId,manuID,manuChannel) combination is a valid one
  
  const AliMpVSegmentation* seg = 
    AliMpSegmentation::Instance()
      ->GetMpSegmentationByElectronics(detElemId,manuID);
  if (!seg) return kFALSE;
  
  if ( seg->PadByLocation(AliMpIntPair(manuID,manuChannel),kFALSE).IsValid() )
  {
    return kTRUE;
  }
  else
  {
    return kFALSE;
  }
}

//_____________________________________________________________________________
TList*
AliMpManuList::ManuList()
{
  /// Create a TList of AliMpIntPair<detElemId,manuID> of all MUON TRK manus
  /// The returned list must be deleted by the client
  
  TList* manuList = new TList;
  
  manuList->SetOwner(kTRUE);
  
  AliMpDEIterator it;
  
  it.First();
  
  while ( !it.IsDone() )
  {
    Int_t detElemId = it.CurrentDEId();
    AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
    if ( stationType != AliMp::kStationTrigger ) 
    {
      for ( Int_t cath = AliMp::kCath0; cath <=AliMp::kCath1 ; ++cath )
      {
        const AliMpVSegmentation* seg 
	  = AliMpSegmentation::Instance()
            ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));
        
        TArrayI manus;
        
        seg->GetAllElectronicCardIDs(manus);
        
        for ( Int_t im = 0; im < manus.GetSize(); ++im )
        {
          manuList->Add(new AliMpIntPair(detElemId,manus[im]));
        }        
      }
    }
    it.Next();
  }
  return manuList;
}

//_____________________________________________________________________________
Int_t 
AliMpManuList::NumberOfChannels(Int_t detElemId, Int_t manuId)
{
  /// Returns the number of channels in that manuID. Answer should be <=64
  /// whatever happens.
  
  const AliMpVSegmentation* seg = 
    AliMpSegmentation::Instance()
      ->GetMpSegmentationByElectronics(detElemId,manuId);
  Int_t n(0);
  for ( Int_t i = 0; i < 64; ++i )
  {
    AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,i),kFALSE);
    if (pad.IsValid()) ++n;
  }
  return n;
}

//_____________________________________________________________________________
Int_t 
AliMpManuList::NumberOfManus(Int_t detElemId)
{
  /// Returns the number of manus contained in the given detection element.
  Int_t n(0);
  for ( Int_t i = AliMp::kCath0; i <= AliMp::kCath1; ++i )
  {
    const AliMpVSegmentation* seg 
      = AliMpSegmentation::Instance()
        ->GetMpSegmentation(detElemId,AliMp::GetCathodType(i));
        
    TArrayI manus;
    seg->GetAllElectronicCardIDs(manus);
    n += manus.GetSize();
  }
  return n;
}

