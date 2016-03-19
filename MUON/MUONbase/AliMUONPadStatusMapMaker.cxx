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
/// \class AliMUONPadStatusMapMaker
/// 
/// Convert a pad statuses into pad status maps.
/// 
/// A pad status is one 32-bits word describing whether this pad pedestal, lv, 
/// hv is correct or not.
///
/// A pad status *map* is one 32-bits (of which 24 only are used)
/// word describing whether this pad neighbours are ok or not
/// (whether a pad is ok or not is determined by applying a given
/// bitmask to the pad status word). Each bit in this word is related to one
/// neighbour, assuming the pad itself is at bit 0
///
/// ----------------
/// |  3 |  5 |  8 |
/// ----------------
/// |  2 |  0 |  7 |
/// ----------------
/// |  1 |  4 |  6 |
/// ----------------
///
/// Note that for instance in NonBending plane of slats, at the boundaries
/// between two pad densities, the pictures is a bit different, e.g.
/// (bits in () are always zero)
///
/// so some care must be taken when designing a mask to be tested ;-) if you
/// want to go farther than immediate neighbours...
///
/// If a pad is at a physical boundary, is will for sure have some bits at 1
/// (i.e. a non-existing neighbour is considered = bad).
///
///
/// add something about the reject list/probabilities here... (LA)
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONPadStatusMapMaker.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONRejectList.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"
#include "AliMpConstants.h"
#include <Riostream.h>
#include <TList.h>
#include "TRandom.h"
#include <cassert>

/// \cond CLASSIMP
ClassImp(AliMUONPadStatusMapMaker)
/// \endcond

Int_t AliMUONPadStatusMapMaker::fgkSelfDead = 1;

//_____________________________________________________________________________
AliMUONPadStatusMapMaker::AliMUONPadStatusMapMaker(const AliMUONPadStatusMaker& padStatusMaker,
                                                   Int_t mask,
                                                   Bool_t deferredInitialization) 
: TObject(),
fkStatusMaker(padStatusMaker),
fMask(mask),
fStatusMap(new AliMUON2DMap(true)),
fRejectProbabilities(new AliMUON2DMap(true)),
fRejectList(0x0),
fComputeOnDemand(deferredInitialization)
{
  /// ctor
  if (!deferredInitialization)
  {
    AliCodeTimerAuto("Computing complete status map at once",0);
    AliMUONVStore* neighboursStore = padStatusMaker.NeighboursStore();
    AliMUONVCalibParam* param;
    TIter next(neighboursStore->CreateIterator());
    while ( ( param = static_cast<AliMUONVCalibParam*>(next()) ) )
    {
      Int_t detElemId = param->ID0();
      Int_t manuId = param->ID1();
      ComputeStatusMap(detElemId,manuId);
    }
  }
  
  /// Whatever the deferred flag is, we *have* to compute the reject 
  /// probabilities here and now, for *all* channels.
  
  AliMUONRejectList* rl = padStatusMaker.CalibrationData().RejectList();
  
  if (rl)
  {
    AliMpManuIterator it;
    Int_t detElemId;
    Int_t manuId;
    
    while ( it.Next(detElemId,manuId) )
    {
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
      Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
      
      AliMUONVCalibParam* param = new AliMUONCalibParamNF(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,0);
      
      Int_t n(0);
      
      for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i ) 
      {
        Float_t proba(0.0);
        
        if ( de->IsConnectedChannel(manuId,i) )
        {
          proba = TMath::Max(rl->DetectionElementProbability(detElemId),rl->BusPatchProbability(busPatchId));
          
          proba = TMath::Max(proba,rl->ManuProbability(detElemId,manuId));
          
          proba = TMath::Max(proba,rl->ChannelProbability(detElemId,manuId,i));
          
          if ( proba > 0 ) 
          {
            ++n;
            param->SetValueAsFloat(i,0,proba);
          }
        }
      }
      
      if ( n > 0 ) 
      {
        fRejectProbabilities->Add(param);
      }
      else
      {
        // no need to add empty stuff...
        delete param;
      }
    }
  
    if ( rl->IsBinary())
    {
      fRejectList = fRejectProbabilities;
      fRejectProbabilities = 0x0;
      AliDebug(1,"RejectList = RejectProbabilities");
      StdoutToAliDebug(1,fRejectList->Print("","MEAN"));
    }
    else
    {
      AliWarning("Will run with non trivial survival probabilities for channels, manus, etc... Better check this is a simulation and not real data !");
      fRejectList = new AliMUON2DMap(true);
    }
  }
  else
  {
    fRejectList = fRejectProbabilities;
    fRejectProbabilities = 0x0;
    AliInfo("No RejectList found, so no RejectList will be used.");
  }
}

//_____________________________________________________________________________
AliMUONPadStatusMapMaker::~AliMUONPadStatusMapMaker()
{
  /// dtor
  delete fStatusMap;
  delete fRejectProbabilities;
  delete fRejectList;
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONPadStatusMapMaker::ComputeStatusMap(Int_t detElemId, Int_t manuId) const
{
  /// Compute the status map for a given manu, and add it to our internal
  /// fStatusMap internal storage
  
  AliCodeTimerAuto("(Int_t,Int_t)",0)
    
  AliMUONVCalibParam* param = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),
                                                      detElemId,manuId,-1);    
                                    
  Bool_t ok = fStatusMap->Add(param);
  if (!ok)
  {
    AliFatal(Form("Could not add manu %d of de %d",manuId,detElemId));
  }
                                  
  AliMUONVCalibParam* neighbours = fkStatusMaker.Neighbours(detElemId,manuId);
  
  AliMUONVCalibParam* statusParam = fkStatusMaker.PadStatus(detElemId,manuId);
  
  Int_t n = neighbours->Dimension();
  
  for ( Int_t manuChannel = 0; manuChannel < param->Size(); ++manuChannel )
  {
    Int_t statusMap(0);
    
    Int_t x = neighbours->ValueAsIntFast(manuChannel,0);
    if ( x < 0 ) 
    {
      // channel is not a valid one (i.e. (manuId,manuChannel) is not an existing pad)
      statusMap = -1;//fgkSelfDead;
      continue;
    }
        
    for ( Int_t i = 0; i < n; ++i )
    {
      // Compute the statusmap related to the status of neighbouring
      // pads. An invalid pad means "outside of edges".
            
      Int_t y = neighbours->ValueAsIntFast(manuChannel,i);      
      Int_t m,c;
      neighbours->UnpackValue(y,m,c);
      if ( c < 0 ) continue;
      Int_t status = 0;
      if ( !m )
      {
        status = -1;
      }
      else
      {
        status = statusParam->ValueAsIntFast(c); //fkStatusMaker.PadStatus(detElemId,m,c);
      }
      if ( ( fMask != 0 ) && ( (status & fMask) != 0 ) )
      {
        statusMap |= (1<<i);
      }
    }    
    param->SetValueAsIntFast(manuChannel,0,statusMap);
  }
  return param;
}

//_____________________________________________________________________________
void 
AliMUONPadStatusMapMaker::RefreshRejectProbabilities()
{
  /// From the (fixed) fRejectProbabilities, compute
  /// a fRejectList that will be valid for one event
  /// If fRejectProbabilities=0x0 it means we're dealing with
  /// trivial probabilities (0 or 1) and those are assumed to be already
  /// in fRejectList then.
  
  if ( !fRejectProbabilities ) return;
  
  AliCodeTimerAuto("",0);
  
  fRejectList->Clear();
  
  TIter next(fRejectProbabilities->CreateIterator());
  AliMUONVCalibParam* paramProba;
  AliMUONVCalibParam* paramReject;
  
  while ( ( paramProba = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    paramReject = new AliMUONCalibParamNF(1,paramProba->Size(),paramProba->ID0(),paramProba->ID1(),0.0);
    
    Int_t n(0);
    
    for ( Int_t i = 0; i < paramProba->Size(); ++i ) 
    {
      Float_t proba = paramProba->ValueAsFloat(i);
      Float_t x(proba);
      
      if ( proba > 0.0 && proba < 1.0 ) 
      {
        x = gRandom->Rndm();
        proba = ( x < proba ) ? 1.0 : 0.0;
      }
      
      if (proba>0.0)
      {
        ++n;
        paramReject->SetValueAsFloat(i,0,proba);
      }
    }
    if (n) fRejectList->Add(paramReject);
  }
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMapMaker::StatusMap(Int_t detElemId, Int_t manuId, 
                                    Int_t manuChannel) const
                                      
{
  /// Get the pad status map
  
  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fStatusMap->FindObject(detElemId,manuId));
  if (!param)
  {
    if ( fComputeOnDemand ) 
    {
      // not yet computed, so do it now
      param = ComputeStatusMap(detElemId,manuId);
    }
    else
    {
      // we're locked. probably a bad manuId ?
      return fgkSelfDead;
    }
  }
  
  Int_t statusMap = param->ValueAsInt(manuChannel);
  
  AliMUONVCalibParam* r = static_cast<AliMUONVCalibParam*>(fRejectList->FindObject(detElemId,manuId));
  
  if (r)
  {
    Float_t v= r->ValueAsFloat(manuChannel);
    
    assert (v==0.0 || v==1.0 ); 

    if ( v > 0 ) 
    {
      statusMap |= fgkSelfDead;
    }
  }
  
  return statusMap;
  
}
