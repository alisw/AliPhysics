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

/// \class AliMUONQAMappingCheck
/// 
/// Helper class for AliMUONQADataMakerRec, which takes care
/// of building an AliMUONVTrackerData to store the location of
/// all the clusters encountered during a run, and all the clusters
/// that have charge on only one cathode (aka mono-cathode clusters).
///
/// This is to easily spot mapping errors and/or region of complete
/// inefficiencies.
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONQAMappingCheck.h"

#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONTrackerData.h"
#include "AliMUONVCluster.h"
#include "AliMUONVDigit.h"
#include "AliMpConstants.h"
#include "AliMpDetElement.h"
#include "AliMpDDLStore.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"

#include "AliMpManuIterator.h"

/// \cond CLASSIMP
ClassImp(AliMUONQAMappingCheck)
/// \endcond

//_____________________________________________________________________________
AliMUONQAMappingCheck::AliMUONQAMappingCheck(Int_t runNumber)
: TObject(),
fStore(new AliMUON2DMap(kTRUE)),
fGeometryTransformer(new AliMUONGeometryTransformer),
fDigitCalibrator(new AliMUONDigitCalibrator(runNumber)),
fNumberOfEvents(0),
fNumberOfClusters(0),
fNumberOfMonoCathodeClusters(0),
fNumberOfLegitimateMonoCathodeClusters(0)
{
  /// Ctor
  
  AliCodeTimerAuto(Form("RUN %d",runNumber),0);
  
  fGeometryTransformer->LoadGeometryData();
  
  // Init the store with all the manus. Note that this is not strictly necessary, 
  // but it helps not to get its growth (that would otherwise happen in
  // AddClusterLocation each time we get a cluster associated with a manu where
  // we got no cluster yet) confused with a memory leak...
  AliMpManuIterator it;
  Int_t detElemId, manuId;
  
  while (it.Next(detElemId,manuId))
  {
    fStore->Add(new AliMUONCalibParamND(4,AliMpConstants::ManuNofChannels(),detElemId,manuId,0.0));
  }
}

//_____________________________________________________________________________
AliMUONQAMappingCheck::~AliMUONQAMappingCheck()
{
  /// Dtor. Report on the global number of clusters encountered
  AliInfo(Form("Nevents %d Nclusters %d Nmono-cathodes %d Nlegitimate-mono-cathodes %d",
               fNumberOfEvents,
               fNumberOfClusters,
               fNumberOfMonoCathodeClusters,
               fNumberOfLegitimateMonoCathodeClusters));
  delete fStore;
  delete fGeometryTransformer;
  delete fDigitCalibrator;
}

//____________________________________________________________________________
void AliMUONQAMappingCheck::AddClusterLocation(Int_t detElemId,
                                               Int_t manuId, Int_t manuChannel, 
                                               Bool_t monoCathode,
                                               Bool_t legitimateMonoCathode)
{
  /// Add one cluster location to our internal store
  if ( manuId > 0 )
  {
    AliMUONVCalibParam* p = static_cast<AliMUONVCalibParam*>(fStore->FindObject(detElemId,manuId));
    if (!p)
    {
      p = new AliMUONCalibParamND(4,AliMpConstants::ManuNofChannels(),detElemId,manuId,0.0);
      fStore->Add(p);
    }
    if ( !monoCathode) 
    {
      p->SetValueAsDouble(manuChannel,0,p->ValueAsDouble(manuChannel,0)+1.0);
    }
    else 
    {
      p->SetValueAsDouble(manuChannel,1,p->ValueAsDouble(manuChannel,1)+1.0); 
      if (!legitimateMonoCathode)
      {
        p->SetValueAsDouble(manuChannel,2,p->ValueAsDouble(manuChannel,2)+1.0); 
      }
    }
  }
}

//____________________________________________________________________________ 
void 
AliMUONQAMappingCheck::NewEvent()
{
  /// Increment number of events seen
  ++fNumberOfEvents;
}

//____________________________________________________________________________ 
AliMUONVTrackerData* 
AliMUONQAMappingCheck::CreateData(const char* name) const
{
  /// Make a trackerData from our internal store
  
  AliMUONVStore* store = static_cast<AliMUONVStore*>(fStore->Clone());
  
  TIter next(store->CreateIterator());
  AliMUONVCalibParam* param;
  
  while ( ( param = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    for ( Int_t i = 0; i < param->Size(); ++i ) 
    {
      param->SetValueAsDouble(i,3,fNumberOfEvents);
    }
  }
  
  AliMUONTrackerData* data = new AliMUONTrackerData(name,name,4,kTRUE);
  data->SetDimensionName(0,"all"); // all clusters
  data->SetDimensionName(1,"mono"); // mono-cathode clusters
  data->SetDimensionName(2,"suspect"); // not legitimate mono-cathode clusters
  data->SetDimensionName(3,"Nevents"); // number of events
  data->DisableChannelLevel();
  
  data->Add(*store);
  
  delete store;
  
  return data;
}

//____________________________________________________________________________ 
void 
AliMUONQAMappingCheck::GetClusterLocation(AliMUONVCluster& cluster, 
                                          Int_t& manuBending, Int_t& manuBendingChannel, 
                                          Int_t& manuNonBending, Int_t& manuNonBendingChannel,
                                          Bool_t& monoCathode, Bool_t& legitimateMonoCathode)
{
  /// Get the pad under the center of the cluster, and whether or not this cluster
  /// has charge on both cathodes
  
  Int_t detElemId = cluster.GetDetElemId();
  
  Double_t x,y,z;
  
  fGeometryTransformer->Global2Local(detElemId,cluster.GetX(),cluster.GetY(),cluster.GetZ(),x,y,z);
  
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
  
  const AliMpVSegmentation* segB = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,de->GetCathodType(AliMp::kBendingPlane));
  const AliMpVSegmentation* segNB = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,de->GetCathodType(AliMp::kNonBendingPlane));
  
  AliMpPad padB = segB->PadByPosition(x,y);
  AliMpPad padNB = segNB->PadByPosition(x,y);
  
  manuBending = padB.GetManuId();
  manuBendingChannel = padB.GetManuChannel();
  
  manuNonBending = padNB.GetManuId();
  manuNonBendingChannel = padNB.GetManuChannel();
  
  Bool_t bending(kFALSE);
  Bool_t nonBending(kFALSE);
  
  for ( Int_t i = 0; i < cluster.GetNDigits(); ++i ) 
//    for ( Int_t i = 0; i < cluster.GetNDigits() && !(bending && nonBending); ++i ) 
  {
    UInt_t digitId = cluster.GetDigitId(i);
    Int_t manuId = AliMUONVDigit::ManuId(digitId);
    if ( manuId > 0 )
    {
      if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) )
      {
        nonBending = kTRUE;
      }
      else
      {
        bending = kTRUE;
      }
    }
  }
  
  monoCathode = ( bending != nonBending );
  
  if ( monoCathode ) 
  {
    legitimateMonoCathode = kFALSE;
    if (!bending) 
    {
      if ( IsManuDead(detElemId,manuBending) ) legitimateMonoCathode = kTRUE;    
    }
    
    if (!nonBending) 
    {
    
      if ( IsManuDead(detElemId,manuNonBending) ) legitimateMonoCathode = kTRUE;
    }
  }
  
  if (!bending) manuBending *= -1;
  if (!nonBending) manuNonBending *= -1;
  
  ++fNumberOfClusters;
  if ( monoCathode ) 
  {
    ++fNumberOfMonoCathodeClusters;
    if ( legitimateMonoCathode ) ++fNumberOfLegitimateMonoCathodeClusters;  
  }
}

//____________________________________________________________________________ 
Bool_t AliMUONQAMappingCheck::IsManuDead(Int_t detElemId, Int_t manuId) const
{
  /// Using the statusmaker, tells if a given manu is to be considered 
  /// as dead (here dead means at least one manas dead)
  
  if ( manuId <= 0 ) return kTRUE;
  
  Int_t n(0);
  
  for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i) 
  {
    if ( IsChannelDead(detElemId,manuId,i) ) ++n;
  }
  return n > 16;
}

//____________________________________________________________________________ 
Bool_t AliMUONQAMappingCheck::IsChannelDead(Int_t detElemId, Int_t manuId, Int_t manuChannel) const
{
  /// Using the statusmaker, tells if a given channel is dead
  
  return ( fDigitCalibrator->StatusMap(detElemId,manuId,manuChannel) & (AliMUONPadStatusMapMaker::SelfDeadMask() != 0) );
}

//____________________________________________________________________________ 
void 
AliMUONQAMappingCheck::Store(AliMUONVCluster& cluster)
{
  /// Store information about a single cluster
                                      
  if ( cluster.GetCharge() < 10 ) return;
  
  Int_t manuBendingId, manuBendingChannel;
  Int_t manuNonBendingId, manuNonBendingChannel;
  Bool_t monoCathode, legitimateMonoCathode;
  
  GetClusterLocation(cluster, manuBendingId, manuBendingChannel,manuNonBendingId, manuNonBendingChannel, monoCathode,legitimateMonoCathode);
  
  AddClusterLocation(cluster.GetDetElemId(),manuBendingId,manuBendingChannel,monoCathode,legitimateMonoCathode);
  AddClusterLocation(cluster.GetDetElemId(),manuNonBendingId,manuNonBendingChannel,monoCathode,legitimateMonoCathode);    
}
