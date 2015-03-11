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

#include "AliMUONSimpleClusterServer.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONCluster.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONPad.h"
#include "AliMUONTriggerTrackToTrackerClusters.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterFinder.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONRecoParam.h"
#include "AliMpArea.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpExMap.h"
#include "AliMpExMapIterator.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include <Riostream.h>
#include <TObjArray.h>
#include <TString.h>
#include <float.h>

/// \class AliMUONSimpleClusterServer
///
/// Implementation of AliMUONVClusterServer interface
/// 
/// 
/// \author Laurent Aphecetche, Subatech

using std::endl;
using std::cout;
/// \cond CLASSIMP  
ClassImp(AliMUONSimpleClusterServer)
/// \endcond

namespace
{
  TString AsString(const AliMpArea& area)
  {
    return Form("(X,Y)=(%7.3f,%7.3f) (DX,DY)=(%7.3f,%7.3f)",
                area.GetPositionX(),
                area.GetPositionY(),
                area.GetDimensionX(),  /// TBCL was Y !!!
                area.GetDimensionY());
  }
}

//_____________________________________________________________________________
AliMUONSimpleClusterServer::AliMUONSimpleClusterServer(AliMUONVClusterFinder* clusterFinder,
                                                       const AliMUONGeometryTransformer& transformer)
: AliMUONVClusterServer(),
  fDigitStore(0x0),
  fClusterFinder(clusterFinder),
  fkTransformer(transformer),
  fPads(),
  fTriggerTrackStore(0x0),
  fBypass(0x0)
{
    /// Ctor
    /// Note that we take ownership of the clusterFinder
    
    fPads[0] = new AliMpExMap;
    fPads[1] = new AliMpExMap;
    
    fPadsIterator[0] = fPads[0]->CreateIterator();
    fPadsIterator[1] = fPads[1]->CreateIterator();
}

//_____________________________________________________________________________
AliMUONSimpleClusterServer::~AliMUONSimpleClusterServer()
{
  /// Dtor
  delete fClusterFinder;
  delete fPads[0];
  delete fPads[1];
  delete fPadsIterator[0];
  delete fPadsIterator[1];
  delete fBypass;
}

//_____________________________________________________________________________
Int_t 
AliMUONSimpleClusterServer::Clusterize(Int_t chamberId,
                                       AliMUONVClusterStore& clusterStore,
                                       const AliMpArea& area,
                                       const AliMUONRecoParam* recoParam)
{
  /// Area is in absolute coordinate. If not valid, means clusterize all
  /// the chamber.
  ///
  /// We first find out the list of DE that have a non-zero overlap with area,
  /// and then use the clusterfinder to find clusters in those areas (and DE).
  
  AliCodeTimerAuto(Form("Chamber %d",chamberId),0);
  
  if ( fTriggerTrackStore && chamberId >= 6 ) 
  {
    return fBypass->GenerateClusters(chamberId,clusterStore);
  }
  
  if (!recoParam) {
    AliError("Reconstruction parameters are missing: unable to clusterize");
    return 0;
  }
  
  AliMpDEIterator it;
  
  it.First(chamberId);
  
  Int_t nofAddedClusters(0);
  Int_t fNCluster = clusterStore.GetSize();

  AliDebug(1,Form("chamberId = %2d NofClusters before = %d searchArea=%s",
                  chamberId,fNCluster,AsString(area).Data()));
  
  while ( !it.IsDone() )
  {
    Int_t detElemId = it.CurrentDEId();
    
    TObjArray* pads[2] = 
    { 
      static_cast<TObjArray*>(fPads[0]->GetValue(detElemId)),
      static_cast<TObjArray*>(fPads[1]->GetValue(detElemId)) 
    };
    
    if ( ( pads[0] && pads[0]->GetLast()>=0 ) || 
         ( pads[1] && pads[1]->GetLast()>=0 ) )
    {
      AliMpArea deArea; // area in DE-local-coordinates
      Bool_t ok(kTRUE);
      
      if ( area.IsValid() ) 
      {
        ok = Overlap(detElemId,area,deArea);
      }
      
      if ( ok ) 
      {      
	const AliMpVSegmentation* seg[2] = 
        { AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::kCath0),
          AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::kCath1)
        };
        
        fClusterFinder->SetChargeHints(recoParam->LowestPadCharge(),
                                       recoParam->LowestClusterCharge());
        
        if ( fClusterFinder->NeedSegmentation() )
        {
          fClusterFinder->Prepare(detElemId,pads,deArea,seg);
        }
        else
        {
          fClusterFinder->Prepare(detElemId,pads,deArea);
        }
        
        AliDebug(1,Form("Clusterizing DE %04d with %3d pads (cath0) and %3d pads (cath1)",
                        detElemId,
                        (pads[0] ? pads[0]->GetLast()+1 : 0),
                        (pads[1] ? pads[1]->GetLast()+1 : 0)));
        
        AliMUONCluster* cluster;
        
        while ( ( cluster = fClusterFinder->NextCluster() ) ) 
        {      
          // add new cluster to the store with information to build its ID
          // increment the number of clusters into the store
          AliMUONVCluster* rawCluster = clusterStore.Add(chamberId, detElemId, fNCluster++);
          
          ++nofAddedClusters;
          
          // fill array of Id of digits attached to this cluster
          Int_t nPad = cluster->Multiplicity();
          if (nPad < 1) AliWarning("no pad attached to the cluster");
          
          for (Int_t iPad=0; iPad<nPad; iPad++) 
          {
            AliMUONPad *pad = cluster->Pad(iPad);
	    
	    // skip virtual pads
	    if (!pad->IsReal()) continue;
            
	    rawCluster->AddDigitId(pad->GetUniqueID());
          }
          
          // fill charge and other cluster informations
          rawCluster->SetCharge(cluster->Charge());
          rawCluster->SetChi2(cluster->Chi2());
          
          Double_t xg, yg, zg;
          fkTransformer.Local2Global(detElemId, 
                                    cluster->Position().X(), cluster->Position().Y(), 
                                    0, xg, yg, zg);
          rawCluster->SetXYZ(xg, yg, zg);
          rawCluster->SetErrXY(recoParam->GetDefaultNonBendingReso(chamberId),recoParam->GetDefaultBendingReso(chamberId));
          
          // Set MC label
          if (fDigitStore && fDigitStore->HasMCInformation()) 
          {
            rawCluster->SetMCLabel(FindMCLabel(*cluster, detElemId, seg));
          }
          
          AliDebug(1,Form("Adding RawCluster detElemId %4d mult %2d charge %e (xl,yl,zl)=(%e,%e,%e) (xg,yg,zg)=(%e,%e,%e) label %d",
                          detElemId,rawCluster->GetNDigits(),rawCluster->GetCharge(),
                          cluster->Position().X(),cluster->Position().Y(),0.0,
                          xg,yg,zg,rawCluster->GetMCLabel()));
        }
      }
    }
    it.Next();
  }
  
  AliDebug(1,Form("chamberId = %2d NofClusters after = %d",chamberId,fNCluster));
  
  return nofAddedClusters;
}


//_____________________________________________________________________________
void
AliMUONSimpleClusterServer::Global2Local(Int_t detElemId, const AliMpArea& globalArea,
                                         AliMpArea& localArea) const
{
  /// Convert a global area in local area for a given DE
  
  Double_t xl,yl,zl;
  
  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
  if ( chamberId < 0 ) {
    AliErrorStream() << "Cannot get chamberId from detElemId=" << detElemId << endl;
    return;
  }  
  Double_t zg = AliMUONConstants::DefaultChamberZ(chamberId);
  
  fkTransformer.Global2Local(detElemId,
                             globalArea.GetPositionX(),globalArea.GetPositionY(),zg,
                             xl,yl,zl);
  
  localArea = AliMpArea(xl,yl, globalArea.GetDimensionX(), globalArea.GetDimensionY());
}

//_____________________________________________________________________________
Bool_t
AliMUONSimpleClusterServer::Overlap(Int_t detElemId,
                                    const AliMpArea& area,
                                    AliMpArea& deArea) const
{
  /// Check whether (global) area overlaps with the given DE.
  /// If it is, set deArea to the overlap region and convert it
  /// in the local coordinate system of that DE.
  
  Bool_t overlap(kFALSE);
  
  AliMpArea* globalDEArea = fkTransformer.GetDEArea(detElemId);
  
  if (!globalDEArea) return kFALSE;
  
  AliMpArea overlapArea;
  
  if ( area.Overlap(*globalDEArea) )
  {
    overlapArea = area.Intersect(*globalDEArea);
    Global2Local(detElemId,overlapArea,deArea);
    overlap = kTRUE;
  }
  else
  {
    deArea = AliMpArea();
  }
  
  AliDebug(1,Form("DE %04d area %s globalDEArea %s overlapArea %s deArea %s overlap=%d",
                  detElemId,
                  AsString(area).Data(),
                  AsString(*globalDEArea).Data(),
                  AsString(overlapArea).Data(),
                  AsString(deArea).Data(),
                  overlap));
                  
  return overlap;
}

//_____________________________________________________________________________
TObjArray* 
AliMUONSimpleClusterServer::PadArray(Int_t detElemId, Int_t cathode) const
{
  /// Return array for given cathode of given DE
  
  return static_cast<TObjArray*>(fPads[cathode]->GetValue(detElemId));
}

//_____________________________________________________________________________
Bool_t 
AliMUONSimpleClusterServer::UseTriggerTrackStore(AliMUONVTriggerTrackStore* trackStore)
{
  /// Tells us to use trigger track store, and thus to bypass St45 clusters
  fTriggerTrackStore = trackStore; // not owner
  delete fBypass;
  fBypass = new AliMUONTriggerTrackToTrackerClusters(fkTransformer,fTriggerTrackStore);
  return kTRUE;
}

//_____________________________________________________________________________
void 
AliMUONSimpleClusterServer::UseDigits(TIter& next, AliMUONVDigitStore* digitStore)
{
  /// Convert digitStore into two arrays of AliMUONPads
  
  fDigitStore = digitStore;
  
  // Clear pads arrays in the maps
  for ( Int_t i=0; i<2; i++ ) {
    fPadsIterator[i]->Reset();
    Int_t key; TObject* obj;
    while ( ( obj = fPadsIterator[i]->Next(key) ) ) {
      //cout << "clearing array for detElemId " << key <<  "  ";
      obj->Clear();
    }  
  }   

  AliMUONVDigit* d;
  while ( ( d = static_cast<AliMUONVDigit*>(next()) ) )
  {
    if ( ! (d->Charge() > 0.) ) continue; // skip void digits.
    if ( ! d->IsTracker() ) continue; // skip trigger digits
    Int_t ix = d->PadX();
    Int_t iy = d->PadY();
    Int_t cathode = d->Cathode();
    Int_t detElemId = d->DetElemId();
    const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->
      GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
    AliMpPad pad = seg->PadByIndices(ix,iy);
    
    TObjArray* padArray = PadArray(detElemId,cathode);
    if (!padArray)
    {
      padArray = new TObjArray(100);
      padArray->SetOwner(kTRUE);
      fPads[cathode]->Add(detElemId,padArray);
    }
    
    AliMUONPad* mpad = new AliMUONPad(detElemId,cathode,
                    ix,iy,pad.GetPositionX(),pad.GetPositionY(),
                    pad.GetDimensionX(),pad.GetDimensionY(),
                    d->Charge());
    if ( d->IsSaturated() ) mpad->SetSaturated(kTRUE);
    mpad->SetUniqueID(d->GetUniqueID());
    padArray->Add(mpad);      
  }
}

//_____________________________________________________________________________
Int_t
AliMUONSimpleClusterServer::FindMCLabel(const AliMUONCluster& cluster, Int_t detElemId, const AliMpVSegmentation* seg[2]) const
{
  /// Find the label of the most contributing MC track (-1 in case of failure)
  /// The data member fDigitStore must be set
  
  // --- get the digit (if any) located at the cluster position on both cathods ---
  Int_t nTracks[2] = {0, 0};
  AliMUONVDigit* digit[2] = {0x0, 0x0};
  for (Int_t iCath = 0; iCath < 2; iCath++) {
    AliMpPad pad 
      = seg[AliMp::GetCathodType(iCath)]->PadByPosition(cluster.Position().X(), cluster.Position().Y(),kFALSE);
    if (pad.IsValid()) {
      digit[iCath] = fDigitStore->FindObject(detElemId, pad.GetManuId(), pad.GetManuChannel(), iCath);
      if (digit[iCath]) nTracks[iCath] = digit[iCath]->Ntracks();
    }
  }
  
  if (nTracks[0] + nTracks[1] == 0) return -1;
  
  // --- build the list of contributing tracks and of the associated charge ---
  Int_t* trackId = new Int_t[nTracks[0] + nTracks[1]];
  Float_t* trackCharge = new Float_t[nTracks[0] + nTracks[1]];
  Int_t nTracksTot = 0;
  
  // fill with contributing tracks on first cathod
  for (Int_t iTrack1 = 0; iTrack1 < nTracks[0]; iTrack1++) {
    trackId[iTrack1] = digit[0]->Track(iTrack1);
    trackCharge[iTrack1] = digit[0]->TrackCharge(iTrack1);
  }
  nTracksTot = nTracks[0];
  
  // complement with contributing tracks on second cathod
  for (Int_t iTrack2 = 0; iTrack2 < nTracks[1]; iTrack2++) {
    Int_t trackId2 = digit[1]->Track(iTrack2);
    // check if track exist
    Bool_t trackExist = kFALSE;
    for (Int_t iTrack1 = 0; iTrack1 < nTracks[0]; iTrack1++) {
      if (trackId2 == trackId[iTrack1]) {
	// complement existing track
	trackCharge[iTrack1] += digit[1]->TrackCharge(iTrack2);
	trackExist = kTRUE;
	break;
      }
    }
    // add the new track
    if (!trackExist) {
      trackId[nTracksTot] = trackId2;
      trackCharge[nTracksTot] = digit[1]->TrackCharge(iTrack2);
      nTracksTot++;
    }
  }
  
  // --- Find the most contributing track ---
  Int_t mainTrackId = -1;
  Float_t maxCharge = 0.;
  for (Int_t iTrack = 0; iTrack < nTracksTot; iTrack++) {
    if (trackCharge[iTrack] > maxCharge) {
      mainTrackId = trackId[iTrack];
      maxCharge = trackCharge[iTrack];
    }
  }
  
  delete[] trackId;
  delete[] trackCharge;
  
  return mainTrackId;
}

//_____________________________________________________________________________
void 
AliMUONSimpleClusterServer::Print(Option_t*) const
{
  /// Printout for debug only
  
  AliMpDEIterator it;
  
  it.First();
  
  while ( !it.IsDone() )
  {
    Int_t detElemId = it.CurrentDEId();
    
    // printout the number of pads / de, and number of used pads / de
    
    if ( ( PadArray(detElemId,0) && PadArray(detElemId,0)->GetLast() >= 0 ) || 
         ( PadArray(detElemId,1) && PadArray(detElemId,1)->GetLast() >= 0 ) )
    {
      cout << Form("---- DE %04d",detElemId) << endl;
      
      for ( Int_t cathode = 0; cathode < 2; ++cathode ) 
      {
        cout << Form("  -- Cathode %1d",cathode) << endl;
        
        TObjArray* padArray = PadArray(detElemId,cathode);
        
        if (!padArray)
        {
          cout << "no pad array" << endl;
        }
        else if ( padArray->GetLast() < 0 ) 
        {
          cout << "no pads" << endl;
        }
        else
        {
          TIter next(padArray);
          AliMUONPad* pad;
          while ( ( pad = static_cast<AliMUONPad*>(next()) ) )
          {
            pad->Print("full");
          }
        }
      }
    }
    it.Next();
  }
}  


