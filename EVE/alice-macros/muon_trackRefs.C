// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file muon_trackRefs.C
/// \brief Macro to visualise trackRef in MUON spectrometer 
/// (both tracker and trigger).
///
/// Use muon_trackRefs(Bool_t showSimClusters) in order to run it
///
/// Needs that alieve_init() is already called
///
/// \author P. Pillot, L. Aphecetche; Subatech

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONClusterStoreV2.h"
#include "AliMUONRawClusterV2.h"
#include "AliMUONVCluster.h"
#include "AliMUONConstants.h"
#include "AliMUONRecoParam.h"
#include "AliMUONCDB.h"

#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliTrackReference.h"

#include "EveBase/AliEveMagField.h"
#include "EveBase/AliEveTrack.h"
#include "EveBase/AliEveEventManager.h"

#include <TEveManager.h>
#include <TEveUtil.h>
#include <TEveTrack.h>
#include <TEvePointSet.h>
#include <TEveVSDStructs.h>
#include <TEveTrackPropagator.h>

#include <TClonesArray.h>
#include <TTree.h>
#include <TParticle.h>
#include <TMath.h>
#include <TROOT.h>

#endif

//______________________________________________________________________________
void muon_trackRef_propagator_setup(TEveTrackPropagator* trkProp, Bool_t showVertex)
{
  // set magnetic field
  trkProp->SetMagFieldObj(new AliEveMagField);
  trkProp->SetStepper(TEveTrackPropagator::kRungeKutta);
  
  // set propagation range
  trkProp->SetMaxR(500.);
  trkProp->SetMaxZ(-AliMUONConstants::DefaultChamberZ(13)+3.5);
  
  // go through pathmarks
  trkProp->SetFitDaughters(kFALSE);
  trkProp->SetFitReferences(kTRUE);
  trkProp->SetFitDecay(kTRUE);
  trkProp->SetFitCluster2Ds(kFALSE);
  
  // Render first vertex if required
  if (showVertex)
  {
    trkProp->SetRnrFV(kTRUE);
    trkProp->RefFVAtt().SetMarkerSize(0.5);
    trkProp->RefFVAtt().SetMarkerColor(kRed);
  }
}

//______________________________________________________________________________
Bool_t isReconstructible(Bool_t* chHit)
{
  // load recoParam from OCDB
  static AliMUONRecoParam* gRecoParam = 0x0;
  static UInt_t gRequestedStationMask = 0;
  static Bool_t gRequest2ChInSameSt45 = kFALSE;
  if (!gRecoParam)
  {
    gRecoParam = AliMUONCDB::LoadRecoParam();
    if (!gRecoParam) exit(-1);
    // compute the mask of requested stations
    gRequestedStationMask = 0;
    for (Int_t i = 0; i < 5; i++) if (gRecoParam->RequestStation(i)) gRequestedStationMask |= ( 1 << i );
    // get whether a track need 2 chambers hit in the same station (4 or 5) or not to be reconstructible
    gRequest2ChInSameSt45 = !gRecoParam->MakeMoreTrackCandidates();
  }
  
  // check which chambers are hit
  UInt_t presentStationMask(0);
  Int_t nChHitInSt4 = 0, nChHitInSt5 = 0;
  for (Int_t ich=0; ich<10; ich++)
  {
    if (!chHit[ich]) continue;
    Int_t ist = ich/2;
    presentStationMask |= ( 1 << ist );
    if (ist == 3) nChHitInSt4++;
    if (ist == 4) nChHitInSt5++;
  }
  
  // at least one cluster per requested station
  if ((gRequestedStationMask & presentStationMask) != gRequestedStationMask) return kFALSE;
  
  // 2 chambers hit in the same station (4 or 5)
  if (gRequest2ChInSameSt45) return (nChHitInSt4 == 2 || nChHitInSt5 == 2);
  // or 2 chambers hit in station 4 & 5 together
  else return (nChHitInSt4+nChHitInSt5 >= 2);
}

//______________________________________________________________________________
void add_muon_trackRefs(AliStack* stack, TTree* treeTR, TEveTrackList* reco, TEveTrackList* other,
			TEvePointSet* RecoClusters, TEvePointSet* OtherClusters)
{
  TClonesArray* trackRefs = 0;
  treeTR->SetBranchAddress("TrackReferences", &trackRefs);
  Int_t nSimTracks = stack->GetNtrack();
  AliMUONClusterStoreV2 clusters;
  TClonesArray clustersTrigger("AliMUONRawClusterV2", 10);
  
  // loop over simulated track
  for (Int_t itr = 0; itr < nSimTracks; itr++)
  {
    treeTR->GetEntry(stack->TreeKEntry(itr));
    Int_t nTrackRefs = trackRefs->GetEntriesFast();
    if (nTrackRefs <= 0) continue;
    
    TEveTrack* track = 0x0;
    Bool_t chHit[10];
    for (Int_t ich=0; ich<10; ich++) chHit[ich] = kFALSE;
    
    // loop over simulated track hits
    for (Int_t itrR = 0; itrR < nTrackRefs; ++itrR)
    {
      AliTrackReference* atr = static_cast<AliTrackReference*>(trackRefs->UncheckedAt(itrR));
      
      // skip trackRefs not in MUON
      if (atr->DetectorId() != AliTrackReference::kMUON) continue;
      
      // record chamber hit
      Int_t detElemId = atr->UserId();
      Int_t chamberId = detElemId / 100 - 1;
      if (chamberId < 0) continue;
      if (chamberId < 10) chHit[chamberId] = kTRUE;
      
      // produce eve track if not already done
      if (!track)
      {
	TParticle* p = stack->Particle(itr);
	track = new AliEveTrack(p, itr, 0x0);
	track->SetName(Form("%s [%d]", p->GetName(), itr));
	track->SetStdTitle();
	track->SetSourceObject(p);
	
	clusters.Clear();
	clustersTrigger.Clear("C");
      }
      
      // add path mark
      track->AddPathMark(TEvePathMark(TEvePathMark::kReference,
				      TEveVector(atr->X(),  atr->Y(),  atr->Z()),
				      TEveVector(atr->Px(), atr->Py(), atr->Pz()),
				      atr->GetTime()));
      
      // produce clusters if required
      if (RecoClusters || OtherClusters)
      {
	// from tracker
	if (chamberId < 10)
	{
	  // produce a new cluster on that DE or update existing one
	  AliMUONVCluster* cl = 0x0;
	  Int_t clNum = -1;
	  do cl = clusters.FindObject(AliMUONVCluster::BuildUniqueID(chamberId, detElemId, ++clNum));
	  while (cl && cl->GetNDigits() == 2);
	  if (cl) {
	    cl->SetXYZ((cl->GetX() + atr->X()) / 2., (cl->GetY() + atr->Y()) / 2., (cl->GetZ() + atr->Z()) / 2.);
	  }
	  else
	  {
	    cl = clusters.Add(chamberId, detElemId, clNum);
	    cl->SetXYZ(atr->X(), atr->Y(), atr->Z());
	  }
	  cl->AddDigitId((UInt_t)itrR);
	}
	else // from trigger
	{
	  AliMUONVCluster* cl = new(clustersTrigger[clustersTrigger.GetLast()+1]) AliMUONRawClusterV2();
	  cl->SetXYZ(atr->X(), atr->Y(), atr->Z());
	}
      }
    }
    
    if (track)
    {
      track->SortPathMarksByTime();
      // stop track propagation at last path mark
      track->RefPathMarks().back().fType = TEvePathMarkT<double>::EType_e(TEvePathMark::kDecay);
      
      // add the track and trackRefs to proper lists
      if (isReconstructible(chHit)) {
	track->SetPropagator(reco->GetPropagator());
	track->SetAttLineAttMarker(reco);
	reco->AddElement(track);
	
	// trackRefs
	if (RecoClusters)
	{
	  // tracker
	  TIter next(clusters.CreateIterator());
	  gROOT->ProcessLine(Form("add_muon_clusters((TIter*)%p, (TEvePointSet*)%p);",&next, RecoClusters));
	  // trigger
	  TIter next2(clustersTrigger.MakeIterator());
	  gROOT->ProcessLine(Form("add_muon_clusters((TIter*)%p, (TEvePointSet*)%p);",&next2, RecoClusters));
	}
      }
      else {
	track->SetPropagator(other->GetPropagator());
	track->SetAttLineAttMarker(other);
	other->AddElement(track);
	
	// trackRefs
	if (OtherClusters)
	{
	  // tracker
	  TIter next(clusters.CreateIterator());
	  gROOT->ProcessLine(Form("add_muon_clusters((TIter*)%p, (TEvePointSet*)%p);",&next, OtherClusters));
	  // trigger
	  TIter next2(clustersTrigger.MakeIterator());
	  gROOT->ProcessLine(Form("add_muon_clusters((TIter*)%p, (TEvePointSet*)%p);",&next2, OtherClusters));
	}
      }
    }
  }
  delete trackRefs;
}

//______________________________________________________________________________
void muon_trackRefs(Bool_t showSimClusters)
{
  // load kinematics and trackRefs
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();
  if (!stack) return;
  rl->LoadTrackRefs();
  TTree* treeTR = rl->TreeTR();  
  if (!treeTR) return;
  
  // track containers
  TEveElementList* trackCont = new TEveElementList("Sim MUON Tracks");
  
  TEveTrackList* reco = new TEveTrackList("reconstructible");
  reco->SetRnrPoints(kFALSE);
  reco->SetRnrLine(kTRUE);
  reco->SetLineColor(kRed);
  reco->SetLineStyle(2);
  muon_trackRef_propagator_setup(reco->GetPropagator(), kTRUE);
  trackCont->AddElement(reco);
  
  TEveTrackList* other = new TEveTrackList("others");
  other->SetRnrPoints(kFALSE);
  other->SetRnrLine(kTRUE);
  other->SetLineColor(kWhite);
  other->SetLineStyle(3);
  muon_trackRef_propagator_setup(other->GetPropagator(), kFALSE);
  trackCont->AddElement(other);
  
  // cluster container
  TEveElementList* clusterCont = 0x0;
  TEvePointSet* RecoClusters = 0x0;
  TEvePointSet* OtherClusters = 0x0;
  if (showSimClusters)
  {
    clusterCont = new TEveElementList("Sim MUON Clusters");
    
    RecoClusters = new TEvePointSet(1000);
    RecoClusters->SetName("Reconstructibles");
    RecoClusters->SetPickable(kFALSE);
    RecoClusters->SetMarkerStyle(2);
    RecoClusters->SetMarkerColor(kRed);
    RecoClusters->SetMarkerSize(0.5);
    clusterCont->AddElement(RecoClusters);
    
    OtherClusters = new TEvePointSet(10000);
    OtherClusters->SetName("Others");
    OtherClusters->SetPickable(kFALSE);
    OtherClusters->SetMarkerStyle(2);
    OtherClusters->SetMarkerColor(kWhite);
    OtherClusters->SetMarkerSize(0.5);
    clusterCont->AddElement(OtherClusters);
    
    TEveUtil::LoadMacro("muon_clusters.C+");
  }
  
  // add tracks to the proper list and propagate them. Add also clusters if required.
  add_muon_trackRefs(stack, treeTR, reco, other, RecoClusters, OtherClusters);
  reco->MakeTracks();
  other->MakeTracks();
  trackCont->SetTitle(Form("N=%d", reco->NumChildren()+other->NumChildren()));
  reco->SetTitle(Form("N=%d",reco->NumChildren()));
  other->SetTitle(Form("N=%d",other->NumChildren()));
  if (showSimClusters)
  {
    clusterCont->SetTitle(Form("N=%d",RecoClusters->GetLastPoint()+OtherClusters->GetLastPoint()+2));
    RecoClusters->SetTitle(Form("N=%d",RecoClusters->GetLastPoint()+1));
    OtherClusters->SetTitle(Form("N=%d",OtherClusters->GetLastPoint()+1));
  }
  
  // add graphic containers
  gEve->DisableRedraw();
  gEve->AddElement(trackCont);
  if (clusterCont) gEve->AddElement(clusterCont);
  gEve->EnableRedraw();
  gEve->Redraw3D();
}
