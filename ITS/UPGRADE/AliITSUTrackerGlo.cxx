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

/*
  TO CHECK
  GetBz usage

 */

//-------------------------------------------------------------------------
//  Implementation of the ITS Upgrade tracker mother class.              
//-------------------------------------------------------------------------
#include <TTree.h>
#include <Riostream.h> 
#include <TMath.h>

#include "AliITSUTrackerGlo.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliITSURecoDet.h"
#include "AliITSURecoSens.h"
#include "AliITSUReconstructor.h"
#include "AliITSReconstructor.h"
#include "AliITSUSeed.h"
#include "AliITSUAux.h"
using namespace AliITSUAux;
using namespace TMath;


ClassImp(AliITSUTrackerGlo)
//_________________________________________________________________________
AliITSUTrackerGlo::AliITSUTrackerGlo(AliITSUReconstructor* rec)
:  fReconstructor(rec)
  ,fITS(0)
  ,fCurrMass(kPionMass)
  ,fSeedsPool("AliITSUSeed",100)
{
  // Default constructor
  if (rec) Init(rec);
}

//_________________________________________________________________________
AliITSUTrackerGlo::~AliITSUTrackerGlo()
{
 // Default destructor
 //  
  delete fITS;
  //
}

//_________________________________________________________________________
void AliITSUTrackerGlo::Init(AliITSUReconstructor* rec)
{
  // init with external reconstructor
  fITS = new AliITSURecoDet(rec->GetGeom(),"ITSURecoInterface");
  for (int ilr=fITS->GetNLayersActive();ilr--;) {
    fITS->GetLayerActive(ilr)->SetClusters(rec->GetClusters(ilr));
  }
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::Clusters2Tracks(AliESDEvent *esdEv)
{
  //
  //
  fITS->ProcessClusters();
  // select ESD tracks to propagate
  int nTrESD = esdEv->GetNumberOfTracks();
  for (int itr=0;itr<nTrESD;itr++) {
    AliESDtrack *esdTr = esdEv->GetTrack(itr);
    FindTrack(esdTr);
  }

  return 0;
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::PropagateBack(AliESDEvent * /*event*/)
{
  //
  // To be implemented 
  //
  
 Info("PropagateBack","To be implemented");
  return 0;
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::RefitInward(AliESDEvent * /*event*/)
{
  //
  // To be implemented 
  //
  
  Info("RefitInward","To be implemented");
  return 0;
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::LoadClusters(TTree * treeRP)
{
  // read from tree (if pointer provided) or directly from the ITS reco interface
  //
  return fReconstructor->LoadClusters(treeRP);
} 

//_________________________________________________________________________
void AliITSUTrackerGlo::UnloadClusters()
{
  //
  // To be implemented 
  //
  
  Info("UnloadClusters","To be implemented");
} 
//_________________________________________________________________________
AliCluster * AliITSUTrackerGlo::GetCluster(Int_t /*index*/) const
{
  //
  // To be implemented 
  //
  
  Info("GetCluster","To be implemented");
  return 0x0;
} 

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::NeedToProlong(AliESDtrack* esdTr)
{
  // do we need to match this track to ITS?
  //
  static double bz = GetBz();
  if (!esdTr->IsOn(AliESDtrack::kTPCin) ||
      esdTr->IsOn(AliESDtrack::kTPCout) ||
      esdTr->IsOn(AliESDtrack::kITSin)  ||
      esdTr->GetKinkIndex(0)>0) return kFALSE;
  //
  if (esdTr->Pt()<AliITSReconstructor::GetRecoParam()->GetMinPtForProlongation()) return kFALSE;
  //
  float dtz[2];
  esdTr->GetDZ(GetX(),GetY(),GetZ(),bz,dtz); 
  // if track is not V0 candidata but has large offset wrt IP, reject it. RS TOCHECK
  if ( !(esdTr->GetV0Index(0)>0 && dtz[0]>AliITSReconstructor::GetRecoParam()->GetMaxDforV0dghtrForProlongation()) 
       && (Abs(dtz[0])>AliITSReconstructor::GetRecoParam()->GetMaxDForProlongation() ||
	   Abs(dtz[1])>AliITSReconstructor::GetRecoParam()->GetMaxDZForProlongation())) return kFALSE;
  //
  return kTRUE;
}

//_________________________________________________________________________
void AliITSUTrackerGlo::FindTrack(AliESDtrack* esdTr)
{
  // find prolongaion candidates finding for single seed
  //
  if (!NeedToProlong(esdTr)) return;  // are we interested in this track?
  if (!InitSeed(esdTr))         return;  // initialize prolongations hypotheses tree
  //
  double rwData[kNTrData];
  AliITSURecoSens *hitSens[AliITSURecoSens::kNNeighbors];
  AliITSUSeed seed0;
  //
  for (int ila=fITS->GetNLayersActive();ila--;) {
    int ilaUp = ila+1;                         // prolong seeds from layer above
    int nSeedsUp = GetNSeeds(ilaUp);
    for (int isd=0;isd<nSeedsUp;isd++) {
      seed0 = *GetSeed(ilaUp,isd);  // copy of the seed on prev.active layer
      // go till next active layer
      if (!TransportToLayer(&seed0, fITS->GetLrIDActive(ilaUp), fITS->GetLrIDActive(ila)) ) {
	//
	// Check if the seed satisfies to track definition
	if (NeedToKill(&seed0,kTransportFailed)) KillSeed(ilaUp,isd); 
	continue; // RS TODO: decide what to do with tracks stopped on higher layers w/o killing
      }
      AliITSURecoLayer* lrA = fITS->GetLayerActive(ila);
      if (!GetRoadWidth(&seed0, ila, rwData)) { // failed to find road width on the layer
	if (NeedToKill(&seed0,kRWCheckFailed)) KillSeed(ilaUp,isd); 
	continue;
      }
      int nsens = lrA->FindSensors(&rwData[kTrPhi0], hitSens);  // find detectors which may be hit by the track (max 4)
      //
      for (int isn=nsens;isn--;) {


      }
      
    }
  }
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::InitSeed(AliESDtrack *esdTr)
{
  // int prolongaion candidates finding for single seed
  fCurrMass = esdTr->GetMass();
  if (fCurrMass<kPionMass*0.9) fCurrMass = kPionMass; // don't trust to mu, e identification from TPCin
 
  return kTRUE;
  // TO DO
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::TransportToLayer(AliITSUSeed* seed, Int_t lFrom, Int_t lTo)
{
  // transport seed from layerFrom to the entrance of layerTo
  //  
  const double kToler = 1e-6; // tolerance for layer on-surface check
  //
  int dir = lTo > lFrom ? 1:-1;
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom); // this can be 0 when extrapolation from TPC to ITS is requested
  Bool_t checkFirst = kTRUE;
  while(lFrom!=lTo) {
    double curR2 = seed->GetX()*seed->GetX() + seed->GetY()*seed->GetY(); // current radius
    if (lrFr) {
      Bool_t doLayer = kTRUE;
      double xToGo = dir>0 ? lrFr->GetRMax() : lrFr->GetRMin();
      if (checkFirst) { // do we need to track till the surface of the current layer ?
	checkFirst = kFALSE;
	if      (dir>0) { if (curR2-xToGo*xToGo>kToler) doLayer = kFALSE; } // on the surface or outside of the layer
	else if (dir<0) { if (xToGo*xToGo-curR2>kToler) doLayer = kFALSE; } // on the surface or outside of the layer
      }
      if (doLayer) {
	if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
	// go via layer to its boundary, applying material correction.
	if (!PropagateTrackTo(seed,xToGo,fCurrMass, lrFr->GetMaxStep(), kFALSE, -1, 0, kTRUE)) return kFALSE;
      }
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom+=dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    // go the entrance of the layer, assuming no materials in between
    double xToGo = dir>0 ? lrTo->GetRMin() : lrTo->GetRMax();
    if (seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
    if (!seed->PropagateTo(xToGo, GetBz())) return kFALSE; // RS: do we need BxByBz?
    lrFr = lrTo;
  }
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::GetRoadWidth(AliITSUSeed* seed, int ilrA, double* rwData)
{
  // calculate road width in terms of phi and z for the track which MUST be on the external radius of the layer
  // as well as some aux info
  double bz = GetBz();
  AliITSURecoLayer* lrA = fITS->GetLayerActive(ilrA);
  seed->GetXYZ(&rwData[kTrXIn]);    // lab position at the entrance from above
  //
  rwData[kTrPhiIn] = ATan2(rwData[kTrYIn],rwData[kTrXIn]);
  if (!seed->Rotate(rwData[kTrPhiIn])) return kFALSE; // go to the frame of the entry point into the layer
  double dr  = lrA->GetDR();                          // approximate X dist at the inner radius
  if (!seed->GetXYZAt(seed->GetX()-dr, bz, rwData + kTrXOut)) {
    // special case: track does not reach inner radius, might be tangential
    double r = seed->GetD(0,0,bz);
    double x;
    if (!seed->GetXatLabR(r,x,bz,-1)) {
      AliError(Form("This should not happen: r=%f",r));
      seed->Print();
      return kFALSE;
    }
    dr = Abs(seed->GetX() - x);
    if (!seed->GetXYZAt(x, bz, rwData + kTrXOut)) {
      AliError(Form("This should not happen: x=%f",x));
      seed->Print();
      return kFALSE;      
    }
  }
  //
  rwData[kTrPhiOut] = ATan2(rwData[kTrYOut],rwData[kTrXOut]);
  double sgy = seed->GetSigmaY2() + dr*seed->GetSigmaSnp2() + AliITSReconstructor::GetRecoParam()->GetSigmaY2(ilrA);
  double sgz = seed->GetSigmaZ2() + dr*seed->GetSigmaTgl2() + AliITSReconstructor::GetRecoParam()->GetSigmaZ2(ilrA);
  sgy = Sqrt(sgy)*AliITSReconstructor::GetRecoParam()->GetNSigmaRoadY();
  sgz = Sqrt(sgz)*AliITSReconstructor::GetRecoParam()->GetNSigmaRoadZ();
  rwData[kTrPhi0] = 0.5*(rwData[kTrPhiOut]+rwData[kTrPhiIn]);
  rwData[kTrZ0]   = 0.5*(rwData[kTrZOut]+rwData[kTrPhiIn]);
  rwData[kTrDPhi] = 0.5*Abs(rwData[kTrPhiOut]-rwData[kTrPhiIn]) + sgy/lrA->GetR();
  rwData[kTrDZ]   = 0.5*Abs(rwData[kTrZOut]-rwData[kTrPhiIn])   + sgz;
  //  
  return kTRUE;
}

//_________________________________________________________________________
AliITSUSeed* AliITSUTrackerGlo::NewSeed(const AliITSUSeed* src)
{
  // create new seed, optionally copying from the source
  return src ? 
    new(fSeedsPool[fSeedsPool.GetEntriesFast()]) AliITSUSeed(*src) :
    new(fSeedsPool[fSeedsPool.GetEntriesFast()]) AliITSUSeed();
}
