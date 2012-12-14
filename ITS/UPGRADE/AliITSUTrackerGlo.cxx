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
#include "AliITSUClusterPix.h"
using namespace AliITSUAux;
using namespace TMath;


//----------------- tmp stuff -----------------
Int_t currLabel=-1;


ClassImp(AliITSUTrackerGlo)
//_________________________________________________________________________
AliITSUTrackerGlo::AliITSUTrackerGlo(AliITSUReconstructor* rec)
:  fReconstructor(rec)
  ,fITS(0)
  ,fCurrESDtrack(0)
  ,fCurrMass(kPionMass)
  ,fSeedsLr(0)
  ,fSeedsPool("AliITSUSeed",0)
  ,fTrCond()
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
  delete[] fSeedsLr;
  //
}

//_________________________________________________________________________
void AliITSUTrackerGlo::Init(AliITSUReconstructor* rec)
{
  // init with external reconstructor
  //
  fITS = new AliITSURecoDet(rec->GetGeom(),"ITSURecoInterface");
  for (int ilr=fITS->GetNLayersActive();ilr--;) {
    fITS->GetLayerActive(ilr)->SetClusters(rec->GetClusters(ilr));
  }
  //
  fSeedsPool.ExpandCreateFast(1000); // RS TOCHECK
  int n = fITS->GetNLayersActive()+1;
  fSeedsLr = new TObjArray[n];
  //
  fTrCond.SetNLayers(fITS->GetNLayersActive());
  fTrCond.AddNewCondition(5);
  fTrCond.AddGroupPattern( (0x1<<0)|(0x1<<1) );
  fTrCond.AddGroupPattern( (0x1<<3)|(0x1<<4) );
  fTrCond.AddGroupPattern( (0x1<<5)|(0x1<<6) );
  //
  fTrCond.AddNewCondition(5);
  fTrCond.AddGroupPattern( (0x1<<0)|(0x1<<2) );
  fTrCond.AddGroupPattern( (0x1<<3)|(0x1<<4) );
  fTrCond.AddGroupPattern( (0x1<<5)|(0x1<<6) );
  //
  fTrCond.AddNewCondition(5);
  fTrCond.AddGroupPattern( (0x1<<1)|(0x1<<2) );
  fTrCond.AddGroupPattern( (0x1<<3)|(0x1<<4) );
  fTrCond.AddGroupPattern( (0x1<<5)|(0x1<<6) );
  //
  printf("Tracking Conditions: ");
  fTrCond.Print();
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::Clusters2Tracks(AliESDEvent *esdEv)
{
  //
  //
  AliITSUReconstructor::GetRecoParam()->Print();

  fITS->ProcessClusters();
  // select ESD tracks to propagate
  int nTrESD = esdEv->GetNumberOfTracks();
  for (int itr=0;itr<nTrESD;itr++) {
    AliESDtrack *esdTr = esdEv->GetTrack(itr);
    AliInfo(Form("Processing track %d | MCLabel: %d",itr,esdTr->GetTPCLabel()));
    currLabel = Abs(esdTr->GetTPCLabel());
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
  if (esdTr->Pt()<AliITSUReconstructor::GetRecoParam()->GetMinPtForProlongation()) return kFALSE;
  //
  float dtz[2];
  esdTr->GetDZ(GetX(),GetY(),GetZ(),bz,dtz); 
  // if track is not V0 candidata but has large offset wrt IP, reject it. RS TOCHECK
  if ( !(esdTr->GetV0Index(0)>0 && dtz[0]>AliITSUReconstructor::GetRecoParam()->GetMaxDforV0dghtrForProlongation()) 
       && (Abs(dtz[0])>AliITSUReconstructor::GetRecoParam()->GetMaxDForProlongation() ||
	   Abs(dtz[1])>AliITSUReconstructor::GetRecoParam()->GetMaxDZForProlongation())) return kFALSE;
  //
  return kTRUE;
}

//_________________________________________________________________________
void AliITSUTrackerGlo::FindTrack(AliESDtrack* esdTr)
{
  // find prolongaion candidates finding for single seed
  //
  if (!NeedToProlong(esdTr)) return;  // are we interested in this track?
  if (!InitSeed(esdTr))      return;  // initialize prolongations hypotheses tree
  //
  AliITSURecoSens *hitSens[AliITSURecoSens::kNNeighbors+1];
  AliITSUSeed seedUC;  // copy of the seed from the upper layer
  AliITSUSeed seedT;   // transient seed between the seedUC and new prolongation hypothesis
  //
  TObjArray clArr; // container for transfer of clusters matching to seed
  //
  for (int ila=fITS->GetNLayersActive();ila--;) {
    int ilaUp = ila+1;                         // prolong seeds from layer above
    int nSeedsUp = GetNSeeds(ilaUp);
    for (int isd=0;isd<nSeedsUp;isd++) {
      AliITSUSeed* seedU = GetSeed(ilaUp,isd);  // seed on prev.active layer to prolong
      seedUC = *seedU;                          // its copy will be prolonged
      seedUC.SetParent(seedU);
      // go till next active layer
      AliInfo(Form("working on Lr:%d Seed:%d of %d",ila,isd,nSeedsUp));
      if (!TransportToLayer(&seedUC, fITS->GetLrIDActive(ilaUp), fITS->GetLrIDActive(ila)) ) {
	//
	AliInfo("Transport failed");
	// Check if the seed satisfies to track definition
	if (NeedToKill(&seedUC,kTransportFailed)) seedU->Kill(); 
	continue; // RS TODO: decide what to do with tracks stopped on higher layers w/o killing
      }
      AliITSURecoLayer* lrA = fITS->GetLayerActive(ila);
      if (!GetRoadWidth(&seedUC, ila)) { // failed to find road width on the layer
	if (NeedToKill(&seedUC,kRWCheckFailed)) seedU->Kill(); 
	continue;
      }
      int nsens = lrA->FindSensors(&fTrImpData[kTrPhi0], hitSens);  // find detectors which may be hit by the track
      AliInfo(Form("Will check %d sensors on lr:%d ",nsens,ila));
      //
      for (int isn=nsens;isn--;) {
	seedT = seedUC;
	AliITSURecoSens* sens = hitSens[isn];
	//
	if (!seedT.Propagate(sens->GetPhiTF(),sens->GetXTF(),GetBz())) continue; // propagation failed, seedT is intact
	int clID0 = sens->GetFirstClusterId();
	for (int icl=sens->GetNClusters();icl--;) {
	  int res = CheckCluster(&seedT,ila,clID0+icl);
	  //
	  if (res==kStopSearchOnSensor) break;     // stop looking on this sensor
	  if (res==kClusterNotMatching) continue;  // cluster does not match
	  // cluster is matching and it was added to the hypotheses tree
	}
      }
      // cluster search is done. Do we need ta have a version of this seed skipping current layer
      seedT.SetLr(ila);
      if (!NeedToKill(&seedT,kMissingCluster)) AddProlongationHypothesis(NewSeedFromPool(&seedT) ,ila);      
    }
    printf(">>> All hypotheses on lr %d: \n",ila);
    for (int ih=0;ih<GetNSeeds(ila);ih++) {
      printf(" #%3d ",ih); GetSeed(ila,ih)->Print();
    }
  }
  //
  ResetSeedTree();
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::InitSeed(AliESDtrack *esdTr)
{
  // init prolongaion candidates finding for single seed
  fCurrMass = esdTr->GetMass();
  fCurrESDtrack = esdTr;
  if (fCurrMass<kPionMass*0.9) fCurrMass = kPionMass; // don't trust to mu, e identification from TPCin
  //
  //
  AliITSUSeed* seed = NewSeedFromPool();
  seed->SetLr(fITS->GetNLayersActive());   // fake layer
  seed->AliExternalTrackParam::operator=(*esdTr);
  seed->SetParent(esdTr);
  AddProlongationHypothesis(seed,fITS->GetNLayersActive());
  return kTRUE;
  // TO DO
}

//_________________________________________________________________________
void AliITSUTrackerGlo::ResetSeedTree()
{
  // reset current hypotheses tree
  for (int i=fITS->GetNLayersActive()+1;i--;) fSeedsLr[i].Clear();
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
	if (!PropagateSeed(seed,xToGo,fCurrMass, lrFr->GetMaxStep())) return kFALSE;
      }
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom+=dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    // go the entrance of the layer, assuming no materials in between
    double xToGo = dir>0 ? lrTo->GetRMin() : lrTo->GetRMax();
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) return kFALSE;
    lrFr = lrTo;
  }
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::GetRoadWidth(AliITSUSeed* seed, int ilrA)
{
  // calculate road width in terms of phi and z for the track which MUST be on the external radius of the layer
  // as well as some aux info
  double bz = GetBz();
  AliITSURecoLayer* lrA = fITS->GetLayerActive(ilrA);
  seed->GetXYZ(&fTrImpData[kTrXIn]);    // lab position at the entrance from above
  //
  fTrImpData[kTrPhiIn] = ATan2(fTrImpData[kTrYIn],fTrImpData[kTrXIn]);
  if (!seed->Rotate(fTrImpData[kTrPhiIn])) return kFALSE; // go to the frame of the entry point into the layer
  double dr  = lrA->GetDR();                              // approximate X dist at the inner radius
  if (!seed->GetXYZAt(seed->GetX()-dr, bz, fTrImpData + kTrXOut)) {
    // special case: track does not reach inner radius, might be tangential
    double r = seed->GetD(0,0,bz);
    double x;
    if (!seed->GetXatLabR(r,x,bz,-1)) {
      AliError(Form("This should not happen: r=%f",r));
      seed->Print();
      return kFALSE;
    }
    dr = Abs(seed->GetX() - x);
    if (!seed->GetXYZAt(x, bz, fTrImpData + kTrXOut)) {
      AliError(Form("This should not happen: x=%f",x));
      seed->Print();
      return kFALSE;      
    }
  }
  //
  fTrImpData[kTrPhiOut] = ATan2(fTrImpData[kTrYOut],fTrImpData[kTrXOut]);
  double sgy = seed->GetSigmaY2() + dr*dr*seed->GetSigmaSnp2() + AliITSUReconstructor::GetRecoParam()->GetSigmaY2(ilrA);
  double sgz = seed->GetSigmaZ2() + dr*dr*seed->GetSigmaTgl2() + AliITSUReconstructor::GetRecoParam()->GetSigmaZ2(ilrA);
  sgy = Sqrt(sgy)*AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadY();
  sgz = Sqrt(sgz)*AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadZ();
  fTrImpData[kTrPhi0] = 0.5*(fTrImpData[kTrPhiOut]+fTrImpData[kTrPhiIn]);
  fTrImpData[kTrZ0]   = 0.5*(fTrImpData[kTrZOut]+fTrImpData[kTrZIn]);
  fTrImpData[kTrDPhi] = 0.5*Abs(fTrImpData[kTrPhiOut]-fTrImpData[kTrPhiIn]) + sgy/lrA->GetR();
  fTrImpData[kTrDZ]   = 0.5*Abs(fTrImpData[kTrZOut]-fTrImpData[kTrZIn])   + sgz;
  //  
  return kTRUE;
}

//_________________________________________________________________________
AliITSUSeed* AliITSUTrackerGlo::NewSeedFromPool(const AliITSUSeed* src)
{
  // create new seed, optionally copying from the source
  return src ? 
    new(fSeedsPool[fSeedsPool.GetEntriesFast()]) AliITSUSeed(*src) :
    new(fSeedsPool[fSeedsPool.GetEntriesFast()]) AliITSUSeed();
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::CheckCluster(AliITSUSeed* track, Int_t lr, Int_t clID) 
{
  // Check if the cluster (in tracking frame!) is matching to track. 
  // The track must be already propagated to sensor tracking frame.
  // Returns:  kStopSearchOnSensor if the search on given sensor should be stopped, 
  //           kClusterMatching    if the cluster is matching
  //           kClusterMatching    otherwise
  //
  // The seed is already propagated to cluster
  const double kTolerX = 5e-4;
  AliCluster *cl = fITS->GetLayerActive(lr)->GetCluster(clID);
  //
  Bool_t goodCl = kFALSE;
  for (int i=0;i<3;i++) if (cl->GetLabel(i)>=0 && cl->GetLabel(i)==currLabel) goodCl = kTRUE;
  //
  if (TMath::Abs(cl->GetX())>kTolerX) { // if due to the misalingment X is large, propagate track only
    if (!track->PropagateParamOnlyTo(track->GetX()+cl->GetX(),GetBz())) {
      if (goodCl) {printf("Loose good cl: Failed propagation. |"); cl->Print();}
      return kStopSearchOnSensor; // propagation failed, seedT is intact
    }
  }
  double dy = cl->GetY()-track->GetY();
  double dz = cl->GetZ()-track->GetZ();
  //
  double dy2 = dy*dy;
  double tol2 = (track->GetSigmaY2() + AliITSUReconstructor::GetRecoParam()->GetSigmaY2(lr))*
    AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadY()*AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadY(); // RS TOOPTIMIZE
  if (dy2>tol2) {                          // the clusters are sorted in Z(col) then in Y(row). 
    if (goodCl) {printf("Loose good cl: dy2=%e > tol2=%e |",dy2,tol2); cl->Print();}
    if (dy>0) return kStopSearchOnSensor;  // No chance that other cluster of this sensor will match (all Y's will be even larger)
    else      return kClusterNotMatching;   // Other clusters may match
  }
  double dz2 = dz*dz;
  tol2 = (track->GetSigmaZ2() + AliITSUReconstructor::GetRecoParam()->GetSigmaZ2(lr))*
    AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadZ()*AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadZ(); // RS TOOPTIMIZE
  if (dz2>tol2) {
    if (goodCl) {printf("Loose good cl: dz2=%e > tol2=%e |",dz2,tol2); cl->Print();}
    return kClusterNotMatching; // Other clusters may match
  }
  //
  // check chi2
  Double_t p[2]={cl->GetY(), cl->GetZ()};
  Double_t cov[3]={cl->GetSigmaY2(), cl->GetSigmaYZ(), cl->GetSigmaZ2()};
  double chi2 = track->GetPredictedChi2(p,cov);
  if (chi2>AliITSUReconstructor::GetRecoParam()->GetMaxTr2ClChi2(lr)) {
    if (goodCl) {
      printf("Loose good cl: Chi2=%e > Chi2Max=%e \n",chi2,AliITSUReconstructor::GetRecoParam()->GetMaxTr2ClChi2(lr)); 
      track->Print("etp");
      cl->Print("");
    }
    return kClusterNotMatching;
  }
  //
  track = NewSeedFromPool(track);  // input track will be reused, use its clone for updates
  if (!track->Update(p,cov)) {
    if (goodCl) {printf("Loose good cl: Failed update |"); cl->Print();}
    return kClusterNotMatching;
  }
  track->SetChi2Cl(chi2);
  track->SetLrClusterID(lr,clID);
  cl->IncreaseClusterUsage();
  //
  AliInfo(Form("AddCl Cl%d lr:%d: dY:%+8.4f dZ:%+8.4f (MC: %5d %5d %5d) |Chi2=%f",clID,lr,dy,dz2,cl->GetLabel(0),cl->GetLabel(1),cl->GetLabel(2), chi2));
  //
  AddProlongationHypothesis(track,lr);
  //
  return kClusterMatching;
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::NeedToKill(AliITSUSeed *seed, Int_t flag)
{
  // check if the seed should not be discarded
  const UShort_t kMask = 0xffff;
  if (flag==kMissingCluster) {
    int lastChecked = seed->GetLayerID();
    UShort_t patt = seed->GetHitsPattern();
    if (lastChecked) patt |= ~(kMask<<lastChecked); // not all layers were checked, complete unchecked once by potential hits
    Bool_t seedOK = fTrCond.CheckPattern(patt);
    return !seedOK;
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliITSUTrackerGlo::PropagateSeed(AliITSUSeed *seed, Double_t xToGo, Double_t mass, Double_t maxStep, Bool_t matCorr) 
{
  // propagate seed to given x applying material correction if requested
  const Double_t kEpsilon = 1e-5;
  Double_t xpos     = seed->GetX();
  Int_t dir         = (xpos<xToGo) ? 1:-1;
  Double_t xyz0[3],xyz1[3],param[7];
  //
  if (matCorr) seed->GetXYZ(xyz1);   //starting global position
  while ( (xToGo-xpos)*dir > kEpsilon){
    Double_t step = dir*TMath::Min(TMath::Abs(xToGo-xpos), maxStep);
    Double_t x    = xpos+step;
    Double_t bz=GetBz();   // getting the local Bz
    if (!seed->PropagateToX(x,bz))  return kFALSE;
    if (matCorr) {
      xyz0[0]=xyz1[0]; // global pos at the beginning of step
      xyz0[1]=xyz1[1];
      xyz0[2]=xyz1[2];
      seed->GetXYZ(xyz1);    //  // global pos at the end of step
      MeanMaterialBudget(xyz0,xyz1,param);	
      Double_t xrho=param[0]*param[4], xx0=param[1];
      if (dir>0) xrho = -xrho; // outward should be negative
      if (!seed->ApplyMaterialCorrection(xx0,xrho,mass,kFALSE)) return kFALSE;
    }
    xpos = seed->GetX();
  }
  return kTRUE;
}
