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

ClassImp(AliITSUTrackerGlo)

const Double_t AliITSUTrackerGlo::fgkToler =  1e-6;// tolerance for layer on-surface check


//_________________________________________________________________________
AliITSUTrackerGlo::AliITSUTrackerGlo(AliITSUReconstructor* rec)
:  fReconstructor(rec)
  ,fITS(0)
  ,fCurrESDtrack(0)
  ,fCurrMass(kPionMass)
  ,fHypStore(100)
  ,fCurrHyp(0)
  ,fSeedsPool("AliITSUSeed",0)
  ,fFreeSeedsID(0)
  ,fNFreeSeeds(0)
  ,fLastSeedID(0)
  ,fTrCond()
  ,fTrackPhase(-1)
  ,fClInfo(0)
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
  delete[] fClInfo;
  //
}

//_________________________________________________________________________
void AliITSUTrackerGlo::Init(AliITSUReconstructor* rec)
{
  // init with external reconstructor
  //
  fITS = new AliITSURecoDet(rec->GetGeom(),"ITSURecoInterface");
  int nLr = fITS->GetNLayersActive();
  fClInfo = new Int_t[nLr<<1];
  for (int ilr=nLr;ilr--;) {
    fITS->GetLayerActive(ilr)->SetClusters(rec->GetClusters(ilr));
  }
  //
  fSeedsPool.ExpandCreateFast(1000); // RS TOCHECK
  fFreeSeedsID.Set(1000);
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
  SetTrackingPhase(kClus2Tracks);
  ResetSeedsPool();
  int nTrESD = esdEv->GetNumberOfTracks();
  AliInfo(Form("Will try to find prolongations for %d tracks",nTrESD));
  AliITSUReconstructor::GetRecoParam()->Print();
  fHypStore.Delete();
  if (fHypStore.GetSize()<nTrESD) fHypStore.Expand(nTrESD+100);
  //
  fITS->ProcessClusters();
  // select ESD tracks to propagate
  for (int itr=0;itr<nTrESD;itr++) {
    AliESDtrack *esdTr = esdEv->GetTrack(itr);
    AliInfo(Form("Processing track %d | MCLabel: %d",itr,esdTr->GetTPCLabel()));
    if (!NeedToProlong(esdTr)) continue;  // are we interested in this track?
    FindTrack(esdTr, itr);
  }
  //
  AliInfo(Form("SeedsPool: %d, BookedUpTo: %d, free: %d",fSeedsPool.GetSize(),fSeedsPool.GetEntriesFast(),fNFreeSeeds));
  fHypStore.Print();
  FinalizeHypotheses();
  //
  return 0;
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::PropagateBack(AliESDEvent *esdEv)
{
  //
  // Do outward fits in ITS
  //
  SetTrackingPhase(kPropBack);
  int nTrESD = esdEv->GetNumberOfTracks();
  AliInfo(Form("Will propagate back %d tracks",nTrESD));
  //
  double bz0 = GetBz();
  Double_t xyzTrk[3],xyzVtx[3]={GetX(),GetY(),GetZ()};
  AliITSUTrackHyp dummyTr,*currTr=0;
  const double kWatchStep=10.; // for larger steps watch arc vs segment difference
  Double_t times[AliPID::kSPECIES];
  //
  //
  for (int itr=0;itr<nTrESD;itr++) {
    AliESDtrack *esdTr = esdEv->GetTrack(itr);
    // Start time integral and add distance from current position to vertex 
    if (esdTr->IsOn(AliESDtrack::kITSout)) continue; // already done
    //
    esdTr->GetXYZ(xyzTrk); 
    Double_t dst = 0.;     // set initial track lenght, tof
    {
      double dxs = xyzTrk[0] - xyzVtx[0];
      double dys = xyzTrk[1] - xyzVtx[1];
      double dzs = xyzTrk[2] - xyzVtx[2];
      // RS: for large segment steps d use approximation of cicrular arc s by
      // s = 2R*asin(d/2R) = d/p asin(p) \approx d/p (p + 1/6 p^3) = d (1+1/6 p^2)
      // where R is the track radius, p = d/2R = 2C*d (C is the abs curvature)
      // Hence s^2/d^2 = (1+1/6 p^2)^2
      dst = dxs*dxs + dys*dys;
      if (dst > kWatchStep*kWatchStep) { // correct circular part for arc/segment factor
	double crv = Abs(esdTr->GetC(bz0));
	double fcarc = 1.+crv*crv*dst/6.;
	dst *= fcarc*fcarc;
      }
      dst += dzs*dzs;
      dst = Sqrt(dst); 
    }
    //    
    esdTr->SetStatus(AliESDtrack::kTIME);
    //
    if (!esdTr->IsOn(AliESDtrack::kITSin)) { // case of tracks w/o ITS prolongation: just set the time integration
      dummyTr.AliExternalTrackParam::operator=(*esdTr);
      dummyTr.StartTimeIntegral();
      dummyTr.AddTimeStep(dst);
      dummyTr.GetIntegratedTimes(times); 
      esdTr->SetIntegratedTimes(times);
      esdTr->SetIntegratedLength(dummyTr.GetIntegratedLength());
      continue;
    }
    //
    currTr = GetTrackHyp(itr);
    currTr->StartTimeIntegral();
    currTr->AddTimeStep(dst);
    printf("Before resetCov: "); currTr->AliExternalTrackParam::Print();
    currTr->ResetCovariance(10000);
    if (RefitTrack(currTr,fITS->GetRMax())) { // propagate to exit from the ITS/TPC screen
      UpdateESDTrack(currTr,AliESDtrack::kITSout);
    }
    else {
      AliInfo(Form("Refit Failed for track %d",itr));
    }
    //
  }
  //
  return 0;
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::RefitInward(AliESDEvent *esdEv)
{
  //
  // refit the tracks inward, using current cov.matrix
  //
  SetTrackingPhase(kRefitInw);
  Int_t nTrESD = esdEv->GetNumberOfTracks();
  AliInfo(Form("Will refit inward %d tracks",nTrESD));
  AliITSUTrackHyp *currTr=0;
  //
  for (int itr=0;itr<nTrESD;itr++) {
    AliESDtrack *esdTr = esdEv->GetTrack(itr);
    // Start time integral and add distance from current position to vertex 
    UInt_t trStat = esdTr->GetStatus();
    if ( !(trStat & AliESDtrack::kITSout) ) continue;
    if (   trStat & AliESDtrack::kITSrefit ) continue; // already done
    if (  (trStat & AliESDtrack::kTPCout) && !(trStat & AliESDtrack::kTPCrefit) ) continue;
    //
    currTr = GetTrackHyp(itr);
    currTr->AliExternalTrackParam::operator=(*esdTr);  // fetch current ESDtrack kinematics
    if (RefitTrack(currTr,fITS->GetRMin())) { // propagate up to inside radius of the beam pipe
      UpdateESDTrack(currTr,AliESDtrack::kITSrefit);
    }
    else {
      AliInfo(Form("Refit Failed for track %d",itr));
    }
  }    
  //
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
void AliITSUTrackerGlo::FindTrack(AliESDtrack* esdTr, Int_t esdID)
{
  // find prolongaion candidates finding for single seed
  //
  AliITSUSeed seedUC;  // copy of the seed from the upper layer
  AliITSUSeed seedT;   // transient seed between the seedUC and new prolongation hypothesis
  //
  if (!InitHypothesis(esdTr,esdID)) return;  // initialize prolongations hypotheses tree
  //
  AliITSURecoSens *hitSens[AliITSURecoSens::kNNeighbors+1];
  //
  TObjArray clArr; // container for transfer of clusters matching to seed
  //
  int nLrActive = fITS->GetNLayersActive();
  for (int ila=nLrActive;ila--;) {
    int ilaUp = ila+1;                         // prolong seeds from layer above
    //
    // for the outermost layer the seed is created from the ESD track
    int nSeedsUp = (ilaUp==nLrActive) ? 1 : fCurrHyp->GetNSeeds(ilaUp);
    //
    for (int isd=0;isd<nSeedsUp;isd++) {
      AliITSUSeed* seedU;
      if (ilaUp==nLrActive) {
	seedU = 0;
	seedUC.InitFromESDTrack(esdTr);
      }
      else {
	seedU = fCurrHyp->GetSeed(ilaUp,isd);  // seed on prev.active layer to prolong	
	seedUC = *seedU;                       // its copy will be prolonged
	seedUC.SetParent(seedU);	
      }
      seedUC.ResetFMatrix();                    // reset the matrix for propagation to next layer
      // go till next active layer
      AliInfo(Form("working on Lr:%d Seed:%d of %d",ila,isd,nSeedsUp));
      if (!TransportToLayer(&seedUC, fITS->GetLrIDActive(ilaUp), fITS->GetLrIDActive(ila)) ) {
	//
	AliInfo("Transport failed");
	// Check if the seed satisfies to track definition
	if (NeedToKill(&seedUC,kTransportFailed) && seedU) KillSeed(seedU,kTRUE);
	continue; // RS TODO: decide what to do with tracks stopped on higher layers w/o killing
      }
      AliITSURecoLayer* lrA = fITS->GetLayerActive(ila);
      if (!GetRoadWidth(&seedUC, ila)) { // failed to find road width on the layer
	if (NeedToKill(&seedUC,kRWCheckFailed) && seedU) KillSeed(seedU,kTRUE);
	continue;
      }
      int nsens = lrA->FindSensors(&fTrImpData[kTrPhi0], hitSens);  // find detectors which may be hit by the track
      AliInfo(Form("Will check %d sensors on lr:%d ",nsens,ila));
      //
      double bz = GetBz();
      for (int isn=nsens;isn--;) {
	seedT = seedUC;
	AliITSURecoSens* sens = hitSens[isn];
	//
	// We need to propagate the seed to sensor on lrA staying the frame of the sensor from prev.layer,
	// since the transport matrix should be defined in this frame.
	double xs; // X in the TF of current seed, corresponding to intersection with sensor plane
	if (!seedT.GetTrackingXAtXAlpha(sens->GetXTF(),sens->GetPhiTF(),bz, xs)) continue;
	if (!PropagateSeed(&seedT,xs,fCurrMass)) continue;
	//	if (!seedT.PropagateToX(xs,bz)) continue;
	//	if (!seedT.Rotate(sens->GetPhiTF())) continue;
	if (!seedT.RotateToAlpha(sens->GetPhiTF())) continue;
	//
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
      if (!NeedToKill(&seedT,kMissingCluster)) {
	AliITSUSeed* seedSkp = NewSeedFromPool(&seedT);
	double penalty = -AliITSUReconstructor::GetRecoParam()->GetMissPenalty(ila);
	// to do: make penalty to account for probability to miss the cluster for good reason
	seedSkp->SetChi2Cl(penalty);
	AddProlongationHypothesis(seedSkp,ila);      
      }
    }
    ((TObjArray*)fCurrHyp->GetLayerSeeds(ila))->Sort();
    printf(">>> All hypotheses on lr %d: \n",ila);
    for (int ih=0;ih<fCurrHyp->GetNSeeds(ila);ih++) {
      printf(" #%3d ",ih); fCurrHyp->GetSeed(ila,ih)->Print();
      //
      /*
      if (ila!=0) continue;
      double vecL[5] = {0};
      double matL[15] = {0};
      AliITSUSeed* sp = fCurrHyp->GetSeed(ila,ih);
      while(sp->GetParent()) {
	sp->Smooth(vecL,matL);
	if (sp->GetLayerID()>=fITS->GetNLayersActive()-1) break;
	sp = (AliITSUSeed*)sp->GetParent();
      }
      */
    }
  }
  //
  SaveCurrentTrackHypotheses();
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::InitHypothesis(AliESDtrack *esdTr, Int_t esdID)
{
  // init prolongaion candidates finding for single seed
  fCurrHyp = GetTrackHyp(esdID);
  if (fCurrHyp) return kTRUE;
  //
  fCurrMass = esdTr->GetMass();
  fCurrESDtrack = esdTr;
  if (fCurrMass<kPionMass*0.9) fCurrMass = kPionMass; // don't trust to mu, e identification from TPCin
  //
  fCurrHyp = new AliITSUTrackHyp(fITS->GetNLayersActive());
  fCurrHyp->SetESDTrack(esdTr);
  fCurrHyp->SetUniqueID(esdID);
  fCurrHyp->SetMass(fCurrMass);
  SetTrackHyp(fCurrHyp,esdID);
  //
  return kTRUE;
  // TO DO
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::TransportToLayer(AliITSUSeed* seed, Int_t lFrom, Int_t lTo)
{
  // transport seed from layerFrom to the entrance of layerTo
  //  
  //
  if (lTo==lFrom) AliFatal(Form("was called with lFrom=%d lTo=%d",lFrom,lTo));
  //
  int dir = lTo > lFrom ? 1:-1;
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom); // this can be 0 when extrapolation from TPC to ITS is requested
  Bool_t checkFirst = kTRUE;
  while(lFrom!=lTo) {
    if (lrFr) {
      if (!GoToExitFromLayer(seed,lrFr,dir,checkFirst)) return kFALSE; // go till the end of current layer
      checkFirst = kFALSE;
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom+=dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    // go the entrance of the layer, assuming no materials in between
    double xToGo = lrTo->GetR(-dir);
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) return kFALSE;
    lrFr = lrTo;
  }
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::TransportToLayer(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo)
{
  // transport track from layerFrom to the entrance of layerTo
  //  
  if (lTo==lFrom) AliFatal(Form("was called with lFrom=%d lTo=%d",lFrom,lTo));
  //
  int dir = lTo > lFrom ? 1:-1;
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom); // this can be 0 when extrapolation from TPC to ITS is requested
  Bool_t checkFirst = kTRUE;
  while(lFrom!=lTo) {
    if (lrFr) {
      if (!GoToExitFromLayer(seed,lrFr,dir,checkFirst)) return kFALSE; // go till the end of current layer
      checkFirst = kFALSE;
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom+=dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    if (!GoToExitFromLayer(seed,lrTo,dir)) return kFALSE; // go the entrance of the layer, assuming no materials in between
    lrFr = lrTo;
  }
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::GoToExitFromLayer(AliITSUSeed* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check)
{
  // go to the exit from lr in direction dir, applying material corrections in steps specific for this layer
  // If check is requested, do this only provided the track has not exited the layer already
  double xToGo = lr->GetR(dir);
  if (check) { // do we need to track till the surface of the current layer ?
    double curR2 = seed->GetX()*seed->GetX() + seed->GetY()*seed->GetY(); // current radius
    //    AliInfo(Form(" dir:%d Cur: %e Tgt: %e",dir,Sqrt(curR2),xToGo));
    if      (dir>0) { if (curR2-xToGo*xToGo>fgkToler) return kTRUE; } // on the surface or outside of the layer
    else if (dir<0) { if (xToGo*xToGo-curR2>fgkToler) return kTRUE; } // on the surface or outside of the layer
  }
  if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
  // go via layer to its boundary, applying material correction.
  if (!PropagateSeed(seed,xToGo,fCurrMass, lr->GetMaxStep())) return kFALSE;
  //
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::GoToExitFromLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check)
{
  // go to the exit from lr in direction dir, applying material corrections in steps specific for this layer
  // If check is requested, do this only provided the track has not exited the layer already
  double xToGo = lr->GetR(dir);
  if (check) { // do we need to track till the surface of the current layer ?
    double curR2 = seed->GetX()*seed->GetX() + seed->GetY()*seed->GetY(); // current radius
    if      (dir>0) { if (curR2-xToGo*xToGo>fgkToler) return kTRUE; } // on the surface or outside of the layer
    else if (dir<0) { if (xToGo*xToGo-curR2>fgkToler) return kTRUE; } // on the surface or outside of the layer
  }
  if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
  // go via layer to its boundary, applying material correction.
  if (!PropagateSeed(seed,xToGo,fCurrMass, lr->GetMaxStep())) return kFALSE;
  //
  return kTRUE;
  //
}


//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::GoToEntranceToLayer(AliITSUSeed* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check)
{
  // go to the entrance of lr in direction dir, w/o applying material corrections.
  // If check is requested, do this only provided the track did not reach the layer already
  double xToGo = lr->GetR(-dir);
  if (check) { // do we need to track till the surface of the current layer ?
    double curR2 = seed->GetX()*seed->GetX() + seed->GetY()*seed->GetY(); // current radius
    if      (dir>0) { if (curR2-xToGo*xToGo>fgkToler) return kTRUE; } // already passed
    else if (dir<0) { if (xToGo*xToGo-curR2>fgkToler) return kTRUE; } // already passed
  }
  if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
  // go via layer to its boundary, applying material correction.
  if (!PropagateSeed(seed,xToGo,fCurrMass, 100, kFALSE)) return kFALSE;
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::GoToEntranceToLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check)
{
  // go to the entrance of lr in direction dir, w/o applying material corrections.
  // If check is requested, do this only provided the track did not reach the layer already
  double xToGo = lr->GetR(-dir);
  if (check) { // do we need to track till the surface of the current layer ?
    double curR2 = seed->GetX()*seed->GetX() + seed->GetY()*seed->GetY(); // current radius
    if      (dir>0) { if (curR2-xToGo*xToGo>fgkToler) return kTRUE; } // already passed
    else if (dir<0) { if (xToGo*xToGo-curR2>fgkToler) return kTRUE; } // already passed
  }
  if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
  // go via layer to its boundary, applying material correction.
  if (!PropagateSeed(seed,xToGo,fCurrMass, 100, kFALSE)) return kFALSE;
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
  static AliExternalTrackParam sc;   // seed copy for manipalitions
  sc = *seed;
  //
  fTrImpData[kTrPhiIn] = ATan2(fTrImpData[kTrYIn],fTrImpData[kTrXIn]);
  if (!sc.Rotate(fTrImpData[kTrPhiIn])) return kFALSE; // go to the frame of the entry point into the layer
  double dr  = lrA->GetDR();                              // approximate X dist at the inner radius
  if (!sc.GetXYZAt(sc.GetX()-dr, bz, fTrImpData + kTrXOut)) {
    // special case: track does not reach inner radius, might be tangential
    double r = sc.GetD(0,0,bz);
    double x;
    if (!sc.GetXatLabR(r,x,bz,-1)) {
      sc.Print();
      AliFatal(Form("This should not happen: r=%f",r));
    }
    dr = Abs(sc.GetX() - x);
    if (!sc.GetXYZAt(x, bz, fTrImpData + kTrXOut)) {
      sc.Print();
      AliFatal(Form("This should not happen: x=%f",x));
    }
  }
  //
  fTrImpData[kTrPhiOut] = ATan2(fTrImpData[kTrYOut],fTrImpData[kTrXOut]);
  double sgy = sc.GetSigmaY2() + dr*dr*sc.GetSigmaSnp2() + AliITSUReconstructor::GetRecoParam()->GetSigmaY2(ilrA);
  double sgz = sc.GetSigmaZ2() + dr*dr*sc.GetSigmaTgl2() + AliITSUReconstructor::GetRecoParam()->GetSigmaZ2(ilrA);
  sgy = Sqrt(sgy)*AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadY();
  sgz = Sqrt(sgz)*AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadZ();
  fTrImpData[kTrPhi0] = 0.5*(fTrImpData[kTrPhiOut]+fTrImpData[kTrPhiIn]);
  fTrImpData[kTrZ0]   = 0.5*(fTrImpData[kTrZOut]+fTrImpData[kTrZIn]);
  fTrImpData[kTrDPhi] = 0.5*Abs(fTrImpData[kTrPhiOut]-fTrImpData[kTrPhiIn]) + sgy/lrA->GetR();
  fTrImpData[kTrDZ]   = 0.5*Abs(fTrImpData[kTrZOut]-fTrImpData[kTrZIn])   + sgz;
  //  
  return kTRUE;
}

//________________________________________
void AliITSUTrackerGlo::ResetSeedsPool()
{
  // mark all seeds in the pool as unused
  AliInfo(Form("CurrentSize: %d, BookedUpTo: %d, free: %d",fSeedsPool.GetSize(),fSeedsPool.GetEntriesFast(),fNFreeSeeds));
  fNFreeSeeds = 0;
  fSeedsPool.Clear(); // seeds don't allocate memory
}


//________________________________________
void AliITSUTrackerGlo::MarkSeedFree(AliITSUSeed *sd) 
{
  // account that this seed is "deleted" 
  int id = sd->GetPoolID();
  if (id<0) {
    AliError(Form("Freeing of seed %p NOT from the pool is requested",sd)); 
    return;
  }
  //  AliInfo(Form("%d %p",id, seed));
  fSeedsPool.RemoveAt(id);
  if (fFreeSeedsID.GetSize()<=fNFreeSeeds) fFreeSeedsID.Set( 2*fNFreeSeeds + 100 );
  fFreeSeedsID.GetArray()[fNFreeSeeds++] = id;
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
  const double kTolerX = 5e-4;
  AliCluster *cl = fITS->GetLayerActive(lr)->GetCluster(clID);
  //
  Bool_t goodCl = kFALSE;
  int currLabel = Abs(fCurrESDtrack->GetTPCLabel());
  //
  if (cl->GetLabel(0)>=0) {
    for (int i=0;i<3;i++) if (cl->GetLabel(i)>=0 && cl->GetLabel(i)==currLabel) {goodCl = kTRUE; break;}
  }
  else goodCl = kTRUE;
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
      printf("Loose good cl: Chi2=%e > Chi2Max=%e |dy: %+.3e dz: %+.3e\n",
	     chi2,AliITSUReconstructor::GetRecoParam()->GetMaxTr2ClChi2(lr),dy,dz); 
      track->Print("etp");
      cl->Print("");
    }
    return kClusterNotMatching;
  }
  //
  track = NewSeedFromPool(track);  // input track will be reused, use its clone for updates
  if (!track->Update()) {
    if (goodCl) {printf("Loose good cl: Failed update |"); cl->Print();}
    MarkSeedFree(track);
    return kClusterNotMatching;
  }
  track->SetChi2Cl(chi2);
  track->SetLrClusterID(lr,clID);
  cl->IncreaseClusterUsage();
  //
  track->SetFake(!goodCl);
  //
  AliInfo(Form("AddCl(%d) Cl%d lr:%d: dY:%+8.4f dZ:%+8.4f (MC: %5d %5d %5d) |Chi2=%f(%c)",
	       goodCl,clID,lr,dy,dz2,cl->GetLabel(0),cl->GetLabel(1),cl->GetLabel(2), chi2, track->IsFake() ? '-':'+'));
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
  Bool_t updTime = dir>0 && seed->IsStartedTimeIntegral();
  if (matCorr || updTime) seed->GetXYZ(xyz1);   //starting global position
  while ( (xToGo-xpos)*dir > kEpsilon){
    Double_t step = dir*TMath::Min(TMath::Abs(xToGo-xpos), maxStep);
    Double_t x    = xpos+step;
    Double_t bz=GetBz();   // getting the local Bz
    if (!seed->PropagateToX(x,bz))  return kFALSE;
    double ds = 0;
    if (matCorr || updTime) {
      xyz0[0]=xyz1[0]; // global pos at the beginning of step
      xyz0[1]=xyz1[1];
      xyz0[2]=xyz1[2];
      seed->GetXYZ(xyz1);    //  // global pos at the end of step
      if (matCorr) {
	MeanMaterialBudget(xyz0,xyz1,param);	
	Double_t xrho=param[0]*param[4], xx0=param[1];
	if (dir>0) xrho = -xrho; // outward should be negative
	if (!seed->ApplyMaterialCorrection(xx0,xrho,mass,kFALSE)) return kFALSE;
	ds = param[4];
      }
       else { // matCorr is not requested but time integral is
	double d0 = xyz1[0]-xyz0[0];
	double d1 = xyz1[1]-xyz0[1];
	double d2 = xyz1[2]-xyz0[2];	
	ds = TMath::Sqrt(d0*d0+d1*d1+d2*d2);
      }     
    }
    if (updTime) seed->AddTimeStep(ds);
    xpos = seed->GetX();
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliITSUTrackerGlo::PropagateSeed(AliExternalTrackParam *seed, Double_t xToGo, Double_t mass, Double_t maxStep, Bool_t matCorr) 
{
  // propagate seed to given x applying material correction if requested
  const Double_t kEpsilon = 1e-5;
  Double_t xpos     = seed->GetX();
  Int_t dir         = (xpos<xToGo) ? 1:-1;
  Double_t xyz0[3],xyz1[3],param[7];
  //
  Bool_t updTime = dir>0 && seed->IsStartedTimeIntegral();
  if (matCorr || updTime) seed->GetXYZ(xyz1);   //starting global position
  while ( (xToGo-xpos)*dir > kEpsilon){
    Double_t step = dir*TMath::Min(TMath::Abs(xToGo-xpos), maxStep);
    Double_t x    = xpos+step;
    Double_t bz=GetBz();   // getting the local Bz
    if (!seed->PropagateTo(x,bz))  return kFALSE;
    double ds = 0;
    if (matCorr || updTime) {
      xyz0[0]=xyz1[0]; // global pos at the beginning of step
      xyz0[1]=xyz1[1];
      xyz0[2]=xyz1[2];
      seed->GetXYZ(xyz1);    //  // global pos at the end of step
      //
      if (matCorr) {
	MeanMaterialBudget(xyz0,xyz1,param);	
	Double_t xrho=param[0]*param[4], xx0=param[1];
	if (dir>0) xrho = -xrho; // outward should be negative
	if (!seed->CorrectForMeanMaterial(xx0,xrho,mass)) return kFALSE;
	ds = param[4];
      }
      else { // matCorr is not requested but time integral is
	double d0 = xyz1[0]-xyz0[0];
	double d1 = xyz1[1]-xyz0[1];
	double d2 = xyz1[2]-xyz0[2];	
	ds = TMath::Sqrt(d0*d0+d1*d1+d2*d2);
      }
    }
    if (updTime) seed->AddTimeStep(ds);
    //
    xpos = seed->GetX();
  }
  return kTRUE;
}

//______________________________________________________________________________
void AliITSUTrackerGlo::SaveCurrentTrackHypotheses()
{
  // RS: shall we clean up killed seeds?
  fCurrHyp = 0;
  // TODO
  
}

//______________________________________________________________________________
void AliITSUTrackerGlo::FinalizeHypotheses()
{
  // select winner for each hypothesis, remove cl. sharing conflicts
  AliInfo("TODO");
  //
  int nh = fHypStore.GetEntriesFast();
  for (int ih=0;ih<nh;ih++) {
    AliITSUTrackHyp* hyp = (AliITSUTrackHyp*) fHypStore.UncheckedAt(ih); 
    if (!hyp || !hyp->DefineWinner()) continue; // TODO
    CookMCLabel(hyp);
    UpdateESDTrack(hyp,AliESDtrack::kITSin);
  }

}

//______________________________________________________________________________
void AliITSUTrackerGlo::UpdateESDTrack(AliITSUTrackHyp* hyp,Int_t flag)
{
  // update ESD track with current best hypothesis
  AliESDtrack* esdTr = hyp->GetESDTrack();
  if (!esdTr) return;
  AliITSUSeed* win = hyp->GetWinner();
  if (!win) return;
  switch (flag) {
  case AliESDtrack::kITSin: 
    esdTr->UpdateTrackParams(hyp,flag); // update kinematics
    // TODO: set cluster info
    break;
    //
  case AliESDtrack::kITSout: 
    esdTr->UpdateTrackParams(hyp,flag); // update kinematics
    // TODO: avoid setting friend
    break;
    //
  case AliESDtrack::kITSrefit: 
    esdTr->UpdateTrackParams(hyp,flag); // update kinematics
    // TODO: avoid setting cluster info
    break;
  default:
    AliFatal(Form("Unknown flag %d",flag));
  }
  //
  // transfer module indices
  // TODO
}

//______________________________________________________________________________
Bool_t AliITSUTrackerGlo::RefitTrack(AliITSUTrackHyp* trc, Double_t rDest, Bool_t stopAtLastCl)
{
  // refit track till radius rDest
  AliITSUTrackHyp tmpTr;
  //
  double rCurr = Sqrt(trc->GetX()*trc->GetX() + trc->GetY()*trc->GetY());
  int dir,lrStart,lrStop;
  //
  dir = rCurr<rDest ? 1 : -1;
  lrStart = fITS->FindFirstLayerID(rCurr,dir);
  lrStop  = fITS->FindLastLayerID(rDest,dir);
  printf("Refit %d: Lr:%d (%f) -> Lr:%d (%f)\n",dir,lrStart,rCurr, lrStop,rDest);
  printf("Before refit: "); trc->AliExternalTrackParam::Print();
  if (lrStop<0 || lrStart<0) AliFatal(Form("Failed to find start(%d) or last(%d) layers. Track from %.3f to %.3f",lrStart,lrStop,rCurr,rDest));
  //
  int ncl = trc->FetchClusterInfo(fClInfo);
  fCurrMass = trc->GetMass();
  tmpTr.AliKalmanTrack::operator=(*trc);
  tmpTr.SetChi2(0);
  int iclLr[2],nclLr,clCount=0;
  //
  for (int ilr=lrStart;ilr!=lrStop;ilr+=dir) {
    AliITSURecoLayer* lr = fITS->GetLayer(ilr);
    if ( dir*(rCurr-lr->GetR(dir))>0) continue; // this layer is already passed
    int ilrA2,ilrA = lr->GetActiveID();
    // passive layer or active w/o hits will be traversed on the way to next cluster
    if (!lr->IsActive() || fClInfo[ilrA2=(ilrA<<1)]<0) continue; 
    //
    if (ilr!=lrStart && !TransportToLayer(&tmpTr,lrStart,ilr)) {
      AliInfo(Form("Failed to transport %d -> %d\n",lrStart,ilr));
      return kFALSE; // go to the entrance to the layer
    }
    lrStart = ilr;
    //
    // select the order in which possible 2 clusters (in case of the overlap) will be traversed and fitted
    nclLr=0;
    if (dir>0) { // clusters are stored in increasing radius order
      iclLr[nclLr++]=fClInfo[ilrA2++];
      if (fClInfo[ilrA2]>=0) iclLr[nclLr++]=fClInfo[ilrA2];
    }
    else {
      if ( fClInfo[ilrA2+1]>=0 ) iclLr[nclLr++]=fClInfo[ilrA2+1];
      iclLr[nclLr++]=fClInfo[ilrA2];
    }
    //
    for (int icl=0;icl<nclLr;icl++) {
      AliITSUClusterPix* clus =  (AliITSUClusterPix*)lr->GetCluster(iclLr[icl]);
      AliITSURecoSens* sens = lr->GetSensorFromID(clus->GetVolumeId());
      if (!tmpTr.Rotate(sens->GetPhiTF())) {
	AliInfo(Form("Failed on rotate to %f",sens->GetPhiTF()));
	return kFALSE;
      }
      //printf("Refit cl:%d on lr %d Need to go %.4f -> %.4f\n",icl,ilrA,tmpTr.GetX(),sens->GetXTF()+clus->GetX());
      if (!PropagateSeed(&tmpTr,sens->GetXTF()+clus->GetX(),fCurrMass)) {
	AliInfo(Form("Failed on propagate to %f",sens->GetXTF()+clus->GetX()));	
	return kFALSE;
      }
      if (!tmpTr.Update(clus)) {
	AliInfo(Form("Failed on Update"));		
	return kFALSE;
      }
      //printf("AfterRefit: "); tmpTr.AliExternalTrackParam::Print();
      if (stopAtLastCl && ++clCount==ncl) return kTRUE; // it was requested to not propagate after last update
    }
    //
  }
  // All clusters were succesfully fitted. Even if the track does not reach rDest, this is enough to validate it.
  // Still, try to go as close as possible to rDest.
  //
  if (lrStart!=lrStop) {
    //printf("Going to last layer %d -> %d\n",lrStart,lrStop);
    if (!TransportToLayer(&tmpTr,lrStart,lrStop)) {
      AliInfo(Form("Failed on TransportToLayer %d->%d",lrStart,lrStop));		
      return kTRUE;
    }    
    if (!GoToExitFromLayer(&tmpTr,fITS->GetLayer(lrStop),dir)) {
      AliFatal(Form("Failed on GoToExitFromLayer %d",lrStop));		
      return kTRUE; // go till the exit from layer
    }
    //
    //printf("On exit from last layer\n");
    tmpTr.AliExternalTrackParam::Print();
    // go to the destination radius
    if (!tmpTr.GetXatLabR(rDest,rDest,GetBz(),dir)) return kTRUE;
    if (!PropagateSeed(&tmpTr,rDest,fCurrMass, 100, kFALSE)) return kTRUE;
  }
  trc->AliKalmanTrack::operator=(tmpTr);
  printf("After refit (now at lr %d): ",lrStart); trc->AliExternalTrackParam::Print();
  return kTRUE;
}

//__________________________________________________________________
void AliITSUTrackerGlo::CookMCLabel(AliITSUTrackHyp* hyp) 
{
  // build MC label
  //
  const int kMaxLbPerCl = 3;
  int lbID[kMaxLayers*6],lbStat[kMaxLayers*6];
  Int_t lr,nLab=0,nCl=0;
  AliITSUSeed *seed = hyp->GetWinner();
  while(seed) {
    int clID = seed->GetLrCluster(lr);
    if (clID>=0) {
      AliCluster *cl = fITS->GetLayerActive(lr)->GetCluster(clID);
      nCl++;
      for (int imc=0;imc<kMaxLbPerCl;imc++) { // labels within single cluster
	int trLb = cl->GetLabel(imc);
	if (imc<0) break;
	// search this mc track in already accounted ones
	int iLab;
	for (iLab=0;iLab<nLab;iLab++) if (lbID[iLab]==trLb) break;
	if (iLab<nLab) lbStat[iLab]++;
	else {
	  lbID[nLab] = trLb;
	  lbStat[nLab++] = 1;
	}
      } // loop over given cluster's labels
    }
    seed = (AliITSUSeed*)seed->GetParent();
  } // loop over clusters
  // 
  if (nCl) {
    int maxLab=0,nTPCok=0;
    AliESDtrack* esdTr = hyp->GetESDTrack();
    int tpcLab = esdTr ? Abs(esdTr->GetTPCLabel()) : -kDummyLabel;
    for (int ilb=nLab;ilb--;) {
      int st = lbStat[ilb];
      if (lbStat[maxLab]<st) maxLab=ilb;
      if (tpcLab==lbID[ilb]) nTPCok += st;
    }
    hyp->SetFakeRatio(1.-float(nTPCok)/nCl);
    hyp->SetLabel( nTPCok==nCl ? tpcLab : -tpcLab);
    hyp->SetITSLabel( lbStat[maxLab]==nCl ? lbStat[maxLab] : -lbStat[maxLab]); // winner label
    return;
  }
  //
  hyp->SetFakeRatio(-1.);
  hyp->SetLabel( kDummyLabel );
  hyp->SetITSLabel( kDummyLabel );
}
