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
#include <TFile.h>

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
  ,fCurrESDtrMClb(kDummyLabel)
  ,fCurrMass(kPionMass)
  ,fCountProlongationTrials(0)
  ,fCountITSin(0)
  ,fCountITSout(0)
  ,fCountITSrefit(0)
  ,fHypStore(100)
  ,fNBranchesAdded(0)
  ,fNCandidatesAdded(0)
  ,fCurrHyp(0)
  ,fSeedsPool("AliITSUSeed",0)
  ,fFreeSeedsID(0)
  ,fNFreeSeeds(0)
  ,fLastSeedID(0)
  ,fDefTrackConds(0)
  ,fCurrTrackCond(0)
  ,fTrackPhase(-1)
  ,fClInfo(0)
#ifdef  _FILL_CONTROL_HISTOS_
  ,fCHistoArr(0)
#endif
{
  // Default constructor
  if (rec) Init(rec);
}

//_________________________________________________________________________
AliITSUTrackerGlo::~AliITSUTrackerGlo()
{
 // Default destructor
 //  
  delete[] fClInfo;
  //
#ifdef  _FILL_CONTROL_HISTOS_
  if (fCHistoArr) {
    TFile* ctrOut = TFile::Open("itsuControlHistos.root","recreate");
    ctrOut->cd();
    AliInfo("Storing control histos");
    fCHistoArr->Print();
    //    ctrOut->WriteObject(fCHistoArr,"controlH","kSingleKey");
    fCHistoArr->Write();
    ctrOut->Close();
    delete ctrOut;
    fCHistoArr = 0;
  }
#endif 
  //
}

//_________________________________________________________________________
void AliITSUTrackerGlo::Init(AliITSUReconstructor* rec)
{
  // init with external reconstructor
  //
  fITS = rec->GetITSInterface();
  int nLr = fITS->GetNLayersActive();
  fClInfo = new Int_t[nLr<<1];
  //
  fSeedsPool.ExpandCreateFast(1000); // RS TOCHECK
  fFreeSeedsID.Set(1000);
  //
}

//_________________________________________________________________________
void AliITSUTrackerGlo::CreateDefaultTrackCond()
{
  // creates default tracking conditions to be used when recoparam does not provide them
  int nLr = fITS->GetNLayersActive();
  fClInfo = new Int_t[nLr<<1];
  //
  AliITSUTrackCond* cond = new AliITSUTrackCond();
  //
  cond->SetNLayers(fITS->GetNLayersActive());
  cond->AddNewCondition(nLr);
  cond->AddGroupPattern( 0xffff ); // require all layers hit
  //
  fDefTrackConds.AddLast(cond);
  //
  AliInfo("Created conditions: ");
  for (int i=0;i<fDefTrackConds.GetEntriesFast();i++) fDefTrackConds[i]->Print();
  //
}


//_________________________________________________________________________
Int_t AliITSUTrackerGlo::Clusters2Tracks(AliESDEvent *esdEv)
{
  //
  SetTrackingPhase(kClus2Tracks);
  //
#ifdef  _FILL_CONTROL_HISTOS_
  if (!fCHistoArr) BookControlHistos();
#endif

  TObjArray *trackConds = 0;
  //
  fCountProlongationTrials = 0;
  fCountITSin = 0;
  fCountITSout = 0;
  fCountITSrefit = 0;
  //
  ResetSeedsPool();
  int nTrESD = esdEv->GetNumberOfTracks();
  AliInfo(Form("Will try to find prolongations for %d tracks",nTrESD));
  int nTrackCond = AliITSUReconstructor::GetRecoParam()->GetNTrackingConditions();
  if (nTrackCond<1) {
    if (!fDefTrackConds.GetEntriesFast()) {
      AliInfo("No tracking conditions found in recoparams, creating default one requesting all layers hit");
      CreateDefaultTrackCond();
    }
    trackConds = &fDefTrackConds;
    nTrackCond = trackConds->GetEntriesFast();
  } 
  else trackConds = AliITSUReconstructor::GetRecoParam()->GetTrackingConditions();
  //
  static Bool_t first = kTRUE;
  if (first) {
    AliITSUReconstructor::GetRecoParam()->Print();
    first = kFALSE;
  }
  fHypStore.Delete();
  if (fHypStore.GetSize()<nTrESD) fHypStore.Expand(nTrESD+100);
  //
  fITS->ProcessClusters();
  //
  FlagSplitClusters(); // tmp RS
  //
  for (int icnd=0;icnd<nTrackCond;icnd++) {
    fCurrTrackCond = (AliITSUTrackCond*)trackConds->UncheckedAt(icnd);
    // select ESD tracks to propagate
    for (int itr=0;itr<nTrESD;itr++) {
      fCurrESDtrack = esdEv->GetTrack(itr);
      fCurrESDtrMClb = fCurrESDtrack->GetLabel();
      //
      if (!NeedToProlong(fCurrESDtrack)) continue;  // are we interested in this track?
      AliDebug(+1,Form("Processing track %d | M=%.3f Pt=%.3f | MCLabel: %d",itr,fCurrESDtrack->GetMass(kTRUE),fCurrESDtrack->Pt(),fCurrESDtrMClb));//RS
      FindTrack(fCurrESDtrack, itr);
    }   
    //
    if (AliDebugLevelClass()>+2) {
      AliInfo(Form("SeedsPool: %d, BookedUpTo: %d, free: %d",fSeedsPool.GetSize(),fSeedsPool.GetEntriesFast(),fNFreeSeeds));
      fHypStore.Print();
    }
    FinalizeHypotheses();
  }
  //
  AliInfo(Form("%d ITSin for %d tried TPC seeds out of %d ESD tracks\n",fCountITSin,fCountProlongationTrials,nTrESD));
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
  AliDebug(1,Form("Will propagate back %d tracks",nTrESD));
  //
  double bz0 = GetBz();
  Double_t xyzTrk[3],xyzVtx[3]={GetX(),GetY(),GetZ()};
  AliITSUTrackHyp dummyTr,*currTr=0;
  const double kWatchStep=10.; // for larger steps watch arc vs segment difference
  Double_t times[AliPID::kSPECIES];
  //
  //
  for (int itr=0;itr<nTrESD;itr++) {
    fCurrESDtrack = esdEv->GetTrack(itr);
    fCurrESDtrMClb = fCurrESDtrack->GetLabel();
    // Start time integral and add distance from current position to vertex 
    if (fCurrESDtrack->IsOn(AliESDtrack::kITSout)) continue; // already done
    //
    fCurrESDtrack->GetXYZ(xyzTrk); 
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
	double crv = Abs(fCurrESDtrack->GetC(bz0));
	double fcarc = 1.+crv*crv*dst/6.;
	dst *= fcarc*fcarc;
      }
      dst += dzs*dzs;
      dst = Sqrt(dst); 
    }
    //    
    fCurrESDtrack->SetStatus(AliESDtrack::kTIME);
    //
    if (!fCurrESDtrack->IsOn(AliESDtrack::kITSin)) { // case of tracks w/o ITS prolongation: just set the time integration
      dummyTr.AliExternalTrackParam::operator=(*fCurrESDtrack);
      dummyTr.StartTimeIntegral();
      dummyTr.AddTimeStep(dst);
      dummyTr.GetIntegratedTimes(times); 
      fCurrESDtrack->SetIntegratedTimes(times);
      fCurrESDtrack->SetIntegratedLength(dummyTr.GetIntegratedLength());
      continue;
    }
    //
    currTr = GetTrackHyp(itr);
    currTr->StartTimeIntegral();
    currTr->AddTimeStep(dst);
    //    printf("Before resetCov: "); currTr->AliExternalTrackParam::Print();
    currTr->ResetCovariance(10000);
    if (RefitTrack(currTr,fITS->GetRMax())) { // propagate to exit from the ITS/TPC screen
      UpdateESDTrack(currTr,AliESDtrack::kITSout);
      fCountITSout++;
    }
    else {
      AliDebug(2,Form("Refit Failed for track %d | ESDtrack#%d (MClb:%d)",itr,fCurrESDtrack->GetID(),fCurrESDtrMClb));
      //currTr->AliExternalTrackParam::Print();
      //currTr->GetWinner()->Print();
    }
    //
  }
  //
  AliInfo(Form("%d ITSout in %d ITSin tracks for %d tried TPC seeds out of %d ESD tracks\n",
	       fCountITSout,fCountITSin,fCountProlongationTrials,nTrESD));
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
  //  AliLog::SetClassDebugLevel("AliITSUTrackerGlo",10);

  AliDebug(1,Form("Will refit inward %d tracks",nTrESD));
  AliITSUTrackHyp *currTr=0;
  //
  for (int itr=0;itr<nTrESD;itr++) {
    fCurrESDtrack = esdEv->GetTrack(itr);
    fCurrESDtrMClb = fCurrESDtrack->GetLabel();
    // Start time integral and add distance from current position to vertex 
    UInt_t trStat = fCurrESDtrack->GetStatus();
    if ( !(trStat & AliESDtrack::kITSout) ) continue;
    if (   trStat & AliESDtrack::kITSrefit ) continue; // already done
    if (  (trStat & AliESDtrack::kTPCout) && !(trStat & AliESDtrack::kTPCrefit) ) continue;
    //
    currTr = GetTrackHyp(itr);
    currTr->AliExternalTrackParam::operator=(*fCurrESDtrack);  // fetch current ESDtrack kinematics
    if (RefitTrack(currTr,fITS->GetRMin())) { // propagate up to inside radius of the beam pipe
      UpdateESDTrack(currTr,AliESDtrack::kITSrefit);
      fCountITSrefit++;
    }
    else {
      AliDebug(2,Form("Refit Failed for track %d |ESDtrack#%d (MClb:%d)",itr,fCurrESDtrack->GetID(),fCurrESDtrMClb));
    }
  }    
  //
  AliInfo(Form("%d ITSrefit in %d ITSout in %d ITSin tracks for %d tried TPC seeds out of %d ESD tracks\n",
	       fCountITSrefit,fCountITSout,fCountITSin,fCountProlongationTrials,nTrESD));
  //
  //  AliLog::SetClassDebugLevel("AliITSUTrackerGlo",0);
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
    int maxNBranches   = fCurrTrackCond->GetMaxBranches(ila);
    int maxNCandidates = fCurrTrackCond->GetMaxCandidates(ila);
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
      AliDebug(2,Form("working on Lr:%d Seed:%d of %d for esdID=%d (MClb:%d) | pT=%.3f",ila,isd,nSeedsUp,esdID,fCurrESDtrMClb,seedUC.Pt()));
      if (!TransportToLayer(&seedUC, fITS->GetLrIDActive(ilaUp), fITS->GetLrIDActive(ila)) ) {
	//
	AliDebug(2,Form("Transport failed | esdID=%d (MClb:%d)",esdID,fCurrESDtrMClb));
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
      AliDebug(2,Form("Will check %d sensors on lr:%d | esdID=%d (MClb:%d)",nsens,ila,esdID,fCurrESDtrMClb));
      //
      seedUC.SetLr(ila);
      //
      double bz = GetBz();
      for (int isn=nsens;isn--;) {
	AliITSURecoSens* sens = hitSens[isn];
	int ncl = sens->GetNClusters();
	if (!ncl) continue;
	seedT = seedUC;
	//
	// We need to propagate the seed to sensor on lrA staying in the frame of the sensor from prev.layer,
	// since the transport matrix should be defined in this frame.
	double xs; // X in the TF of current seed, corresponding to intersection with sensor plane
	if (!seedT.GetTrackingXAtXAlpha(sens->GetXTF(),sens->GetPhiTF(),bz, xs)) {
	  if (AliDebugLevelClass()>2) {
	    printf("Failed on GetTrackingXAtXAlpha: X=%.4f alp=%.3f\n",sens->GetXTF(),sens->GetPhiTF());
	    seedT.Print("etp");
	  }
	  continue;
	}
	if (xs<seedT.GetX()) {
	  if (!PropagateSeed(&seedT,xs,fCurrMass)) continue;
	}
	else { // some low precision tracks may hit the sensor plane outside of the layer radius
	  if (AliDebugLevelClass()>2) {
	    if (!seedT.ContainsFake()) {
	      printf("WRONG PROP on L%d, sens%d of %d: %.4f > %.4f\n",ila,isn,nsens,xs,seedT.GetX());
	      sens->Print();
	      seedT.Print("etp");	    
	    }
	  }
	  if (!seedT.PropagateParamOnlyTo(xs,bz)) continue;
	}
	//	if (!seedT.PropagateToX(xs,bz)) continue;
	//	if (!seedT.Rotate(sens->GetPhiTF())) continue;
	if (!seedT.RotateToAlpha(sens->GetPhiTF())) continue;
	//
	int clID0 = sens->GetFirstClusterId();
	for (int icl=ncl;icl--;) {
	  //	  AliLog::SetClassDebugLevel("AliITSUTrackerGlo",10);
	  int res = CheckCluster(&seedT,ila,clID0+icl);
	  //	  AliLog::SetClassDebugLevel("AliITSUTrackerGlo", 0);
	  //
	  if (res==kStopSearchOnSensor) break;     // stop looking on this sensor
	  if (res==kClusterNotMatching) continue;  // cluster does not match
	  // cluster is matching and it was added to the hypotheses tree
	}
      }
      // cluster search is done. Do we need to have a version of this seed skipping current layer
      if (!NeedToKill(&seedUC,kMissingCluster)) {
	AliITSUSeed* seedSkp = NewSeedFromPool(&seedUC);
	double penalty = -AliITSUReconstructor::GetRecoParam()->GetMissPenalty(ila);
	// to do: make penalty to account for probability to miss the cluster for good reason
	seedSkp->SetChi2Cl(penalty);
	AddProlongationHypothesis(seedSkp,ila);      
      }
      // transfer the new branches of the seed to the hypothesis container
      if (fNBranchesAdded) ValidateAllowedBranches(maxNBranches);
      //
    }
    if (fNCandidatesAdded) ValidateAllowedCandidates(ila,maxNCandidates);
    //    ((TObjArray*)fCurrHyp->GetLayerSeeds(ila))->Sort();
    if (AliDebugLevelClass()>2) { //RS
      printf(">>> All hypotheses on lr %d: \n",ila);
      for (int ih=0;ih<fCurrHyp->GetNSeeds(ila);ih++) {printf(" #%3d ",ih); fCurrHyp->GetSeed(ila,ih)->Print();}
    }
    //
    /*
    for (int ih=0;ih<fCurrHyp->GetNSeeds(ila);ih++) {
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
  fCountProlongationTrials++;
  //
  fCurrMass = esdTr->GetMass();
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
  // transport seed from layerFrom to the entrance (along the direction of the track) of layerTo
  //  
  //
  if (lTo==lFrom) AliFatal(Form("was called with lFrom=%d lTo=%d",lFrom,lTo));
  //
  int dir = lTo > lFrom ? 1:-1;
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom); // this can be 0 when extrapolation from TPC to ITS is requested
  Bool_t checkFirst = kTRUE;
  while(lFrom!=lTo) {
    if (lrFr) {
      if (!GoToExitFromLayer(seed,lrFr,dir,checkFirst)) {
	//printf("FailHere0 Dir=%d\n",dir);
	//seed->Print("etp");
	return kFALSE; // go till the end of current layer
      }
      checkFirst = kFALSE;
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom+=dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    // go the entrance of the layer, assuming no materials in between
    double xToGo = lrTo->GetR(-dir);
    //    double xts = xToGo;
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) {
      //printf("FailHere1: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      //
      
      return kFALSE;
    }
    AliDebug(2,Form("go in dir=%d to R=%.4f(X:%.4f)",dir,lrTo->GetR(-dir), xToGo));
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) {
      //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      return kFALSE;
    }
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
    // go the entrance of the layer, assuming no materials in between
    double xToGo = lrTo->GetR(-dir);
    //    double xts = xToGo;
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) {
      //      printf("FailHere1: %f %f %d\n",xts,xToGo,dir);
      //      seed->Print("etp");
      return kFALSE;
    }
    AliDebug(2,Form("go in dir=%d to R=%.4f(X:%.4f)",dir,lrTo->GetR(-dir), xToGo));
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) {
      //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      return kFALSE;
    }
    lrFr = lrTo;
  }
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::TransportToLayerX(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo, Double_t xStop)
{
  // transport track from layerFrom to the entrance of layerTo but do not pass control parameter X
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
    double xToGo = lrTo->GetR(-dir); // R of the entrance to layer
    //
    //    double xts = xToGo;
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) {
      //      printf("FailHere1: %f %f %d\n",xts,xToGo,dir);
      //      seed->Print("etp");
      return kFALSE;
    }
    if ( (dir>0&&xToGo>xStop) || (dir<0&&xToGo<xStop) ) xToGo = xStop;
    //
    AliDebug(2,Form("go in dir=%d to R=%.4f(X:%.4f)",dir,lrTo->GetR(-dir), xToGo));
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) {
      //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      return kFALSE;
    }
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
    AliDebug(3,Form(" dir:%d Cur: %e Tgt: %e",dir,Sqrt(curR2),xToGo));
    if      (dir>0) { if (curR2-xToGo*xToGo>-fgkToler) return kTRUE; } // on the surface or outside of the layer
    else if (dir<0) { if (xToGo*xToGo-curR2>-fgkToler) return kTRUE; } // on the surface or outside of the layer
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
    if      (dir>0) { if (curR2-xToGo*xToGo>-fgkToler) return kTRUE; } // on the surface or outside of the layer
    else if (dir<0) { if (xToGo*xToGo-curR2>-fgkToler) return kTRUE; } // on the surface or outside of the layer
  }
  AliDebug(2,Form(" dir=%d : from R=%.4f -> R=%.4f",dir,Sqrt(seed->GetX()*seed->GetX() + seed->GetY()*seed->GetY()), xToGo));
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
    if      (dir>0) { if (curR2-xToGo*xToGo>-fgkToler) return kTRUE; } // already passed
    else if (dir<0) { if (xToGo*xToGo-curR2>-fgkToler) return kTRUE; } // already passed
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
    if      (dir>0) { if (curR2-xToGo*xToGo>-fgkToler) return kTRUE; } // already passed
    else if (dir<0) { if (xToGo*xToGo-curR2>-fgkToler) return kTRUE; } // already passed
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
  static AliExternalTrackParam sc;   // seed copy for manipulations
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
  AliDebug(-1,Form("CurrentSize: %d, BookedUpTo: %d, free: %d",fSeedsPool.GetSize(),fSeedsPool.GetEntriesFast(),fNFreeSeeds));
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
  int currLabel = Abs(fCurrESDtrMClb);
  //
  if (cl->GetLabel(0)>=0) {
    for (int i=0;i<3;i++) if (cl->GetLabel(i)>=0 && cl->GetLabel(i)==currLabel) {goodCl = kTRUE; break;}
  }
  Bool_t goodSeed = !track->ContainsFake();
  //
  if (TMath::Abs(cl->GetX())>kTolerX) { // if due to the misalingment X is large, propagate track only
    if (!track->PropagateParamOnlyTo(track->GetX()+cl->GetX(),GetBz())) {
      if (goodCl&&goodSeed && AliDebugLevelClass()>2 ) {
	AliDebug(2,Form("Lost good cl on L:%d failed propagation. |ESDtrack#%d (MClb:%d)",lr,fCurrESDtrack->GetID(),fCurrESDtrMClb)); 
	track->Print("etp");
	cl->Print();
      }
      return kStopSearchOnSensor; // propagation failed, seedT is intact
    }
  }
  double dy = cl->GetY()-track->GetY();
  double dz = cl->GetZ()-track->GetZ();
  //
#ifdef  _FILL_CONTROL_HISTOS_
  int hcOffs = fTrackPhase*kHistosPhase + lr;
  double htrPt=-1;
  if (goodCl&&goodSeed && fCHistoArr /*&& track->GetChi2Penalty()<1e-5*/) {
    htrPt = track->Pt();
    ((TH2*)fCHistoArr->At(kHResY+hcOffs))->Fill(htrPt,dy);
    ((TH2*)fCHistoArr->At(kHResZ+hcOffs))->Fill(htrPt,dz);
    double errY = track->GetSigmaY2();
    double errZ = track->GetSigmaZ2();
    if (errY>0) ((TH2*)fCHistoArr->At(kHResYP+hcOffs))->Fill(htrPt,dy/Sqrt(errY));
    if (errZ>0) ((TH2*)fCHistoArr->At(kHResZP+hcOffs))->Fill(htrPt,dz/Sqrt(errZ));
  }
#endif  
  //
  double dy2 = dy*dy;
  double tol2 = (track->GetSigmaY2() + AliITSUReconstructor::GetRecoParam()->GetSigmaY2(lr))*
    AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadY()*AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadY(); // RS TOOPTIMIZE
  if (dy2>tol2) {                          // the clusters are sorted in Z(col) then in Y(row). 
    if (goodCl&&goodSeed && AliDebugLevelClass()>2) {    
      AliDebug(2,Form("Lost good cl: dy2=%e > tol2=%e |ESDtrack#%d (MClb:%d)",dy2,tol2,fCurrESDtrack->GetID(),fCurrESDtrMClb)); 
      track->Print("etp");
      cl->Print();
      PrintSeedClusters(track);
    }
    if (dy>0) return kStopSearchOnSensor;  // No chance that other cluster of this sensor will match (all Y's will be even larger)
    else      return kClusterNotMatching;   // Other clusters may match
  }
  double dz2 = dz*dz;
  tol2 = (track->GetSigmaZ2() + AliITSUReconstructor::GetRecoParam()->GetSigmaZ2(lr))*
    AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadZ()*AliITSUReconstructor::GetRecoParam()->GetNSigmaRoadZ(); // RS TOOPTIMIZE
  if (dz2>tol2) {
    if (goodCl&&goodSeed && AliDebugLevelClass()>2) {
      AliDebug(2,Form("Lost good cl on L:%d : dz2=%e > tol2=%e |ESDtrack#%d (MClb:%d)",lr,dz2,tol2,fCurrESDtrack->GetID(),fCurrESDtrMClb)); 
      track->Print("etp");
      cl->Print();
      PrintSeedClusters(track);
    }
    return kClusterNotMatching; // Other clusters may match
  }
  //
  // check chi2
  Double_t p[2]={cl->GetY(), cl->GetZ()};
  Double_t cov[3]={cl->GetSigmaY2(), cl->GetSigmaYZ(), cl->GetSigmaZ2()};
  double chi2 = track->GetPredictedChi2(p,cov);
  //
#ifdef  _FILL_CONTROL_HISTOS_
  if (htrPt>0) {
    ((TH2*)fCHistoArr->At(kHChi2Cl+hcOffs))->Fill(htrPt,chi2);
  }
#endif
  //
  if (chi2>AliITSUReconstructor::GetRecoParam()->GetMaxTr2ClChi2(lr)) {
    if (goodCl&&goodSeed && AliDebugLevelClass()>2) {
      AliDebug(2,Form("Lost good cl on L:%d : Chi2=%e > Chi2Max=%e |dy: %+.3e dz: %+.3e |ESDtrack#%d (MClb:%d)\n",
		      lr,chi2,AliITSUReconstructor::GetRecoParam()->GetMaxTr2ClChi2(lr),dy,dz,fCurrESDtrack->GetID(),fCurrESDtrMClb)); 
      track->Print("etp");
      cl->Print("");
      PrintSeedClusters(track);
    }
    return kClusterNotMatching;
  }
  //
  track = NewSeedFromPool(track);  // input track will be reused, use its clone for updates
  if (!track->Update()) {
    if (goodCl&&goodSeed && AliDebugLevelClass()>2) {
      AliDebug(2,Form("Lost good cluster on L:%d : failed to update",lr));
      track->Print("etp");
      cl->Print("");
      PrintSeedClusters(track);
    }
    MarkSeedFree(track);
    return kClusterNotMatching;
  }
  track->SetChi2Cl(chi2);
  track->SetLrClusterID(lr,clID);
  cl->IncreaseClusterUsage();
  //
  track->SetFake(!goodCl);
  //
  AliDebug(2,Form("AddCl(%d) Cl%d lr:%d: dY:%+8.4f dZ:%+8.4f (MC: %5d %5d %5d) |Chi2=%f(%c)| ",
	       goodCl,clID,lr,dy,dz2,cl->GetLabel(0),cl->GetLabel(1),cl->GetLabel(2), chi2, track->IsFake() ? '-':'+'));
  //
  AddSeedBranch(track);
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
    Bool_t seedOK = fCurrTrackCond->CheckPattern(patt);
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
    if (!seed->PropagateToX(x,bz)) return kFALSE;
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
	if (!seed->ApplyMaterialCorrection(xx0,xrho,mass,kFALSE)) {AliInfo("Fail2"); seed->Print("etp"); return kFALSE;}
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
  if (AliDebugLevelClass()>1) {
    AliDebug(2,Form("Before Propagate to X=%f with M=%.3f MaxStep=%.4f MatCorr=%d",xToGo,mass,maxStep,matCorr));
    seed->AliExternalTrackParam::Print();
  }
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
  if (AliDebugLevelClass()>1) {
    AliDebug(2,Form("After Propagate to X=%f",xToGo));
    seed->AliExternalTrackParam::Print();
  }
  return kTRUE;
}

//______________________________________________________________________________
void AliITSUTrackerGlo::SaveCurrentTrackHypotheses()
{
  // RS: shall we clean up killed seeds?
  ((TObjArray*)fCurrHyp->GetLayerSeeds(0))->Sort();
  fCurrHyp = 0;
  // TODO
  
}

//______________________________________________________________________________
void AliITSUTrackerGlo::FinalizeHypotheses()
{
  // select winner for each hypothesis, remove cl. sharing conflicts
  AliDebug(2,"TODO");
  //
  int nh = fHypStore.GetEntriesFast();
  for (int ih=0;ih<nh;ih++) {
    AliITSUTrackHyp* hyp = (AliITSUTrackHyp*) fHypStore.UncheckedAt(ih); 
    if (!hyp || !hyp->DefineWinner()) continue; // TODO
    if (hyp->GetESDTrack()->IsOn(AliESDtrack::kITSin)) continue;
    fCountITSin++;
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
  if (AliDebugLevelClass()>2) {
    printf("Refit %d: Lr:%d (%f) -> Lr:%d (%f)\n",dir,lrStart,rCurr, lrStop,rDest);
    printf("Before refit: "); trc->AliExternalTrackParam::Print();
  }
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
    Bool_t transportedToLayer = kFALSE;
    for (int icl=0;icl<nclLr;icl++) {
      AliITSUClusterPix* clus =  (AliITSUClusterPix*)lr->GetCluster(iclLr[icl]);
      AliITSURecoSens* sens = lr->GetSensorFromID(clus->GetVolumeId());
      if (!tmpTr.Rotate(sens->GetPhiTF())) {
	AliDebug(2,Form("Failed on rotate to %f | ESDtrack#%d (MClb:%d)",sens->GetPhiTF(),fCurrESDtrack->GetID(),fCurrESDtrMClb));
	return kFALSE;
      }
      //
      double xClus = sens->GetXTF()+clus->GetX();
      if (!transportedToLayer) {
	if (ilr!=lrStart && !TransportToLayerX(&tmpTr,lrStart,ilr,xClus)) {
	  AliDebug(2,Form("Failed to transport %d -> %d | ESDtrack#%d (MClb:%d)\n",lrStart,ilr,fCurrESDtrack->GetID(),fCurrESDtrMClb));
	  //tmpTr.AliExternalTrackParam::Print(); trc->GetWinner()->Print("etp");
	  return kFALSE; // go to the entrance to the layer
	}
	lrStart = ilr;
	transportedToLayer = kTRUE;
      }
      //
      if (AliDebugLevelClass()>1) {
	AliDebug(2,Form("Propagate to cl:%d on lr %d Need to go %.4f -> %.4f",icl,ilrA,tmpTr.GetX(),xClus));
	//	printf("Before: "); tmpTr.AliExternalTrackParam::Print();
      }
      //
      if (!PropagateSeed(&tmpTr,xClus,fCurrMass)) {
	AliDebug(2,Form("Failed on propagate to %f (dir=%d) | ESDtrack#%d (MClb:%d)",xClus,dir,fCurrESDtrack->GetID(),fCurrESDtrMClb));
	//tmpTr.AliExternalTrackParam::Print(""); trc->GetWinner()->Print("etp");
	return kFALSE;
      }
      //
#ifdef  _FILL_CONTROL_HISTOS_
      int hcOffs = fTrackPhase*kHistosPhase + ilrA;
      double htrPt=-1;
      if (fCHistoArr && trc->GetLabel()>=0/* && trc->Charge()>0*/) {
	htrPt = tmpTr.Pt();
	double dy = clus->GetY()-tmpTr.GetY();
	double dz = clus->GetZ()-tmpTr.GetZ();
	((TH2*)fCHistoArr->At(kHResY+hcOffs))->Fill(htrPt,dy);
	((TH2*)fCHistoArr->At(kHResZ+hcOffs))->Fill(htrPt,dz);
	double errY = tmpTr.GetSigmaY2();
	double errZ = tmpTr.GetSigmaZ2();
	if (errY>0) ((TH2*)fCHistoArr->At(kHResYP+hcOffs))->Fill(htrPt,dy/Sqrt(errY));
	if (errZ>0) ((TH2*)fCHistoArr->At(kHResZP+hcOffs))->Fill(htrPt,dz/Sqrt(errZ));
      }
#endif  
      //
      double chi2;
      if ( (chi2=tmpTr.Update(clus))<0 ) {
	AliDebug(2,Form("Failed on Update | ESDtrack#%d (MClb:%d)",fCurrESDtrack->GetID(),fCurrESDtrMClb));		
	return kFALSE;
      }
#ifdef  _FILL_CONTROL_HISTOS_
      if (htrPt>0) {
	((TH2*)fCHistoArr->At(kHChi2Cl+hcOffs))->Fill(htrPt,chi2);
      }
#endif 
      if (AliDebugLevelClass()>1) {
	printf("AfterRefit: "); tmpTr.AliExternalTrackParam::Print();
      }
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
      AliDebug(1,Form("Failed on TransportToLayer %d->%d | ESDtrack#%d (MClb:%d)",lrStart,lrStop,fCurrESDtrack->GetID(),fCurrESDtrMClb));
      //tmpTr.AliExternalTrackParam::Print();
      //trc->GetWinner()->Print("etp");
      return kTRUE;
    }    
    if (!GoToExitFromLayer(&tmpTr,fITS->GetLayer(lrStop),dir)) {
      AliDebug(3,Form("Failed on GoToExitFromLayer %d | ESDtrack#%d (MClb:%d)",lrStop,fCurrESDtrack->GetID(),fCurrESDtrMClb));
      return kTRUE; // go till the exit from layer
    }
    //
    //printf("On exit from last layer\n");
    //tmpTr.AliExternalTrackParam::Print();
    // go to the destination radius
    if (!tmpTr.GetXatLabR(rDest,rDest,GetBz(),dir)) return kTRUE;
    if (!PropagateSeed(&tmpTr,rDest,fCurrMass, 100, kFALSE)) return kTRUE;
  }
  trc->AliKalmanTrack::operator=(tmpTr);
  if (AliDebugLevelClass()>2) {
    printf("After refit (now at lr %d): ",lrStop); trc->AliExternalTrackParam::Print();
  }
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
	if (trLb<0) break;
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
  if (nCl && nLab) {
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

//__________________________________________________________________
Bool_t AliITSUTrackerGlo::AddSeedBranch(AliITSUSeed* seed)
{
  // add new prolongation branch for the seed to temporary array. The seeds are sorted according to Compare f-n
  if (fNCandidatesAdded+fNBranchesAdded>=AliITSUTrackCond::kMaxCandidates ) {
    AliError(Form("Number of candidates for current seed reached limit of AliITSUTrackCond::kMaxCandidates=%d, increase it",AliITSUTrackCond::kMaxCandidates));
    return kFALSE;
  }
  AliITSUSeed** branches = &fLayerCandidates[fNCandidatesAdded]; // note: fNCandidatesAdded is incremented after adding all branches of current seed
  int slot=fNBranchesAdded++;
  for (int slotF=slot;slotF--;) { // slotF is always slot-1
    AliITSUSeed* si = branches[slotF];
    if (si->Compare(seed)<0) break; // found the last seed with better quality
    // otherwise, shift the worse seed to the next slot
    branches[slot] = si;
    slot = slotF; // slot should be slotF+1
  }
  // if needed, move worse seeds
  branches[slot] = seed;
  return kTRUE;
  //
}

//__________________________________________________________________
void AliITSUTrackerGlo::ValidateAllowedBranches(Int_t acceptMax)
{
  // keep allowed number of branches for current seed and disable extras
  int nb = Min(fNBranchesAdded,acceptMax);
  //  if (nb<fNBranchesAdded) printf("ValidateAllowedBranches: %d of %d (%d) on lr %d\n",nb,fNBranchesAdded,fNCandidatesAdded,ilr);
  // disable unused branches
  AliITSUSeed** branches = &fLayerCandidates[fNCandidatesAdded];
  for (int ib=nb;ib<fNBranchesAdded;ib++) {
    /*
    if (!branches[ib]->ContainsFake()) {
      printf("Suppress good branch as %d of %d |",ib,fNBranchesAdded); branches[ib]->Print();
      printf("Survivors : \n");
      for (int j=0;j<nb;j++) branches[j]->Print();
    }
    */
    MarkSeedFree(branches[ib]);
  }
  fNCandidatesAdded += nb; // update total candidates counter
  fNBranchesAdded = 0; // reset branches counter
  //
}

//__________________________________________________________________
void AliITSUTrackerGlo::ValidateAllowedCandidates(Int_t ilr, Int_t acceptMax)
{
  // transfer allowed number of branches to hypothesis container
  //
  // sort candidates in increasing order of chi2
  Int_t index[AliITSUTrackCond::kMaxCandidates];
  Float_t chi2[AliITSUTrackCond::kMaxCandidates];
  for (int i=fNCandidatesAdded;i--;) chi2[i] = fLayerCandidates[i]->GetChi2GloNrm();
  Sort(fNCandidatesAdded,chi2,index,kFALSE);
  //
  int nb = Min(fNCandidatesAdded,acceptMax);
  //  if (nb<fNCandidatesAdded) printf("ValidateAllowedCandidates: %d of %d on lr %d\n",nb,fNCandidatesAdded,ilr);
  //
  for (int ib=0;ib<nb;ib++) AddProlongationHypothesis(fLayerCandidates[index[ib]],ilr);
  // disable unused candidates
  for (int ib=nb;ib<fNCandidatesAdded;ib++) {
    /*
    if (!fLayerCandidates[index[ib]]->ContainsFake()) {
      printf("Suppress good candidate as %d of %d |",index[ib],fNCandidatesAdded); fLayerCandidates[index[ib]]->Print();
    }
    */
    MarkSeedFree(fLayerCandidates[index[ib]]);    
  }
  fNCandidatesAdded = 0; // reset candidates counter
  //
}

//__________________________________________________________________
void AliITSUTrackerGlo::FlagSplitClusters()
{
  // set special bit on split clusters using MC info
  for (int ilr=fITS->GetNLayersActive();ilr--;) {
    int nsplit=0;
    AliITSURecoLayer* lr = fITS->GetLayerActive(ilr);
    for (int isn=lr->GetNSensors();isn--;) {
      AliITSURecoSens* sens = lr->GetSensor(isn);
      int nCl = sens->GetNClusters();
      if (!nCl) continue;
      int cl0 = sens->GetFirstClusterId();
      for (int icl=nCl;icl--;) {
	AliITSUClusterPix *cl = (AliITSUClusterPix*)lr->GetCluster(cl0+icl);
	for (int icl1=icl;icl1--;) {
	  AliITSUClusterPix *cl1 = (AliITSUClusterPix*)lr->GetCluster(cl0+icl1);
	  if (cl->HasCommonTrack(cl1)) {
	    if (!cl->IsSplit())  nsplit++;
	    if (!cl1->IsSplit()) nsplit++;
	    cl->SetSplit();
	    cl1->SetSplit();
	  }
	}
      }
    }
    AliInfo(Form("%4d out of %4d clusters are split",nsplit,lr->GetNClusters()));
  }
  //
}

//__________________________________________________________________
Bool_t AliITSUTrackerGlo::ContainsSplitCluster(const AliITSUSeed* seed, Int_t maxSize)
{
  // check if the seed contains split cluster with size < maxSize
  int lrID,clID;
  if ( (clID=seed->GetLrCluster(lrID))>=0 ) {
    AliITSUClusterPix* cl = (AliITSUClusterPix*)fITS->GetLayerActive(lrID)->GetCluster(clID);
    if (cl->IsSplit() && cl->GetNPix()<maxSize ) return kTRUE;
  }
  return seed->GetParent() ? ContainsSplitCluster((AliITSUSeed*)seed->GetParent(),maxSize) : kFALSE;
  //
}

//__________________________________________________________________
void AliITSUTrackerGlo::PrintSeedClusters(const AliITSUSeed* seed, Option_t* option)
{
  // print seeds clusters
  int lrID,clID;
  if ( (clID=seed->GetLrCluster(lrID))>=0 ) {
    AliITSUClusterPix* cl = (AliITSUClusterPix*)fITS->GetLayerActive(lrID)->GetCluster(clID);
    cl->Print(option);
  }
  if (seed->GetParent()) PrintSeedClusters((AliITSUSeed*)seed->GetParent(), option);
  //
}

#ifdef  _FILL_CONTROL_HISTOS_
//__________________________________________________________________
void AliITSUTrackerGlo::BookControlHistos()
{
  // book special control histos
  if (!fCHistoArr) { // create control histos
    const int kNResDef=7;
    const double kResDef[kNResDef]={0.05,0.05,0.3, 0.05,1,0.5,1.5};
    fCHistoArr = new TObjArray();
    fCHistoArr->SetOwner(kTRUE);
    const double ptMax=10;
    const double plMax=10;
    const double chiMax=100;
    const int nptbins=50;
    const int nresbins=400;
    const int nplbins=50;
    const int nchbins=200;
    int nblr = fITS->GetNLayersActive();
    TString ttl;
    for (int stp=0;stp<kNTrackingPhases;stp++) {
      for (int ilr=0;ilr<nblr;ilr++) {
	int hoffs = stp*kHistosPhase + ilr;
	double mxdf = ilr>=kNResDef ? kResDef[kNResDef-1] : kResDef[ilr];
	ttl = Form("S%d_residY%d",stp,ilr);
	TH2F* hdy = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax,nresbins,-mxdf,mxdf);
	fCHistoArr->AddAtAndExpand(hdy,hoffs + kHResY);
	hdy->SetDirectory(0);
	//
	ttl = Form("S%d_residYPull%d",stp,ilr);	
	TH2F* hdyp = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax,nplbins,-plMax,plMax);
	fCHistoArr->AddAtAndExpand(hdyp,hoffs + kHResYP);
	hdyp->SetDirectory(0);
	//
	ttl = Form("S%d_residZ%d",stp,ilr);	
	TH2F* hdz = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax,nresbins,-mxdf,mxdf);
	fCHistoArr->AddAtAndExpand(hdz,hoffs + kHResZ);
	hdz->SetDirectory(0);
	//
	ttl = Form("S%d_residZPull%d",stp,ilr);		
	TH2F* hdzp = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax,nplbins,-plMax,plMax);
	hdzp->SetDirectory(0);
	fCHistoArr->AddAtAndExpand(hdzp,hoffs + kHResZP);
	//
	ttl = Form("S%d_chi2Cl%d",stp,ilr);		
	TH2F* hchi = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax, nchbins,0.,chiMax);
	hchi->SetDirectory(0);
	fCHistoArr->AddAtAndExpand(hchi,hoffs + kHChi2Cl);
      }
    }
  }
  //  
}
#endif
