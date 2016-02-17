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
#include "AliITSUSeed.h"
#include "AliITSUClusterPix.h"
#include "AliITSUGeomTGeo.h"
#include "AliCodeTimer.h"
#include "AliRefArray.h"
using namespace AliITSMFTAux;
using namespace TMath;

//----------------- tmp stuff -----------------

ClassImp(AliITSUTrackerGlo)

const Double_t AliITSUTrackerGlo::fgkToler =  1e-6;// tolerance for layer on-surface check


//_________________________________________________________________________
AliITSUTrackerGlo::AliITSUTrackerGlo(AliITSUReconstructor* rec)
:  fReconstructor(rec)
  ,fITS(0)
  ,fMatLUT(0)
  ,fCurrESDtrack(0)
  ,fCurrESDtrMClb(kDummyLabel)
  ,fCurrMass(kPionMass)
  ,fCountProlongationTrials(0)
  ,fCountITSin(0)
  ,fCountITSout(0)
  ,fCountITSrefit(0)
  ,fNTracksESD(0)
  ,fHypStore(100)
  ,fLayerMaxCandidates(1000)
  ,fLayerCandidates(0)
  ,fNBranchesAdded(0)
  ,fNCandidatesAdded(0)
  ,fCurrHyp(0)
  ,fWorkHyp(0)
  ,fSeedsPool("AliITSUSeed",0)
  ,fFreeSeedsID(0)
  ,fESDIndex(0)
  ,fNFreeSeeds(0)
  ,fLastSeedID(0)
  ,fNLrActive(0)
  ,fDefTrackConds(0)
  ,fCurrTrackCond(0)
  ,fCurrActLrID(-1)
  ,fCurrLayer(0)
  ,fTrackPhaseID(-1)
  ,fCurrPassID(-1)
  ,fUseMatLUT(kFALSE)
#ifdef  _ITSU_TUNING_MODE_
  ,fCHistoArrCorr(0)
  ,fCHistoArrFake(0)
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
  delete[] fLayerCandidates;
  if (fWorkHyp) fWorkHyp->SetTPCSeed(0); // this hypothesis does not own the seed
  delete fWorkHyp;
  delete fMatLUT;
  //
#ifdef  _ITSU_TUNING_MODE_
  if (fCHistoArrCorr || fCHistoArrFake) {
    TFile* ctrOut = TFile::Open("itsuControlHistos.root","recreate");
    ctrOut->cd();
    AliInfo("Storing control histos");
    //    ctrOut->WriteObject(fCHistoArr,"controlH","kSingleKey");
    if (fCHistoArrCorr) {
      fCHistoArrCorr->Write();
      delete fCHistoArrCorr;
    }
    if (fCHistoArrFake) {
      fCHistoArrFake->Write();
      delete fCHistoArrFake;
    }
    ctrOut->Close();
    delete ctrOut;
    fCHistoArrCorr = 0;
    fCHistoArrFake = 0;    
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
  fNLrActive = fITS->GetNLayersActive();
  fWorkHyp = new AliITSUTrackHyp(fNLrActive);
  // create material lookup table
  const int kNTest = 1000;
  const double kStepsPerCM=5;
  fMatLUT  = new AliITSUMatLUT(fITS->GetRMin(),fITS->GetRMax(),Nint(kStepsPerCM*(fITS->GetRMax()-fITS->GetRMin())));
  double zmn = 1e6;
  for (int ilr=fITS->GetNLayers();ilr--;) {
    AliITSURecoLayer* lr = fITS->GetLayer(ilr);
    if (zmn>Abs(lr->GetZMin())) zmn = Abs(lr->GetZMin());
    if (zmn>Abs(lr->GetZMax())) zmn = Abs(lr->GetZMax());
  }
  fMatLUT->FillData(kNTest,-zmn,zmn);
  //
  if (fLayerMaxCandidates<1) fLayerMaxCandidates = 1000;
  fLayerCandidates = new AliITSUSeed*[fLayerMaxCandidates];
  fSeedsPool.ExpandCreateFast(1000); // RS TOCHECK
  fFreeSeedsID.Set(1000);
  fESDIndex.Set(1000);
  //
}

//_________________________________________________________________________
void AliITSUTrackerGlo::CreateDefaultTrackCond()
{
  // creates default tracking conditions to be used when recoparam does not provide them

  AliITSUTrackCond* cond = new AliITSUTrackCond();
  //
  cond->SetNLayers(fNLrActive);
  cond->AddNewCondition(fNLrActive);
  cond->AddGroupPattern( 0xffff, fNLrActive ); // require all layers hit
  cond->Init();
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
  AliCodeTimerAuto("",0);
  SetTrackingPhase(kClus2Tracks);
  //
#ifdef  _ITSU_TUNING_MODE_
  if (!fCHistoArrCorr) fCHistoArrCorr = BookControlHistos("Corr");
  if (!fCHistoArrFake) fCHistoArrFake = BookControlHistos("Fake");
#endif
  static int evID = 0;
  static TArrayF esdTrPt(fESDIndex.GetSize()); 
  //
  TObjArray *trackConds = 0;
  //
  fCountProlongationTrials = 0;
  fCountITSin = 0;
  fCountITSout = 0;
  fCountITSrefit = 0;
  //
  ResetSeedsPool();
  fNTracksESD = esdEv->GetNumberOfTracks();
  AliInfo(Form("Will try to find prolongations for %d tracks",fNTracksESD));
  int nTrackCond = AliITSUReconstructor::GetRecoParam()->GetNTrackingConditions();
  fUseMatLUT = AliITSUReconstructor::GetRecoParam()->GetUseMatLUT(fTrackPhaseID);
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
  if (fHypStore.GetSize()<fNTracksESD) fHypStore.Expand(fNTracksESD+100);
  //
  fITS->SortClusters(AliITSUClusterPix::SortModeIdTrkYZ());
  fITS->ProcessClusters();
  //
#ifdef  _ITSU_TUNING_MODE_
  FlagSplitClusters(); // tmp RS
#endif
  //
  // the tracks will be reconstructed in decreasing pt order, sort them
  if (fESDIndex.GetSize()<fNTracksESD) {
    fESDIndex.Set(fNTracksESD+200);
    esdTrPt.Set(fNTracksESD+200);
  }
  int* esdTrackIndex = fESDIndex.GetArray();
  float* trPt = esdTrPt.GetArray();
  for (int itr=fNTracksESD;itr--;) {
    AliESDtrack* esdTr = esdEv->GetTrack(itr);
    esdTr->SetStatus(AliESDtrack::kITSupg);   // Flag all tracks as participating in the ITSup reco
    trPt[itr] = esdTr->Pt();
  }
  Sort(fNTracksESD,trPt,esdTrackIndex,kTRUE);    
  //
  for (int icnd=0;icnd<nTrackCond;icnd++) {
    fCurrPassID = icnd;
    fCurrTrackCond = (AliITSUTrackCond*)trackConds->UncheckedAt(icnd);
    if (!fCurrTrackCond->IsInitDone()) fCurrTrackCond->Init();
    // select ESD tracks to propagate
    for (int itr=0;itr<fNTracksESD;itr++) {
      int trID = esdTrackIndex[itr];
      fCurrESDtrack = esdEv->GetTrack(trID);
      fCurrESDtrMClb = fCurrESDtrack->GetLabel();
      //
      
      if (!NeedToProlong(fCurrESDtrack,trID)) continue;  // are we interested in this track?
      /*
	// if specific tracks should be checked, create a global array int wtc[] = {evITS*10000+trID, ...};
      Bool_t dbg = kFALSE;
      int nwtc = sizeof(wtc)/sizeof(int);

      for (int k=0;k<nwtc;k++) {
	if (wtc[k]==evID*10000+trID) {
	  dbg = kTRUE;
	  printf("At esdTr: %d MC: %d\n",wtc[k],fCurrESDtrMClb);
	  break;
	}
      }
      AliLog::SetClassDebugLevel("AliITSUTrackerGlo",dbg ? 10:0);
      */
#ifdef _ITSU_DEBUG_
      AliDebug(1,Form("Processing track %d(esd%d) | M=%.3f Pt=%.3f | MCLabel: %d",itr,trID,fCurrESDtrack->GetMassForTracking(kTRUE),fCurrESDtrack->Pt(),fCurrESDtrMClb));//RS
#endif
      FindTrack(fCurrESDtrack, trID);
    }   
    //
#ifdef _ITSU_DEBUG_
    if (AliDebugLevelClass()>+2) {
      AliInfo(Form("SeedsPool: %d, BookedUpTo: %d, free: %d",fSeedsPool.GetSize(),fSeedsPool.GetEntriesFast(),fNFreeSeeds));
      fHypStore.Print();
    }
#endif
    FinalizeHypotheses();
#ifdef  _ITSU_TUNING_MODE_
    CheckClusterUsage(); //!!RS
#endif
    //
    //AliLog::SetClassDebugLevel("AliITSUTrackerGlo",0);     // in case wtc array was defined, uncomment this
  }
  //
  AliInfo(Form("%d ITSin for %d tried TPC seeds out of %d ESD tracks\n",fCountITSin,fCountProlongationTrials,fNTracksESD));
  evID++;
  return 0;
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::PropagateBack(AliESDEvent *esdEv)
{
  //
  // Do outward fits in ITS
  AliCodeTimerAuto("",0);
  //
  SetTrackingPhase(kPropBack);
  fNTracksESD = esdEv->GetNumberOfTracks();
  fUseMatLUT = AliITSUReconstructor::GetRecoParam()->GetUseMatLUT(fTrackPhaseID);
#ifdef _ITSU_DEBUG_
  AliDebug(1,Form("Will propagate back %d tracks",fNTracksESD));
#endif
  //
  double bz0 = GetBz();
  Double_t xyzTrk[3],xyzVtx[3]={GetX(),GetY(),GetZ()};
  AliITSUTrackHyp dummyTr;
  const double kWatchStep=10.; // for larger steps watch arc vs segment difference
  Double_t times[AliPID::kSPECIESC];
  //
  for (int itr=0;itr<fNTracksESD;itr++) {
    fCurrESDtrack = esdEv->GetTrack(itr);
    fCurrESDtrMClb = fCurrESDtrack->GetLabel();
    // Start time integral and add distance from current position to vertex 
    if (fCurrESDtrack->IsOn(AliESDtrack::kITSout)) continue; // already done
    if (!fCurrESDtrack->IsOn(AliESDtrack::kTPCin)) continue; // skip ITS s.a.
    //
    fCurrESDtrack->GetXYZ(xyzTrk); 
    Double_t dst = 0.;     // set initial track lenght, tof
    {
      double dxs = xyzTrk[0] - xyzVtx[0];
      double dys = xyzTrk[1] - xyzVtx[1];
      //      double dzs = xyzTrk[2] - xyzVtx[2];
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
      // RS: temporary hack since we don't have SPD vertex:      
      //      dst += dzs*dzs;
      dst *= 1+fCurrESDtrack->GetTgl()*fCurrESDtrack->GetTgl();
      dst = Sqrt(dst); 
    }
    //    
    fCurrESDtrack->SetStatus(AliESDtrack::kTIME);
    //
    if (!fCurrESDtrack->IsOn(AliESDtrack::kITSin)) { // case of tracks w/o ITS prolongation: just set the time integration
      dummyTr.AliExternalTrackParam::operator=(*fCurrESDtrack);
      dummyTr.StartTimeIntegral();
      dummyTr.AddTimeStep(dst);
      dummyTr.GetIntegratedTimes(times,AliPID::kSPECIESC); 
      fCurrESDtrack->SetIntegratedTimes(times);
      fCurrESDtrack->SetIntegratedLength(dummyTr.GetIntegratedLength());
      continue;
    }
    //
    fCurrHyp = GetTrackHyp(itr);
    fCurrMass = fCurrHyp->GetMass();
    fCurrHyp->StartTimeIntegral();
    fCurrHyp->AddTimeStep(dst);
    fCurrHyp->ResetCovariance(10000);
    int nclFit = 0;
    double chi2 = RefitTrack(fCurrHyp,fITS->GetRMax(),nclFit);
    if (chi2>0) { // propagate to exit from the ITS/TPC screen
      int ndf = nclFit*2-5;
      if (ndf>0) chi2 /= ndf;
      fCurrHyp->SetChi2(chi2);
      UpdateESDTrack(fCurrHyp,AliESDtrack::kITSout);
    }
    else {
#ifdef _ITSU_DEBUG_
      AliDebug(2,Form("Refit Failed for track %d | ESDtrack#%d (MClb:%d)",itr,fCurrESDtrack->GetID(),fCurrESDtrMClb));
#endif
      //fCurrHyp->AliExternalTrackParam::Print();
      //fCurrHyp->GetWinner()->Print();
    }
    //
  }
  //
  AliInfo(Form("%d ITSout in %d ITSin tracks for %d tried TPC seeds out of %d ESD tracks\n",
	       fCountITSout,fCountITSin,fCountProlongationTrials,fNTracksESD));
  //
  return 0;
}

//_________________________________________________________________________
Int_t AliITSUTrackerGlo::RefitInward(AliESDEvent *esdEv)
{
  //
  // refit the tracks inward, using current cov.matrix
  AliCodeTimerAuto("",0);
  //
  SetTrackingPhase(kRefitInw);
  fUseMatLUT = AliITSUReconstructor::GetRecoParam()->GetUseMatLUT(fTrackPhaseID);
  fNTracksESD = esdEv->GetNumberOfTracks();
  //  AliLog::SetClassDebugLevel("AliITSUTrackerGlo",10);

  AliDebug(1,Form("Will refit inward %d tracks",fNTracksESD));
  //  Bool_t uselogMS = AliExternalTrackParam::GetUseLogTermMS();
  //  AliExternalTrackParam::SetUseLogTermMS(kTRUE);
  //
  for (int itr=0;itr<fNTracksESD;itr++) {
    fCurrESDtrack = esdEv->GetTrack(itr);
    fCurrESDtrMClb = fCurrESDtrack->GetLabel();
    // Start time integral and add distance from current position to vertex 
    UInt_t trStat = fCurrESDtrack->GetStatus();
    if ( !(trStat & AliESDtrack::kITSout) ) continue;
    if (   trStat & AliESDtrack::kITSrefit ) continue; // already done
    if ( !(trStat & AliESDtrack::kTPCin)   ) continue; // skip ITS s.a.
    if (  (trStat & AliESDtrack::kTPCout) && !(trStat & AliESDtrack::kTPCrefit) ) continue;
    //
    fCurrHyp  = GetTrackHyp(itr);
    fCurrMass = fCurrHyp->GetMass();
    *fCurrHyp = *fCurrESDtrack;  // fetch current ESDtrack kinematics
    //
    int nclFit = 0;
    //    fCurrHyp->ResetCovariance(1000); // in case we want to do SA refit
    double chi2 = RefitTrack(fCurrHyp,fITS->GetRMin(),nclFit);
    if (chi2>0) { // propagate up to inside radius of the beam pipe      
      int ndf = nclFit*2-5;
      if (ndf>0) chi2 /= ndf;
      fCurrHyp->SetChi2(chi2);
      UpdateESDTrack(fCurrHyp,AliESDtrack::kITSrefit);
    }
    else {
#ifdef _ITSU_DEBUG_
      AliDebug(2,Form("Refit Failed for track %d |ESDtrack#%d (MClb:%d)",itr,fCurrESDtrack->GetID(),fCurrESDtrMClb));
#endif
    }
  }    
  //  AliExternalTrackParam::SetUseLogTermMS(uselogMS);
  //
  AliInfo(Form("%d ITSrefit in %d ITSout in %d ITSin tracks for %d tried TPC seeds out of %d ESD tracks\n",
	       fCountITSrefit,fCountITSout,fCountITSin,fCountProlongationTrials,fNTracksESD));
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
  // Remove clusters from the memory 
  //
  /*
  AliITSURecoDet *det=fReconstructor->GetITSInterface();
  Int_t nlayers=det->GetNLayersActive();
  for (Int_t i=0; i<nlayers; i++) {
      TClonesArray *clusters=*(det->GetLayerActive(i)->GetClustersAddress());
      clusters->Delete();
  }
  */
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
Bool_t AliITSUTrackerGlo::NeedToProlong(AliESDtrack* esdTr, Int_t esdID)
{
  // do we need to match this track to ITS?
  //
  static double bz = GetBz();
  //
  // is it already prolonged
  if (esdTr->IsOn(AliESDtrack::kITSin)) return kFALSE;
  AliITSUTrackHyp* hyp = GetTrackHyp(esdID); // in case there is unfinished hypothesis from prev.pass, clean it
  if (hyp) {
    CleanHypothesis(hyp);
    if (hyp->GetSkip()) return kFALSE; // need to skip this
    return kTRUE;
  }
  //
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
  //  if (esdTr->Pt()<3) return kFALSE;//RS
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
  //  AliLog::SetClassDebugLevel("AliITSUTrackerGlo",10);
  AliITSUTrackHyp* hypTr = InitHypothesis(esdTr,esdID);  // initialize prolongations hypotheses tree
  if (!hypTr || hypTr->GetSkip()) return;
  //
  double bz = GetBz();
  fCurrHyp = fWorkHyp;
  fCurrHyp->InitFrom(hypTr);
  //
  AliITSURecoSens *hitSens[AliITSURecoLayer::kMaxSensMatching];
  //
  int ilaUp = fNLrActive;     // previous active layer really checked (some may be excluded!)
  //
  for (int ila=fNLrActive;ila--;ilaUp=fCurrActLrID) {
    if (fCurrTrackCond->IsLayerExcluded(ila)) continue;
    fCurrActLrID = ila;
    fCurrLayer = fITS->GetLayerActive(ila);
    Bool_t noClSharing = fCurrTrackCond->GetClSharing(ila)==0;
    //
    // for the outermost layer the seed is created from the ESD track
    int nSeedsUp = (ilaUp==fNLrActive) ? 1 : fCurrHyp->GetNSeeds(ilaUp);
    int maxNBranches   = fCurrTrackCond->GetMaxBranches(ila);
    int maxNCandidates = fCurrTrackCond->GetMaxCandidates(ila);
    //
    for (int isd=0;isd<nSeedsUp;isd++) {
      AliITSUSeed* seedU;
      if (ilaUp==fNLrActive) {
	seedU = 0;
	seedUC.InitFromSeed(fCurrHyp->GetTPCSeed()); // init with TPC seed from ref.R
      }
      else {
	seedU = fCurrHyp->GetSeed(ilaUp,isd);  // seed on prev.active layer to prolong	
	if (seedU->IsKilled()) continue;
	seedUC = *seedU;                       // its copy will be prolonged
	seedUC.SetParent(seedU);	
      }
      seedUC.ResetFMatrix();                    // reset the matrix for propagation to next layer
      // go till next active layer
#ifdef _ITSU_DEBUG_
      AliDebug(2,Form("working on Lr:%d Seed:%d of %d for esdID=%d (MClb:%d) | pT=%.3f",ila,isd,nSeedsUp,esdID,fCurrESDtrMClb,seedUC.Pt()));
#endif
      if (!TransportToLayer(&seedUC, fITS->GetLrIDActive(ilaUp), fITS->GetLrIDActive(ila)) ) { // external seed already prolonged
	//
#ifdef _ITSU_DEBUG_
	AliDebug(2,Form("Transport failed | esdID=%d (MClb:%d)",esdID,fCurrESDtrMClb));
#endif
	// Check if the seed satisfies to track definition
	if (NeedToKill(&seedUC,kTransportFailed) && seedU) KillSeed(seedU,kTRUE);
	continue; // RS TODO: decide what to do with tracks stopped on higher layers w/o killing
      }
      if (!GetRoadWidth(&seedUC, ila)) { // failed to find road width on the layer
	if (NeedToKill(&seedUC,kRWCheckFailed) && seedU) KillSeed(seedU,kTRUE);
	continue;
      }
      /*
      //RS toremove
      int mcquest = -1;
      if (!seedUC.ContainsFake() && AliDebugLevelClass()>2) {
	mcquest = fCurrESDtrMClb;
	seedUC.Print("etp");
	printf("FSParams: "); for (int ip=0;ip<kNTrImpData;ip++) printf("%+e ",fTrImpData[ip]); printf("\n");
      }
      //
      int nsens = fCurrLayer->FindSensors(&fTrImpData[kTrPhi0], hitSens, mcquest);  // find detectors which may be hit by the track
      */
      int nsens = fCurrLayer->FindSensors(&fTrImpData[kTrPhi0], hitSens);  // find detectors which may be hit by the track
#ifdef _ITSU_DEBUG_
      AliDebug(2,Form("Will check %d sensors on lr:%d | esdID=%d (MClb:%d)",nsens,ila,esdID,fCurrESDtrMClb));
#endif
      //
      seedUC.SetLr(ila);
      //
      for (int isn=nsens;isn--;) {
	AliITSURecoSens* sens = hitSens[isn];
	int ncl = sens->GetNClusters();
	if (!ncl) continue;
	seedT = seedUC;
	//
	// We need to propagate the seed to sensor on fCurrLayer staying in the frame of the sensor from prev.layer,
	// since the transport matrix should be defined in this frame.
	double xs; // X in the TF of current seed, corresponding to intersection with sensor plane
	if (!seedT.GetTrackingXAtXAlpha(sens->GetXTF(),sens->GetPhiTF(),bz, xs)) {
#ifdef _ITSU_DEBUG_
	  if (AliDebugLevelClass()>2) {
	    printf("Failed on GetTrackingXAtXAlpha: X=%.4f alp=%.3f\n",sens->GetXTF(),sens->GetPhiTF());
	    seedT.Print("etp");
	  }
#endif
	  continue;
	}
	if (xs<seedT.GetX()) {
	  if (!PropagateSeed(&seedT,xs,fCurrMass)) continue;
	}
	else { // some low precision tracks may hit the sensor plane outside of the layer radius
#ifdef _ITSU_DEBUG_
	  if (AliDebugLevelClass()>2) {
	    if (!seedT.ContainsFake()) {
	      printf("WRONG PROP on L%d, sens%d of %d: %.4f > %.4f\n",ila,isn,nsens,xs,seedT.GetX());
	      sens->Print();
	      seedT.Print("etp");	    
	    }
	  }
#endif
	  if (!seedT.PropagateParamOnlyTo(xs,bz)) continue;
	}
	//	if (!seedT.PropagateToX(xs,bz)) continue;
	//	if (!seedT.Rotate(sens->GetPhiTF())) continue;
	if (!seedT.RotateToAlpha(sens->GetPhiTF())) continue;
	//
	int clID0 = sens->GetFirstClusterId();
	for (int icl=ncl;icl--;) { // don't change the order, clusters are sorted
	  //AliLog::SetClassDebugLevel("AliITSUTrackerGlo",10);
	  int clID = clID0 + icl;
	  if (noClSharing && fCurrLayer->GetCluster(clID)->IsClusterUsed()) continue;
	  int res = CheckCluster(&seedT,ila,clID0+icl);
	  //AliLog::SetClassDebugLevel("AliITSUTrackerGlo", 0);
	  //
	  if (res==kStopSearchOnSensor) break;     // stop looking on this sensor
	  if (res==kClusterNotMatching) continue;  // cluster does not match
	  // cluster is matching and it was added to the hypotheses tree
	}
      }
      // cluster search is done. Do we need to have a version of this seed skipping current layer
      if (!NeedToKill(&seedUC,kMissingCluster)) {
	AliITSUSeed* seedSkp = NewSeedFromPool(&seedUC);
	if (nsens) { // to do: make penalty to account for probability to miss the cluster for good reason
	  double penalty = -fCurrTrackCond->GetMissPenalty(ila);
	  seedSkp->SetChi2Cl(penalty);
	}
	if (seedSkp->GetChi2GloNrm()>fCurrTrackCond->GetMaxChi2GloNrm(ila)) {
	  MarkSeedFree(seedSkp);
	}
	else AddSeedBranch(seedSkp);
      }
      // transfer the new branches of the seed to the hypothesis container
      if (fNBranchesAdded) ValidateAllowedBranches(maxNBranches);
      //
    }
    if (fNCandidatesAdded) ValidateAllowedCandidates(ila,maxNCandidates);
    //    ((TObjArray*)fCurrHyp->GetLayerSeeds(ila))->Sort();
#ifdef _ITSU_DEBUG_
    if (AliDebugLevelClass()>2) { //RS
      printf(">>> All hypotheses on lr %d: \n",ila);
      for (int ih=0;ih<fCurrHyp->GetNSeeds(ila);ih++) {printf(" #%3d ",ih); fCurrHyp->GetSeed(ila,ih)->Print();}
    }
#endif
    //
    /*
    for (int ih=0;ih<fCurrHyp->GetNSeeds(ila);ih++) {
      if (ila!=0) continue;
      double vecL[5] = {0};
      double matL[15] = {0};
      AliITSUSeed* sp = fCurrHyp->GetSeed(ila,ih);
      while(sp->GetParent()) {
	sp->Smooth(vecL,matL);
	if (sp->GetLayerID()>=fNLrActive-1) break;
	sp = (AliITSUSeed*)sp->GetParent();
      }
      */
  }
  //
  SaveReducedHypothesesTree(hypTr);
#ifdef _ITSU_DEBUG_
  if (AliDebugLevelClass()>1) {
    printf("\nSaved hypotheses for esdTrack %d (MCLab:%d)\n",esdID,fCurrESDtrMClb);
    hypTr->Print("l");
  }
#endif
  fCurrHyp = 0;
  //  AliLog::SetClassDebugLevel("AliITSUTrackerGlo",0);
  //
}

//_________________________________________________________________________
AliITSUTrackHyp* AliITSUTrackerGlo::InitHypothesis(AliESDtrack *esdTr, Int_t esdID)
{
  // init prolongaion candidates finding for single seed
  AliITSUTrackHyp* hyp = GetTrackHyp(esdID);
  if (hyp) return hyp;
  //
  fCountProlongationTrials++;
  //
  fCurrMass = esdTr->GetMassForTracking();
  //
  hyp = new AliITSUTrackHyp(fNLrActive);
  hyp->SetESDTrack(esdTr);
  hyp->SetUniqueID(esdID);
  hyp->SetMass(fCurrMass);
  hyp->SetTPCSeed( new AliExternalTrackParam(*esdTr) );
  SetTrackHyp(hyp,esdID);
  Bool_t lut = fUseMatLUT; // for propagation from TPC use TGeo material budget always
  fUseMatLUT = kFALSE;
  if (!TransportToLayer(hyp->GetTPCSeed(),fITS->GetNLayers(), fITS->GetLrIDActive(fNLrActive-1), fITS->GetRITSTPCRef())) hyp->SetSkip(); // propagate to outer R of ITS
  fUseMatLUT = lut;
  //
  return hyp;
  // TO DO
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::TransportToLayer(AliITSUSeed* seed, Int_t lFrom, Int_t lTo, Double_t rLim)
{
  // transport seed from layerFrom to the entrance (along the direction of the track) of layerTo or to rLim (if>0), wathever is closer
  //  
  //
  if (lTo==lFrom) AliFatal(Form("was called with lFrom=%d lTo=%d",lFrom,lTo));
  //
  int dir = lTo > lFrom ? 1:-1;
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom); // this can be 0 when extrapolation from TPC to ITS is requested
  Bool_t checkFirst = kTRUE;
  Bool_t limReached = kFALSE;
  while(lFrom!=lTo) {
    if (lrFr) {
      if (!GoToExitFromLayer(seed,lrFr,dir,checkFirst))	return kFALSE; // go till the end of current layer
      checkFirst = kFALSE;
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom+=dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    // go the entrance of the layer, assuming no materials in between
    double xToGo = lrTo->GetR(-dir);
    if (rLim>0) {
      if (dir>0) {
	if (rLim<xToGo) {xToGo = rLim; limReached = kTRUE;}
      }
      else {
	if (rLim>xToGo) {xToGo = rLim; limReached = kTRUE;}
      }
    }
    //    double xts = xToGo;
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) {
      //printf("FailHere1: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      //      
      return kFALSE;
    }
#ifdef _ITSU_DEBUG_
    AliDebug(2,Form("go in dir=%d to R=%.4f(X:%.4f)",dir,lrTo->GetR(-dir), xToGo));
#endif
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) {
      //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      return kFALSE;
    }
    lrFr = lrTo;
    if (limReached) break;
  }
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSUTrackerGlo::TransportToLayer(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo, Double_t rLim)
{
  // transport track from layerFrom to the entrance of layerTo or to rLim (if>0), wathever is closer
  //  
  if (lTo==lFrom) AliFatal(Form("was called with lFrom=%d lTo=%d",lFrom,lTo));
  //
  int dir = lTo > lFrom ? 1:-1;
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom); // this can be 0 when extrapolation from TPC to ITS is requested
  Bool_t checkFirst = kTRUE;
  Bool_t limReached = kFALSE;
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
    if (rLim>0) {
      if (dir>0) {
	if (rLim<xToGo) {xToGo = rLim; limReached = kTRUE;}
      }
      else {
	if (rLim>xToGo) {xToGo = rLim; limReached = kTRUE;}
      }
    }
    //    double xts = xToGo;
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) {
      //      printf("FailHere1: %f %f %d\n",xts,xToGo,dir);
      //      seed->Print("etp");
      return kFALSE;
    }
#ifdef _ITSU_DEBUG_
    AliDebug(2,Form("go in dir=%d to R=%.4f(X:%.4f)",dir,lrTo->GetR(-dir), xToGo));
#endif
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) {
      //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      return kFALSE;
    }
    lrFr = lrTo;
    if (limReached) break;
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
#ifdef _ITSU_DEBUG_
    AliDebug(2,Form("go in dir=%d to R=%.4f(X:%.4f)",dir,lrTo->GetR(-dir), xToGo));
#endif
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
#ifdef _ITSU_DEBUG_
    AliDebug(3,Form(" dir:%d Cur: %e Tgt: %e",dir,Sqrt(curR2),xToGo));
#endif
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
#ifdef _ITSU_DEBUG_
  AliDebug(2,Form(" dir=%d : from R=%.4f -> R=%.4f",dir,Sqrt(seed->GetX()*seed->GetX() + seed->GetY()*seed->GetY()), xToGo));
#endif
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
  BringTo02Pi(fTrImpData[kTrPhiIn]);
  double dr  = lrA->GetDR();                              // approximate X dist at the inner radius
  if (!sc.GetXYZAt(sc.GetX()-dr, bz, fTrImpData + kTrXOut)) {
    // special case: track does not reach inner radius, might be tangential
    double r = sc.GetD(0,0,bz);
    double x;
    if (!sc.GetXatLabR(r,x,bz,-1)) return kFALSE;
    dr = Abs(sc.GetX() - x);
    if (!sc.GetXYZAt(x, bz, fTrImpData + kTrXOut)) return kFALSE;
  }
  //
  fTrImpData[kTrPhiOut] = ATan2(fTrImpData[kTrYOut],fTrImpData[kTrXOut]);
  BringTo02Pi(fTrImpData[kTrPhiOut]);
  double sgy = sc.GetSigmaY2() + dr*dr*sc.GetSigmaSnp2() + AliITSUReconstructor::GetRecoParam()->GetSigmaY2(ilrA);
  double sgz = sc.GetSigmaZ2() + dr*dr*sc.GetSigmaTgl2() + AliITSUReconstructor::GetRecoParam()->GetSigmaZ2(ilrA);
  sgy = Sqrt(sgy)*fCurrTrackCond->GetNSigmaRoadY(ilrA);
  sgz = Sqrt(sgz)*fCurrTrackCond->GetNSigmaRoadZ(ilrA);
  double phi0  = MeanPhiSmall(fTrImpData[kTrPhiOut],fTrImpData[kTrPhiIn]);
  double dphi0 = DeltaPhiSmall(fTrImpData[kTrPhiOut],fTrImpData[kTrPhiIn]);
  //
  fTrImpData[kTrPhi0] = phi0;
  fTrImpData[kTrZ0]   = 0.5*(fTrImpData[kTrZOut]+fTrImpData[kTrZIn]);
  dphi0 += sgy/lrA->GetR();
  fTrImpData[kTrDPhi] =  dphi0<PiOver2() ? dphi0 : PiOver2();
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
  //  AliCluster *cl = fITS->GetLayerActive(lr)->GetCluster(clID);
  AliCluster *cl = fCurrLayer->GetCluster(clID);
  //
  Bool_t goodCl = kFALSE;
  int currLabel = Abs(fCurrESDtrMClb);
  //
  if (cl->GetLabel(0)>=0) {
    for (int i=0;i<3;i++) if (cl->GetLabel(i)>=0 && cl->GetLabel(i)==currLabel) {goodCl = kTRUE; break;}
  }
  //
  if (TMath::Abs(cl->GetX())>kTolerX) { // if due to the misalingment X is large, propagate track only
    if (!track->PropagateParamOnlyTo(track->GetX()+cl->GetX(),GetBz())) {
      //
#ifdef _ITSU_DEBUG_
      if (AliDebugLevelClass()>2 && goodCl && !track->ContainsFake()) {
	AliDebug(2,Form("Lost good cl on L:%d failed propagation. |ESDtrack#%d (MClb:%d)",lr,fCurrESDtrack->GetID(),fCurrESDtrMClb)); 
	track->Print("etp");
	cl->Print();
      }
#endif
      return kStopSearchOnSensor; // propagation failed, seedT is intact
    }
  }
  double dy = cl->GetY()-track->GetY();
  double dz = cl->GetZ()-track->GetZ();
  //
  double dy2 = dy*dy;
  double tol2 = (track->GetSigmaY2() + AliITSUReconstructor::GetRecoParam()->GetSigmaY2(lr))*
    fCurrTrackCond->GetNSigmaRoadY(lr)*fCurrTrackCond->GetNSigmaRoadY(lr); // RS TOOPTIMIZE
  if (dy2>tol2) {                          // the clusters are sorted in Z(col) then in Y(row). 
    //
#ifdef _ITSU_DEBUG_
    if (AliDebugLevelClass()>2 && goodCl &&  !track->ContainsFake()) {    
      AliDebug(2,Form("Lost good cl: dy2=%e > tol2=%e |ESDtrack#%d (MClb:%d)",dy2,tol2,fCurrESDtrack->GetID(),fCurrESDtrMClb)); 
      track->Print("etp");
      cl->Print();
      PrintSeedClusters(track);
    }
#endif
    // clusters are sorted in increasing Y and the loop goes from last (highers Y) to 1st.
    // if the track has too large y for this cluster (dy<0) it will be so for all other clusters of the sensor
    if (dy<0) return kStopSearchOnSensor;  // No chance that other cluster of this sensor will match (all Y's will be even larger)
    else      return kClusterNotMatching;   // Other clusters may match
  }
  double dz2 = dz*dz;
  tol2 = (track->GetSigmaZ2() + AliITSUReconstructor::GetRecoParam()->GetSigmaZ2(lr))*
    fCurrTrackCond->GetNSigmaRoadZ(lr)*fCurrTrackCond->GetNSigmaRoadZ(lr); // RS TOOPTIMIZE
  if (dz2>tol2) {
#ifdef _ITSU_DEBUG_
    if (AliDebugLevelClass()>2 && goodCl &&  !track->ContainsFake()) {
      AliDebug(2,Form("Lost good cl on L:%d : dz2=%e > tol2=%e |ESDtrack#%d (MClb:%d)",lr,dz2,tol2,fCurrESDtrack->GetID(),fCurrESDtrMClb)); 
      track->Print("etp");
      cl->Print();
      PrintSeedClusters(track);
    }
#endif
    return kClusterNotMatching; // Other clusters may match
  }
  //
  // check chi2
  Double_t p[2]={cl->GetY(), cl->GetZ()};
  Double_t cov[3]={cl->GetSigmaY2(), cl->GetSigmaYZ(), cl->GetSigmaZ2()};
  double chi2 = track->GetPredictedChi2(p,cov);
  //
  if (chi2>fCurrTrackCond->GetMaxTr2ClChi2(lr)) {
#ifdef _ITSU_DEBUG_
    if (AliDebugLevelClass()>2 && goodCl &&  !track->ContainsFake()) {
      AliDebug(2,Form("Lost good cl on L:%d : Chi2=%e > Chi2Max=%e |dy: %+.3e dz: %+.3e |ESDtrack#%d (MClb:%d)\n",
		      lr,chi2,fCurrTrackCond->GetMaxTr2ClChi2(lr),dy,dz,fCurrESDtrack->GetID(),fCurrESDtrMClb)); 
      track->Print("etp");
      cl->Print("");
      PrintSeedClusters(track);
    }
#endif
    return kClusterNotMatching;
  }
  //
  track = NewSeedFromPool(track);  // input track will be reused, use its clone for updates
  if (!track->Update()) {
#ifdef _ITSU_DEBUG_
    if (AliDebugLevelClass()>2 && goodCl &&  !track->ContainsFake()) {
      AliDebug(2,Form("Lost good cluster on L:%d : failed to update",lr));
      track->Print("etp");
      cl->Print("");
      PrintSeedClusters(track);
    }
#endif
    MarkSeedFree(track);
    return kClusterNotMatching;
  }
  track->SetChi2Cl(chi2);
  track->SetLrClusterID(lr,clID);
  //
  if (track->GetChi2GloNrm()>fCurrTrackCond->GetMaxChi2GloNrm(lr)) {
#ifdef _ITSU_DEBUG_
    if (AliDebugLevelClass()>2 && goodCl &&  !track->ContainsFake()) {
      AliDebug(2,Form("Lost good cl on L:%d : Chi2Glo=%e > Chi2Max=%e |dy: %+.3e dz: %+.3e |ESDtrack#%d (MClb:%d)\n",
		      lr,track->GetChi2GloNrm(),fCurrTrackCond->GetMaxChi2GloNrm(lr),dy,dz,fCurrESDtrack->GetID(),fCurrESDtrMClb)); 
      track->Print("etp");
      cl->Print("");
      PrintSeedClusters(track);
    }
#endif
    MarkSeedFree(track);
    return kClusterNotMatching;
  }
  //  cl->IncreaseClusterUsage(); // do this only for winners
  //
  track->SetFake(!goodCl);
  //
#ifdef _ITSU_DEBUG_
  AliDebug(2,Form("AddCl(%d) Cl%d lr:%d: dY:%+8.4f dZ:%+8.4f (MC: %5d %5d %5d) |Chi2=%f(%c)| ",
	       goodCl,clID,lr,dy,dz2,cl->GetLabel(0),cl->GetLabel(1),cl->GetLabel(2), chi2, track->IsFake() ? '-':'+'));
#endif
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
    if (lastChecked) patt |= (~(kMask<<lastChecked))&fCurrTrackCond->GetAllowedLayers(); // not all layers were checked, complete unchecked ones by potential hits
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
  Double_t xyz0[3],xyz1[3];
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
	Double_t xrho,xx0;
	ds = GetMaterialBudget(xyz0,xyz1,xx0,xrho);	
	if (dir>0) xrho = -xrho; // outward should be negative
	if (!seed->ApplyMaterialCorrection(xx0,xrho,mass,kFALSE)) {AliInfo("Fail2"); seed->Print("etp"); return kFALSE;}
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
  Double_t xyz0[3],xyz1[3];
  //
#ifdef _ITSU_DEBUG_
  if (AliDebugLevelClass()>1) {
    AliDebug(2,Form("Before Propagate to X=%f with M=%.3f MaxStep=%.4f MatCorr=%d",xToGo,mass,maxStep,matCorr));
    seed->AliExternalTrackParam::Print();
  }
#endif
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
	Double_t xrho,xx0;
	ds = GetMaterialBudget(xyz0,xyz1,xx0,xrho);	
	if (dir>0) xrho = -xrho; // outward should be negative
	if (!seed->CorrectForMeanMaterial(xx0,xrho,mass)) return kFALSE;
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
#ifdef _ITSU_DEBUG_
  if (AliDebugLevelClass()>1) {
    AliDebug(2,Form("After Propagate to X=%f",xToGo));
    seed->AliExternalTrackParam::Print();
  }
#endif
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliITSUTrackerGlo::FinalizeHypothesis(AliITSUTrackHyp* hyp)
{
  // finalize single TPC track hypothesis
  if (!hyp || hyp->GetESDTrack()->IsOn(AliESDtrack::kITSin)) return kFALSE;
  AliITSUSeed* winner = 0;
  fCurrHyp = hyp;
  fCurrMass = hyp->GetMass();
  if (!(winner=hyp->DefineWinner(fCurrTrackCond->GetActiveLrInner()))) return kFALSE;
  CookMCLabel(hyp);
  //
#ifdef  _ITSU_TUNING_MODE_  // fill tuning histos
  TObjArray* dest = hyp->GetLabel()>=0 ? fCHistoArrCorr : fCHistoArrFake;
  if (dest) {
    AliITSUSeed* sd = winner;
    double htrPt = hyp->Pt();
    do {
      int clID,lrID;
      if ( (clID=sd->GetLrCluster(lrID))<0 ) continue;
      ((TH2*)dest->At( GetHistoID(lrID,kHChi2Nrm     ,fCurrPassID) ))->Fill(htrPt,sd->GetChi2GloNrm());
      ((TH2*)dest->At( GetHistoID(lrID,kHBestInBranch,fCurrPassID) ))->Fill(htrPt,sd->GetOrdBranch());
      ((TH2*)dest->At( GetHistoID(lrID,kHBestInCand  ,fCurrPassID) ))->Fill(htrPt,sd->GetOrdCand());
      //
      if (dest==fCHistoArrFake && !sd->IsFake()) continue; // for the fake seeds fill only fake clusters part
      //
      ((TH2*)dest->At( GetHistoID(lrID,kHResY        ,fCurrPassID) ))->Fill(htrPt,sd->GetResidY());
      ((TH2*)dest->At( GetHistoID(lrID,kHResZ        ,fCurrPassID) ))->Fill(htrPt,sd->GetResidZ());
      ((TH2*)dest->At( GetHistoID(lrID,kHResYP       ,fCurrPassID) ))->Fill(htrPt,sd->GetPullY());
      ((TH2*)dest->At( GetHistoID(lrID,kHResZP       ,fCurrPassID) ))->Fill(htrPt,sd->GetPullZ());
      ((TH2*)dest->At( GetHistoID(lrID,kHChi2Cl      ,fCurrPassID) ))->Fill(htrPt,sd->GetChi2Cl());
      //
    } while((sd=(AliITSUSeed*)sd->GetParent()));
    //
    ((TH2*)dest->At( GetHistoID(-1,kHChiMatch,fCurrPassID) ))->Fill(htrPt,winner->GetChi2ITSTPC());
    ((TH2*)dest->At( GetHistoID(-1,kHChiITSSA,fCurrPassID) ))->Fill(htrPt,winner->GetChi2ITSSA());
  }
  //
#endif
  //
  CheckClusterSharingConflicts(hyp);
  return hyp->GetWinner();  // winner might change of disappear after resolving conflicts
  //
}

//______________________________________________________________________________
void AliITSUTrackerGlo::FinalizeHypotheses()
{
  // select winner for each hypothesis, remove cl. sharing conflicts
  AliCodeTimerAuto("",0); 
  //
  int* esdTrackIndex = fESDIndex.GetArray();
  for (int ih=0;ih<fNTracksESD;ih++) FinalizeHypothesis(GetTrackHyp(esdTrackIndex[ih])); // finlize in decreasing pt order
  //
  for (int ih=fNTracksESD;ih--;) UpdateESDTrack(GetTrackHyp(esdTrackIndex[ih]),AliESDtrack::kITSin); 
  //
}

//______________________________________________________________________________
void AliITSUTrackerGlo::CheckClusterSharingConflicts(AliITSUTrackHyp* hyp)
{
  // remove eventual cluster sharing conflicts
  AliITSUSeed* winner = 0;
  if (!(winner=hyp->GetWinner())) return;
  UShort_t idH = (UShort_t)hyp->GetUniqueID()+1;
  int lrID,clID;
  int inLr = fCurrTrackCond->GetActiveLrInner();
  AliITSUSeed* winSD=winner;
  do {
    if ( (clID=winSD->GetLrCluster(lrID))<0 ) continue;
    AliITSUClusterPix* cl = (AliITSUClusterPix*)fITS->GetLayerActive(lrID)->GetCluster(clID);
    Int_t refID = cl->GetRecoInfo();                // was it already referred by some track (with refID-1 if >0)
    if (!refID) {cl->SetRecoInfo(idH); continue;}   // if not, refer from track to cluster
    if (refID==idH) continue;                          // ignore reference to itself
    // the cluster is already used by other track, need to resolve the conflict
    AliITSUTrackHyp* hypC = GetTrackHyp(refID-1);   // competitor
    //    printf("Ref to %d (%p) from %d Cl %d %d\n",refID-1,hypC, idH-1,clID,lrID);
    AliITSUSeed* winnerC = hypC->GetWinner();
    if (!winnerC) {
      printf("Missing winner, ERROR\n");
      continue;
    }
    //
#ifdef _ITSU_DEBUG_
    if (AliDebugLevelClass()>1) { //------------------------ DEBUG -------------------
      AliInfo(Form("Shared cl#%4d on Lr:%d: Hyp%3d/Q:%6.3f vs Hyp%3d/Q:%6.3f",
		   clID,lrID,idH-1,winner->GetQualityVar(),refID-1,winnerC->GetQualityVar()));
      // dump winner of hyp
      printf("#%4d Pt:%.3f Chi:%6.3f/%6.3f/%6.3f Ncl:%d MCits%+5d MCtpc:%+5d |",
	     idH-1,winner->Pt(),winner->GetChi2GloNrm(),winner->GetChi2ITSTPC(),winner->GetChi2ITSSA(),
	     winner->GetNLayersHit(),hyp->GetITSLabel(),hyp->GetESDTrack()->GetTPCLabel());
      int prevL=-1;
      AliITSUSeed* sd = winner;
      do {
	int lrs;
	int clIDt = sd->GetLrCluster(lrs);
	if (clIDt<0) continue;
	while( lrs>++prevL ) printf("%4s        ","----");
	printf("%4d (%5.1f)",clIDt,sd->GetChi2Cl());
      } while ((sd=(AliITSUSeed*)sd->GetParent()));
      printf("|\n");
      // dump winner of hypC
      printf("#%4d Pt:%.3f Chi:%6.3f/%6.3f/%6.3f Ncl:%d MCits%+5d MCtpc:%+5d |",
	     refID-1,winnerC->Pt(),winnerC->GetChi2GloNrm(),winnerC->GetChi2ITSTPC(),winnerC->GetChi2ITSSA(),
	     winnerC->GetNLayersHit(),hypC->GetITSLabel(),hypC->GetESDTrack()->GetTPCLabel());
      prevL=-1;
      sd = winnerC;
      do {
	int lrs;
	int clIDt = sd->GetLrCluster(lrs);
	if (clIDt<0) continue;
	while( lrs>++prevL ) printf("%4s        ","----");
	printf("%4d (%5.1f)",clIDt,sd->GetChi2Cl());
      } while ((sd=(AliITSUSeed*)sd->GetParent()));
      printf("|\n");
    }  //-------------------------------------------------- DEBUG -------------------
#endif
    //
    if (winner->GetQualityVar()<winnerC->GetQualityVar()) { // current tracks is better then competitor track
      FlagSeedClusters(winnerC,kFALSE,refID); // unflag cluster usage by loser
      cl->SetRecoInfo(idH);
      //
      //if ( hypC->GetLabel()>=0) AliInfo(Form("KILLING CORRECT: %d",refID-1));
      winnerC->Kill();
      if (hypC->DefineWinner(inLr)) {             // find new winner instead of suppressed candidate
	CookMCLabel(hypC);
	CheckClusterSharingConflicts(hypC);   // and check its sharing conflicts again
	if (winner->IsKilled()) break;        // the current winner might have been killed during check of new winner of hypC!
      }
    }    
    else { // competitor hypC is better than the hyp
      FlagSeedClusters(winner,kFALSE,idH); // unflag cluster usage by loser
      winner->Kill();
      //if ( hyp->GetLabel()>=0) AliInfo(Form("KILLING CORRECT: %d",idH-1));
      if (hyp->DefineWinner(inLr)) {
	CookMCLabel(hyp);
	CheckClusterSharingConflicts(hyp);
      }
      break;
    }
    //   
  } while ((winSD=(AliITSUSeed*)winSD->GetParent()));
  //
  return;
}

//______________________________________________________________________________
void AliITSUTrackerGlo::UpdateESDTrack(AliITSUTrackHyp* hyp,Int_t flag)
{
  // update ESD track with current best hypothesis
  if (!hyp) return;
  AliESDtrack* esdTr = hyp->GetESDTrack();
  if (!esdTr || esdTr->IsOn(flag)) return;
  AliITSUSeed* win = hyp->GetWinner();
  if (!win || win->IsKilled()) return;
  double chiSav;
  //
  switch (flag) {
  case AliESDtrack::kITSin: 
    esdTr->UpdateTrackParams(hyp,flag); // update kinematics
    fCountITSin++;
    // set fakes cluster info
    {
      UShort_t clfake = 0;
      Int_t    clSplit = 0;
      AliITSUSeed* sd = win;
      int ip=0;
      do {
	int lr, clID = sd->GetLrCluster(lr);
	//	if (clID<0 && sd->IsFake()) {printf("Fake non-cluster at lr%d",lr); win->Print("parent etp");}
	if (clID>=0) {
	  if (sd->IsFake()) clfake |= 0x1<<lr;
	  esdTr->SetITSModuleIndex(ip++, sd->GetLrClusterID());
	  AliITSUClusterPix *cl = (AliITSUClusterPix*)fITS->GetLayerActive(lr)->GetCluster(clID);
#ifdef  _ITSU_TUNING_MODE_
	  if (cl->IsSplit()) clSplit |= 0x1<<lr;
#endif
	}
      } while ((sd=(AliITSUSeed*)sd->GetParent()));
      //
      // RS: Temporary set special flag for tracks from the afterburner
      if (fCurrPassID>0) clfake |= 0x1<<7;
      //
      esdTr->SetITSSharedMap(clfake);
      esdTr->SetITSModuleIndex(10,clSplit);
    }
    // TEMPORARY: store iteration id
    esdTr->SetITSModuleIndex(11,fCurrPassID);
    break;
    //
  case AliESDtrack::kITSout: 
    // here the stored chi2 will correspond to backward ITS-SA fit
    esdTr->UpdateTrackParams(hyp,flag); // update kinematics
    fCountITSout++;
    // TODO: avoid setting friend
    break;
    //
  case AliESDtrack::kITSrefit: 

    // at the moment fill as a chi2 the TPC-ITS matching chi2
    chiSav = hyp->GetChi2();
    hyp->SetChi2(win->GetChi2ITSTPC());
    esdTr->UpdateTrackParams(hyp,flag); // update kinematics
    hyp->SetChi2(chiSav);
    fCountITSrefit++;
    // TODO: avoid setting cluster info
    break;
  default:
    AliFatal(Form("Unknown flag %d",flag));
  }
  //
  esdTr->SetITSLabel(hyp->GetITSLabel());
  // transfer chip indices
  // TODO
}

//______________________________________________________________________________
Double_t AliITSUTrackerGlo::RefitTrack(AliITSUTrackHyp* trc, Double_t rDest, Int_t &nclFit, Int_t stopCond)
{
  // refit track till radius rDest. The cluster,mass info is taken from fCurrHyp (and its winner seed)
  // if stopCond<0 : propagate till last cluster then stop
  // if stopCond==0: propagate till last cluster then try to go till limiting rDest, don't mind if fail
  // if stopCond>0 : rDest must be reached
  //
  double rCurr = Sqrt(trc->GetX()*trc->GetX() + trc->GetY()*trc->GetY());
  int dir,lrStart,lrStop;
  //
  dir = rCurr<rDest ? 1 : -1;
  lrStart = fITS->FindFirstLayerID(rCurr,dir);
  lrStop  = fITS->FindLastLayerID(rDest,dir); // lr id before which we have to stop
  //
#ifdef _ITSU_DEBUG_
  if (AliDebugLevelClass()>2) {
    printf("Refit %d: Lr:%d (%f) -> Lr:%d (%f)\n",dir,lrStart,rCurr, lrStop,rDest);
    printf("Before refit: "); trc->Print();
  }
#endif
  if (lrStop<0 || lrStart<0) AliFatal(Form("Failed to find start(%d) or last(%d) layers. Track from %.3f to %.3f",lrStart,lrStop,rCurr,rDest));
  //
  Int_t clInfo[2*AliITSMFTAux::kMaxLayers];
  Int_t nCl = fCurrHyp->FetchClusterInfo(clInfo);
  fCurrMass = fCurrHyp->GetMass();
  AliITSUTrackHyp tmpTr(*(AliKalmanTrack*)trc);
  double chi2 = 0;
  int iclLr[2],nclLr;
  nclFit = 0;
  //
  int lrStop1 = lrStop+dir;
  for (int ilr=lrStart;ilr!=lrStop1;ilr+=dir) {
    AliITSURecoLayer* lr = fITS->GetLayer(ilr);
    if ( dir*(rCurr-lr->GetR(dir))>0) continue; // this layer is already passed
    int ilrA2,ilrA = lr->GetActiveID();
    // passive layer or active w/o hits will be traversed on the way to next cluster
    if (!lr->IsActive() || clInfo[ilrA2=(ilrA<<1)]<0) continue; 
    //
    // select the order in which possible 2 clusters (in case of the overlap) will be traversed and fitted
    nclLr=0;
    if (dir>0) { // clusters are stored in increasing radius order
      iclLr[nclLr++]=clInfo[ilrA2++];
      if (clInfo[ilrA2]>=0) iclLr[nclLr++]=clInfo[ilrA2];
    }
    else {
      if ( clInfo[ilrA2+1]>=0 ) iclLr[nclLr++]=clInfo[ilrA2+1];
      iclLr[nclLr++]=clInfo[ilrA2];
    }
    //
    Bool_t transportedToLayer = kFALSE;
    for (int icl=0;icl<nclLr;icl++) {
      AliITSUClusterPix* clus =  (AliITSUClusterPix*)lr->GetCluster(iclLr[icl]);
      AliITSURecoSens* sens = lr->GetSensorFromID(clus->GetVolumeId());
      if (!tmpTr.Rotate(sens->GetPhiTF())) {
#ifdef _ITSU_DEBUG_
	AliESDtrack* trESD = fCurrHyp->GetESDTrack();
	AliDebug(2,Form("Failed on rotate to %f | ESDtrack#%d (MClb:%d)",sens->GetPhiTF(),trESD->GetID(),trESD->GetTPCLabel()));
#endif
	return -1;
      }
      //
      double xClus = sens->GetXTF()+clus->GetX();
      if (!transportedToLayer) {
	if (ilr!=lrStart && !TransportToLayerX(&tmpTr,lrStart,ilr,xClus)) {
#ifdef _ITSU_DEBUG_
	  AliESDtrack* trESD = fCurrHyp->GetESDTrack();
	  AliDebug(2,Form("Failed to transport %d -> %d | ESDtrack#%d (MClb:%d)\n",lrStart,ilr,trESD->GetID(),trESD->GetTPCLabel()));
#endif
	  return -1; // go to the entrance to the layer
	}
	lrStart = ilr;
	transportedToLayer = kTRUE;
      }
      //
#ifdef _ITSU_DEBUG_
      if (AliDebugLevelClass()>1) {
	AliDebug(2,Form("Propagate to cl:%d on lr %d Need to go %.4f -> %.4f",icl,ilrA,tmpTr.GetX(),xClus));
      }
#endif
      //
      if (!PropagateSeed(&tmpTr,xClus,fCurrMass)) {
#ifdef _ITSU_DEBUG_
	AliESDtrack* trESD = fCurrHyp->GetESDTrack();
	AliDebug(2,Form("Failed on propagate to %f (dir=%d) | ESDtrack#%d (MClb:%d)",xClus,dir,trESD->GetID(),trESD->GetTPCLabel()));
#endif
	return -1;
      }
      //
      Double_t p[2]={clus->GetY(), clus->GetZ()};
      Double_t cov[3]={clus->GetSigmaY2(), clus->GetSigmaYZ(), clus->GetSigmaZ2()};
      double chi2cl = tmpTr.GetPredictedChi2(p,cov);
      chi2 += chi2cl;
      //
#ifdef  _ITSU_TUNING_MODE_
      TObjArray* dest = trc->GetLabel()>=0 ? fCHistoArrCorr : fCHistoArrFake;
      if (dest && fTrackPhaseID>kClus2Tracks) {
	//
	double htrPt = tmpTr.Pt();
	double dy = p[0]-tmpTr.GetY();
	double dz = p[1]-tmpTr.GetZ();
	((TH2*)dest->At( GetHistoID(ilrA,kHResY,0,fTrackPhaseID) ))->Fill(htrPt,dy);
	((TH2*)dest->At( GetHistoID(ilrA,kHResZ,0,fTrackPhaseID) ))->Fill(htrPt,dz);
	double errY = tmpTr.GetSigmaY2();
	double errZ = tmpTr.GetSigmaZ2();
	if (errY>0) ((TH2*)dest->At( GetHistoID(ilrA,kHResYP,0,fTrackPhaseID) ))->Fill(htrPt,dy/Sqrt(errY));
	if (errZ>0) ((TH2*)dest->At( GetHistoID(ilrA,kHResZP,0,fTrackPhaseID) ))->Fill(htrPt,dz/Sqrt(errZ));
	((TH2*)dest->At( GetHistoID(ilrA,kHChi2Cl,0,fTrackPhaseID) ))->Fill(htrPt,chi2cl);
      }
#endif  
      //      
      if ( !tmpTr.Update(p,cov) ) {
#ifdef _ITSU_DEBUG_
	AliESDtrack* trESD = fCurrHyp->GetESDTrack();
	AliDebug(2,Form("Failed on Update | ESDtrack#%d (MClb:%d)",trESD->GetID(),trESD->GetTPCLabel()));
#endif
	return -1;
      }
#ifdef _ITSU_DEBUG_
      if (AliDebugLevelClass()>1) {
	printf("AfterRefit: "); tmpTr.Print();
      }
#endif
      if (++nclFit==nCl && stopCond<0) {
	*trc = (AliKalmanTrack&)tmpTr;
	return chi2; // it was requested to not propagate after last update
      }
    }
    //
  }
  // All clusters were succesfully fitted. Even if the track does not reach rDest, this is enough to validate it.
  // Still, try to go as close as possible to rDest.
  //
  if (lrStart!=lrStop) {
    if (!TransportToLayer(&tmpTr,lrStart,lrStop)) {
#ifdef _ITSU_DEBUG_
      AliDebug(1,Form("Failed on TransportToLayer %d->%d | ESDtrack#%d (MClb:%d), X=%f",lrStart,lrStop,fCurrESDtrack->GetID(),fCurrESDtrMClb,tmpTr.GetX()));
#endif
      return (stopCond>0) ? -chi2 : chi2; // rDest was obligatory
    }    
    if (!GoToExitFromLayer(&tmpTr,fITS->GetLayer(lrStop),dir)) {
#ifdef _ITSU_DEBUG_
      AliDebug(1,Form("Failed on GoToExitFromLayer %d | ESDtrack#%d (MClb:%d), X=%f",lrStop,fCurrESDtrack->GetID(),fCurrESDtrMClb,tmpTr.GetX()));
#endif
      return (stopCond>0) ? -chi2 : chi2; // rDest was obligatory
    }
  }
  // go to the destination radius. Note that here we don't select direction to avoid precision problems
  if (!tmpTr.GetXatLabR(rDest,rDest,GetBz(),0) || !PropagateSeed(&tmpTr,rDest,fCurrMass, 100, kFALSE)) {
    return (stopCond>0) ? -chi2 : chi2; // rDest was obligatory
  }
  *trc = (AliKalmanTrack&)tmpTr;
#ifdef _ITSU_DEBUG_
  if (AliDebugLevelClass()>2) {
    printf("After refit (now at lr %d): ",lrStop); trc->Print();
  }
#endif
  return chi2;
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
  AliESDtrack* esdTr = hyp->GetESDTrack();
  int tpcLab = esdTr ? Abs(esdTr->GetTPCLabel()) : -kDummyLabel;
  if (nCl && nLab) {
    int maxLab=0,nTPCok=0;
    for (int ilb=nLab;ilb--;) {
      int st = lbStat[ilb];
      if (lbStat[maxLab]<st) maxLab=ilb;
      if (tpcLab==lbID[ilb]) nTPCok += st;
    }
    hyp->SetFakeRatio(1.-float(nTPCok)/nCl);
    hyp->SetLabel( nTPCok==nCl ? tpcLab : -tpcLab);
    hyp->SetITSLabel( lbStat[maxLab]==nCl ? lbID[maxLab] : -lbID[maxLab]); // winner label
    return;
  }
  //
  hyp->SetFakeRatio(-1.);
  hyp->SetLabel( -tpcLab );
  hyp->SetITSLabel( kDummyLabel );
}

//__________________________________________________________________
Bool_t AliITSUTrackerGlo::AddSeedBranch(AliITSUSeed* seed)
{
  // add new prolongation branch for the seed to temporary array. The seeds are sorted according to Compare f-n
  if (fNCandidatesAdded+fNBranchesAdded>=fLayerMaxCandidates ) {
    AliInfo(Form("Number of candidates at layer %d for current seed reached %d, increasing buffer",fCurrActLrID,fLayerMaxCandidates));
    fLayerMaxCandidates = 2*(fLayerMaxCandidates+1);
    AliITSUSeed** tmpArr = fLayerCandidates;
    fLayerCandidates = new AliITSUSeed*[fLayerMaxCandidates];
    memcpy(fLayerCandidates,tmpArr,(fNCandidatesAdded+fNBranchesAdded)*sizeof(AliITSUSeed*));
    delete[] tmpArr; // delete only array, not objects
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
  AliCodeTimerAuto("",0); 
  int nb = Min(fNBranchesAdded,acceptMax);
  //  if (nb<fNBranchesAdded) printf("ValidateAllowedBranches: %d of %d (%d) on lr %d\n",nb,fNBranchesAdded,fNCandidatesAdded,ilr);
  // disable unused branches
  AliITSUSeed** branches = &fLayerCandidates[fNCandidatesAdded];
  //
#ifdef  _ITSU_TUNING_MODE_
  for (int ib=0;ib<fNBranchesAdded;ib++) branches[ib]->SetOrdBranch(ib);
#endif
  //
  for (int ib=nb;ib<fNBranchesAdded;ib++) MarkSeedFree(branches[ib]);
  fNCandidatesAdded += nb; // update total candidates counter
  fNBranchesAdded = 0; // reset branches counter
  //
}

//__________________________________________________________________
void AliITSUTrackerGlo::ValidateAllowedCandidates(Int_t ilr, Int_t acceptMax)
{
  // transfer allowed number of branches to hypothesis container
  //
  AliCodeTimerAuto("",0); 
  // sort candidates in increasing order of chi2
  static int lastSize = 0;
  static int *index = 0;
  static float *chi2 = 0;
  if (fLayerMaxCandidates>lastSize) {
    lastSize = fLayerMaxCandidates;
    delete[] index;
    delete[] chi2;
    index = new int[lastSize];
    chi2  = new float[lastSize];
  }
  for (int i=fNCandidatesAdded;i--;) chi2[i] = fLayerCandidates[i]->GetChi2GloNrm();
  Sort(fNCandidatesAdded,chi2,index,kFALSE);
  //
  int nacc=0,nb=0;
  if (ilr>fCurrTrackCond->GetActiveLrInner()) { // just take 1st acceptMax candidates
    nb = Min(fNCandidatesAdded,acceptMax);
    for (int ib=0;ib<nb;ib++) AddProlongationHypothesis(fLayerCandidates[index[ib]],ilr);
  }
  else { // for innermost layer test ITS SA backward chi2 and matching to TPC
    AliITSUSeed* wn0 = fCurrHyp->GetWinner(); // in principle, should be NULL
    for (nb=0;nb<fNCandidatesAdded;nb++) {
      AliITSUSeed* sd = fLayerCandidates[index[nb]];    
      fCurrHyp->SetWinner(sd);
      if (!CheckBackwardMatching(sd)) MarkSeedFree(sd); // discard bad backward tracks
      else {
	AddProlongationHypothesis(sd,ilr);
	if (++nacc==acceptMax) {nb++; break;} // the rest will be discarded
      }
    }
    fCurrHyp->SetWinner(wn0); // restore original winner (NULL?)
  }
  // disable unused candidates
  for (int ib=nb;ib<fNCandidatesAdded;ib++) MarkSeedFree(fLayerCandidates[index[ib]]);    
  fNCandidatesAdded = 0; // reset candidates counter
  //
}

//__________________________________________________________________
Bool_t AliITSUTrackerGlo::CheckBackwardMatching(AliITSUSeed* seed)
{
  // check seed backward propagation chi2 and matching to TPC 
  double bz0 = GetBz();
  double rDest = fITS->GetRITSTPCRef(); // reference radius for comparison
  AliITSUTrackHyp trback;
  trback.AliExternalTrackParam::operator=(*seed);
  trback.ResetCovariance(10000);
  int nclFit = 0;
  double chi2sa = RefitTrack(&trback,rDest,nclFit,1);
  if (chi2sa<0) return kFALSE;
  int ndf = nclFit*2-5;
  if (ndf>0) chi2sa /= ndf;
  if (chi2sa>fCurrTrackCond->GetMaxITSSAChi2(nclFit)) return kFALSE;
  //
  // relate to TPC track at outer layer
  AliExternalTrackParam* tpcSeed = fCurrHyp->GetTPCSeed();
  if (!trback.Rotate(tpcSeed->GetAlpha()) || !trback.PropagateParamOnlyTo(tpcSeed->GetX(),bz0)) {
#ifdef _ITSU_DEBUG_
    if (AliDebugLevelClass()>+1 && !seed->ContainsFake()) {
      AliInfo(Form("Failed to propagate seed to TPC track @ X:%.3f Alpha:%.3f",tpcSeed->GetX(),tpcSeed->GetAlpha()));
      seed->Print("etp");
      trback.Print();
    }
#endif
    return kFALSE; 
  }
  double chi2Match = trback.GetPredictedChi2(tpcSeed)/5;
  if (chi2Match>fCurrTrackCond->GetMaxITSTPCMatchChi2()) return kFALSE;
  //
  seed->SetChi2ITSSA(chi2sa);
  seed->SetChi2ITSTPC(chi2Match);
  return kTRUE;
}

//__________________________________________________________________
void AliITSUTrackerGlo::SaveReducedHypothesesTree(AliITSUTrackHyp* dest)
{
  // remove those hypothesis seeds which dont lead to candidates at final layer
  //
  // 1st, flag the seeds to save
  int lr0 = fCurrTrackCond->GetActiveLrInner();
  for (int isd=0;isd<fCurrHyp->GetNSeeds(lr0);isd++) {
    AliITSUSeed* seed = fCurrHyp->RemoveSeed(lr0,isd);
    if (!seed) continue;
    seed->FlagTree(AliITSUSeed::kSave);
    dest->AddSeed(seed,lr0);
  }
  for (int ilr=0;ilr<fNLrActive;ilr++) {
    int nsd = fCurrHyp->GetNSeeds(ilr);
    for (int isd=0;isd<nsd;isd++) {
      AliITSUSeed* seed = fCurrHyp->RemoveSeed(ilr,isd);
      if (!seed) continue; // already discarded or saved
      if (seed->IsSaved()) dest->AddSeed(seed,ilr);
      else MarkSeedFree(seed);
    }
  }
  //
  //  AliInfo(Form("SeedsPool: %d, BookedUpTo: %d, free: %d",fSeedsPool.GetSize(),fSeedsPool.GetEntriesFast(),fNFreeSeeds));
}

//__________________________________________________________________
void AliITSUTrackerGlo::CleanHypothesis(AliITSUTrackHyp* hyp)
{
  // clean hypothesis
  hyp->SetWinner(0);
  for (int ilr=0;ilr<fNLrActive;ilr++) {
    int nsd = hyp->GetNSeeds(ilr);
    for (int isd=0;isd<nsd;isd++) {
      AliITSUSeed* seed = hyp->RemoveSeed(ilr,isd);
      if (!seed) continue; // already discarded or saved
      MarkSeedFree(seed);
    }
  }
}

//__________________________________________________________________
void AliITSUTrackerGlo::FlagSplitClusters()
{
  // set special bit on split clusters using MC info
  for (int ilr=fNLrActive;ilr--;) {
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

//__________________________________________________________________
void AliITSUTrackerGlo::FlagSeedClusters(const AliITSUSeed* seed, Bool_t flg, UShort_t ref)
{
  // mark used clusters
  int lrID,clID;
  while (seed) {
    if ( (clID=seed->GetLrCluster(lrID))>=0 ) {
      AliITSUClusterPix* cl = (AliITSUClusterPix*)fITS->GetLayerActive(lrID)->GetCluster(clID);
      if (ref) { // do we need to set or delete cluster->track ref?
	if (flg) {
	  if (!cl->GetRecoInfo()) cl->SetRecoInfo(ref);  // set ref only if cluster already does not refer to other track, inc.counter
	}
	else {
	  if (cl->GetRecoInfo()==ref) cl->SetRecoInfo(0); // unset reference only if it refers to ref, decrease counter	    
	}
      }
    }
    seed = (AliITSUSeed*)seed->GetParent();
  }
  //
}

//__________________________________________________________________
void AliITSUTrackerGlo::CheckClusterUsage()
{
  // check cluster usage info
  printf("ClusterUsage at pass %d\n",fCurrPassID);
  for (int ilr=0;ilr<fNLrActive;ilr++) {
    AliITSURecoLayer* lr = fITS->GetLayerActive(ilr);
    int ncl = lr->GetNClusters();
    int nclUsed = 0;
    int nclUs0CntShr1 = 0;
    int nclUs1Cnt0 = 0;
    int nclUs1Shr1 = 0;
    int nclUseL = 0;
    for (int icl=0;icl<ncl;icl++) {
      AliITSUClusterPix *cl = (AliITSUClusterPix*)lr->GetCluster(icl);
      int usc = cl->GetClUsage();
      if (usc<1) {
	if (cl->IsClusterUsed() || cl->IsClusterShared()) {
	  nclUs0CntShr1++;
	  printf("Lr%d Cl%4d: UseCnt=0 but UseFlg=%d UseShr=%d\n",ilr,icl,cl->IsClusterUsed(),cl->IsClusterShared());
	}
	continue;
      }
      nclUsed++;
      if (!cl->IsClusterUsed()) {
	nclUs1Cnt0++;
	printf("Lr%d Cl%4d: UseCnt=%d but UseFlg=%d UseShr=%d\n",ilr,icl,usc, cl->IsClusterUsed(),cl->IsClusterShared());
      }
      if (cl->IsClusterShared()) {
	nclUs1Shr1++;
	printf("Lr%d Cl%4d: UseCnt=%d but UseFlg=%d UseShr=%d\n",ilr,icl,usc, cl->IsClusterUsed(),cl->IsClusterShared());
      }
      if (usc>1) {
	nclUseL++;
	printf("Lr%d Cl%4d: UseCnt=%d!! but UseFlg=%d UseShr=%d\n",ilr,icl,usc, cl->IsClusterUsed(),cl->IsClusterShared());
      }
    }
    printf("ClUsage Lr%d: Ncl=%5d Nusd=%5d (%5.2f) Use0CS1:%4d Use1C0:%4d Use1S1:%4d UseL:%d\n",
	   ilr,ncl,nclUsed,ncl? float(nclUsed)/ncl : 0, nclUs0CntShr1,nclUs1Cnt0,nclUs1Shr1,nclUseL);
  }
}


#ifdef  _ITSU_TUNING_MODE_
//__________________________________________________________________
TObjArray* AliITSUTrackerGlo::BookControlHistos(const char* pref)
{
  // book special control histos
  TString prefS = pref;
  prefS.ToLower();
  TObjArray* dest = new TObjArray;
  dest->SetOwner(kTRUE);
  //
  const int kNResDef=7;
  const double kResDef[kNResDef]={0.05,0.05,0.3, 0.05,1,0.5,1.5};
  const double ptMax=10;
  const double plMax=10;
  const double chiMax=100;
  const int nptbins=50;
  const int nresbins=400;
  const int nplbins=50;
  const int nchbins=200;
  const int maxBr  = 15;
  const int maxCand = 200;
  TString ttl;
  for (int phase=0;phase<kNTrackingPhases;phase++) {
    for (int pass=0;pass<AliITSUReconstructor::GetRecoParam()->GetNTrackingConditions();pass++) {
      //
      for (int ilr=0;ilr<fNLrActive;ilr++) {
	//
	// ----------------- These are histos to be filled in Cluster2Tracks of each pass. 
	// PropagateBack and RefitInward will be stored among the histos of 1st pass
	if (pass>0 && phase!=kClus2Tracks) continue;
	//
	double mxdf = ilr>=kNResDef ? kResDef[kNResDef-1] : kResDef[ilr];
	ttl = Form("Pass%d_S%d_residY%d_%s",pass,phase,ilr,pref);
	TH2F* hdy = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax,nresbins,-mxdf,mxdf);
	dest->AddAtAndExpand(hdy,GetHistoID(ilr,kHResY,pass,phase));
	hdy->SetDirectory(0);
	//
	ttl = Form("Pass%d_S%d_residYPull%d_%s",pass,phase,ilr,pref);	
	TH2F* hdyp = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax,nplbins,-plMax,plMax);
	dest->AddAtAndExpand(hdyp,GetHistoID(ilr,kHResYP,pass,phase));
	hdyp->SetDirectory(0);
	//
	ttl = Form("Pass%d_S%d_residZ%d_%s",pass,phase,ilr,pref);	
	TH2F* hdz = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax,nresbins,-mxdf,mxdf);
	dest->AddAtAndExpand(hdz,GetHistoID(ilr,kHResZ,pass,phase));
	hdz->SetDirectory(0);
	//
	ttl = Form("Pass%d_S%d_residZPull%d_%s",pass,phase,ilr,pref);		
	TH2F* hdzp = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax,nplbins,-plMax,plMax);
	hdzp->SetDirectory(0);
	dest->AddAtAndExpand(hdzp,GetHistoID(ilr,kHResZP,pass,phase));
	//
	ttl = Form("Pass%d_S%d_chi2Cl%d_%s",pass,phase,ilr,pref);		
	TH2F* hchi = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax, nchbins,0.,chiMax);
	hchi->SetDirectory(0);
	dest->AddAtAndExpand(hchi,GetHistoID(ilr,kHChi2Cl,pass,phase));
	//
	// ------------------- These histos are filled for Clusters2Tracks only
	if (phase!=kClus2Tracks) continue;
	//
	ttl = Form("Pass%d_S%d_chi2Nrm%d_%s",pass,phase,ilr,pref);		
	TH2F* hchiN = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax, nchbins,0.,chiMax);
	hchiN->SetDirectory(0);
	dest->AddAtAndExpand(hchiN,GetHistoID(ilr,kHChi2Nrm,pass,phase));
	//
	ttl = Form("Pass%d_S%d_bestInBranch%d_%s",pass,phase,ilr,pref);		
	TH2* hnbr = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax, maxBr,-0.5,maxBr-0.5);
	hnbr->SetDirectory(0);
	dest->AddAtAndExpand(hnbr,GetHistoID(ilr,kHBestInBranch,pass,phase));
	//
	ttl = Form("Pass%d_S%d_bestInCands%d_%s",pass,phase,ilr,pref);		
	TH2* hncn = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax, maxCand,-0.5,maxCand-0.5);
	hncn->SetDirectory(0);
	dest->AddAtAndExpand(hncn,GetHistoID(ilr,kHBestInCand,pass,phase));
	//
      } // loop over layers
      //
      // these are special histos, filled not per layer but in the end of track fit in Clusters2Tracks in EVERY pass
      //  
      if (phase==kClus2Tracks) {
	TH2* hchiMatch = 0; 
	ttl = Form("Pass%d_Chi2Match_%s",pass,pref);
	hchiMatch = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax, nchbins,0.,chiMax);
	hchiMatch->SetDirectory(0);
	dest->AddAtAndExpand(hchiMatch,GetHistoID(-1,kHChiMatch,pass,phase));
	// 
	TH2* hchiSA = 0; 
	ttl = Form("Pass%d_Chi2ITSSA_%s",pass,pref);
	hchiSA = new TH2F(ttl.Data(),ttl.Data(),nptbins,0,ptMax, nchbins/2,0.,chiMax/2);
	hchiSA->SetDirectory(0);
	dest->AddAtAndExpand(hchiSA,GetHistoID(-1,kHChiITSSA,pass,phase));
	// 	
      }
    } // loop over tracking passes
  }// loop over tracking phases
  //
  return dest;
}

//__________________________________________________________________
Int_t AliITSUTrackerGlo::GetHistoID(Int_t lr, Int_t hid, Int_t pass, Int_t phase)
{
  // get id for the requested histo
  if (lr<0) lr=-1;
  lr++;
  return pass*kHistosPass + phase*kHistosPhase + lr*kMaxHID + hid;
  //
}

#endif
