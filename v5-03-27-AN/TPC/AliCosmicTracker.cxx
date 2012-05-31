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
//
//
//
//
//-------------------------------------------------------
//          Implementation of the Cosmic tracker
//
//   Origin:  Xianguo Lu lu@physi.uni-heidelberg.de  Xianguo.Lu@cern.ch
//
//=========================================================================================
// Motivation:
// 
// In the default reconstruction in the ALICE the cosmic tracks are found  as two independent particles.
//
// In general any of  subtracks can be used for the physics studies. In order to avoid the double counting,
// the track from the upper hemisphere can be used. 
//
// The momentum resolution is determined by the lever arm (1/L^2) and by the number of clusters
// used for the track fitting (1/sqrt(Ncl)). 
// Combining/refitting  the two segments together significantly better momentum resolution can be obtained.
//    sigma_{1/pt} ~ 8x10^-3  - defaul tracking (e.g only upper track)
//    sigma_{1/pt} ~ 8x10^-4  - combined tracking
//===========================================================================================
// 
// Interface/Implementation:
// The class AliCosmicTracker provides functionality to find and refit the cosmic tracks. As a starting point, the events reconstruccted using standard tracking are used. 
// Input:  AliESDEvent
// Output: TClonesArray of the AliESDCosmicTrack
//         The array is stored as a data member of the tracker. 
//
// The cosmic tracker can be called in the user analysis code (standard analisys train, using the AliAnalysisTask, 
//    see e.g. AliAnalysisTaskCosmicTracker.h). In oreder to make an analysis simpler and faster it is planned to use the tracker already in the standard reconstruction (To be done).
//   
//===========================================================================================
// Algorithm:
// 1. Reads an ESD event      -  SetESDEvent() function
// 2. Loop over single tracks -  Process() function 
//    cuts are applied for individual ESD tracks (see function ESDtrckCut()). Only ESD tracks with TPCrefit, 
//    no kink and with ESDfriends will be selected. 
//    User defined cuts (as a pointer to the static function) can be used. (Expert usage)  

// 3. Double loop over tracks -  Process() function
//      a.) if not pair ( see function IsPair() ) continue; 
//         To accept the pair the tracks should be close together in the parameter space (AliExternalTrackParam - fP[0]-[4], also cut on ESD-phi and -theta)
//         Absolute, and relative (pull) cut are used
//         The cuts can be modified beyond default values via SetCut***().
//                  
//      b.) Each pair is fit via AliTPCCosmicTrackfit::CombineESDtracks
//      c.) For each pair one AliESDCosmicTrack object is stored in the fTrackStack, which can be passed out via TClonesArray *arr = fCosmicTracker->GetTrackStack();
//
//
//===========================================================================================
// Algorithm numerical debugging:
// The AliCosmicTracker can be used in the different debug/verbose level (see fDebugLevel) 
// Several intermediate variables can be stored in the trees, printout, or draw.
// Given functionality (dumping of variables to the tree) was also used for the tuning of the pair
// selection criterias, and for validation of the fit functionality. 
//
//===========================================================================================
// Usage:
// AliCosmicTracker *fCosmicTracker = new AliCosmicTracker(debuglevel, tag);
// fCosmicTracker->SetESDEvent(fESDEvent); //fTrackStack will be automatically cleared/emptied
// Int_t npair = fCosmicTracker->Process(processtag, kprint); //processtag only relavant if (debuglevel & 4) to draw the tracks in png; number of cosmic candidates are returned; if kprint the event is draw to png
//
//
// Advanced usage:
// fUserCut can be assigned externally so that additional ESDtrack cut can be applied in the very beginning together with those in ESDtrackCut()
//
// Example:
/*
//define static (important!!)  cut function in analysis task, e.g. AliAnalysisTaskCosmicTracker
//1) in AliAnalysisTaskCosmicTracker.h
static Bool_t TrackCut(AliESDtrack *trk);

//2) in AliAnalysisTaskCosmicTracker.cxx
Bool_t AliAnalysisTaskCosmicTracker::TrackCut(AliESDtrack *trk)
{
  //
  //external track cut in addition to the one in AliCosmicTracker (example)
  //
  if(!trk->GetTRDncls())
    return kFALSE;

  return kTRUE;
}
//set user cut function
fCosmicTracker = new AliCosmicTracker;
fCosmicTracker->SetUserESDtrackCut(AliAnalysisTaskCosmicTracker::TrackCut);
*/ 

#include <TTreeStream.h>

#include "AliESDEvent.h"
#include "AliTPCseed.h"
#include "AliTrackerBase.h"

#include "AliESDCosmicTrack.h"
#include "AliCosmicTracker.h"
#include "AliTPCCosmicUtils.h"
#include "AliTPCCosmicTrackfit.h"

AliCosmicTracker::AliCosmicTracker(const Int_t dlev, const TString tag): 
  fUserCut(0x0)
  , fStreamer(0x0), fDebugLevel(dlev)
  , fESDEvent(0x0)
  , fCosmicTrackfit(0x0)
  , fTrackStack(0x0)
  , fTrack0()
  , fTrack1()
  , fRawVtx(-999,-999,-999)
  , fRawDCA(-999)
  , fdPhi(-999)
  , fCutdPhi(-999)
  , fdTheta(-999)
  , fCutdTheta(-999)
  , fErrFlagESDtrackCut(-999)
  , fErrFlagIsPair(-999)
  , fErrFlagCosmicTrackfit(-999)
{
  //
  //constructor
  //

  if(fDebugLevel & 1)
    fStreamer = new TTreeSRedirector(Form("CosmicTracker_%s.root", tag.Data()));

  fCosmicTrackfit = new AliTPCCosmicTrackfit(0, "AliCosmicTracker");
  fTrackStack = new TClonesArray("AliESDCosmicTrack",100);
 
  for(Int_t ii=0; ii<5; ii++){
    fDelta[ii] = -999;
    fPull[ii] = -999;
  }

  fCutdPhi   = 19e-3*5;
  fCutdTheta = 10e-3*5;
  
  fCutPull[0] = 1.9 *10;
  fCutPull[1] = 1.5 *1e10;
  fCutPull[2] = 1.9 *10;//bug-fixed!
  fCutPull[3] = 0.4 *1e10;
  fCutPull[4] = 3.6 *10;

  fCutDelta[0] = 0.8   * 10;
  fCutDelta[1] = 2.7   * 10;
  fCutDelta[2] = 0.012 * 10;//bug-fixed!
  fCutDelta[3] = 0.007 * 10;
  fCutDelta[4] = 0.05  * 10;
}

AliCosmicTracker::~AliCosmicTracker()
{
  //
  //destructor
  //
  delete fStreamer;
  delete fCosmicTrackfit;
  delete fTrackStack;
}

void AliCosmicTracker::SetESDEvent(AliESDEvent *esd)
{
  //
  //set esd event
  //
  fESDEvent = esd;
  fTrackStack->Clear();
}

Int_t AliCosmicTracker::Process(const TString tag, const Bool_t kprint)
{
  //
  //double loop over combinations of esd tracks, cosmic event candidates sotred in fTrackStack
  //

  Int_t npair=0;
  const Int_t ntrk = fESDEvent->GetNumberOfTracks();
  Int_t trkcounter[ntrk];
  for(Int_t ii=0; ii<ntrk; ii++){
    trkcounter[ii]=0;
  }

  Double_t findabler0 = -999;
  Double_t findabler1 = -999;

  fErrFlagESDtrackCut = 0;
  fErrFlagIsPair = 0;
  fErrFlagCosmicTrackfit = 0;

  for(Int_t itrk=0; itrk<ntrk; itrk++){
    if(!ESDtrackCut(fESDEvent->GetTrack(itrk), findabler0)){
      continue;
    }

    for(Int_t jtrk=itrk+1; jtrk<ntrk; jtrk++){
      if(!ESDtrackCut(fESDEvent->GetTrack(jtrk), findabler1))
        continue;

      AliESDtrack * trk0 = fESDEvent->GetTrack(itrk);
      AliESDtrack * trk1 = fESDEvent->GetTrack(jtrk);
      if( IsPair(trk0, trk1) ){
        const Bool_t kfit = fCosmicTrackfit->CombineESDtracks(trk0, trk1);
        fErrFlagCosmicTrackfit = fCosmicTrackfit->GetStatus();

        if(kfit){
          fRawVtx = fCosmicTrackfit->ImpactParameter3D();
          fRawDCA = fCosmicTrackfit->ImpactParameter2D().Mag();

          const Int_t ncls              = fCosmicTrackfit->GetFitNcls();
          const Double_t leverarm       = fCosmicTrackfit->GetFitLeverArm();
          const Double_t chi2percluster = fCosmicTrackfit->GetChi2PerCluster();
          const Double_t impactD        = fCosmicTrackfit->GetImpactD();
          const Double_t impactZ        = fCosmicTrackfit->GetImpactZ();

          const Double_t findableratio  = TMath::Min(findabler0, findabler1);

          trkcounter[itrk]++;
          trkcounter[jtrk]++;
          const Bool_t isreuse          = (trkcounter[itrk]>1 || trkcounter[jtrk]>1);

          const TVector3 icU = fCosmicTrackfit->GetInnerClusterUp();
          const TVector3 icD = fCosmicTrackfit->GetInnerClusterLow();

          Int_t idup = itrk;
          Int_t idlow = jtrk;
          if(fCosmicTrackfit->IsSwap()){
            const Int_t idtmp = idup;
            idup = idlow;
            idlow = idtmp;

            AliExternalTrackParam tptmp = fTrack0;
            fTrack0 = fTrack1;
            fTrack1 = tptmp;
          }

          if(
             (fDebugLevel & 4) && 
             ( (isreuse && ntrk<=4) || kprint )
             ){
            AliESDtrack * trks[]={fESDEvent->GetTrack(idup), fESDEvent->GetTrack(idlow)};
            AliTPCCosmicUtils::DrawTracks(trks, Form("reuse_%03d_%03d_%03d_%s", ntrk, itrk, jtrk, tag.Data()));
          }
          if(
             (fDebugLevel & 8) &&
             impactD > 160 && findableratio < 0.56
             ){
            AliESDtrack * trks[]={fESDEvent->GetTrack(idup), fESDEvent->GetTrack(idlow)};
            AliTPCCosmicUtils::DrawTracks(trks, Form("largevtd_%.f_%.f_%03d_%03d_%03d_%s", impactD, findableratio, ntrk, itrk, jtrk, tag.Data()));
          }

          AliESDCosmicTrack costrk(idup, idlow, fCosmicTrackfit->GetTrackParamUp(), fCosmicTrackfit->GetTrackParamLow(), &fTrack0, &fTrack1, ncls, leverarm, chi2percluster, impactD, impactZ, isreuse, findableratio, icU, icD);
          new((*fTrackStack)[npair]) AliESDCosmicTrack(costrk);
          npair++;

          if(fDebugLevel & 1)
            WriteStreamer(ntrk, &costrk);

        }
        else{
          if(fDebugLevel & 16){
            if(ntrk==2){
              AliESDtrack * trks[]={trk0, trk1};
              AliTPCCosmicUtils::DrawTracks(trks, Form("failCosmicFit_%02d_%03d_%03d_%03d_%s", fCosmicTrackfit->GetStatus(), ntrk, itrk, jtrk, tag.Data()));
            }
          }
        }
      }
      else{
        if(fDebugLevel & 32){
          if(ntrk==2){
            AliESDtrack * trks[]={trk0, trk1};
            AliTPCCosmicUtils::DrawTracks(trks, Form("failIsPair_%02d_%03d_%03d_%03d_%s", fErrFlagIsPair, ntrk, itrk, jtrk, tag.Data()));
          }
        }
      }
    }
  }
  
  return npair;
}

Bool_t AliCosmicTracker::IsPair(AliESDtrack * trk0, AliESDtrack * trk1)
{
  //
  //check whether the two tracks come from one cosmic ray
  //

  //dphi + pi should = 0
  fdPhi   = AliTPCCosmicUtils::AngleInRange(trk0->Phi()   - trk1->Phi()   + TMath::Pi());
  if( TMath::Abs(fdPhi) > fCutdPhi ){
    fErrFlagIsPair = 1;
    return kFALSE;
  }

  fdTheta = AliTPCCosmicUtils::AngleInRange(trk0->Theta() + trk1->Theta() + TMath::Pi());
  if( TMath::Abs(fdTheta) > fCutdTheta ){
    fErrFlagIsPair = 2;
    return kFALSE;
  }

  //use fIp, the comments on the web is wrong (M. Ivanov)
  if(!trk0->GetInnerParam()){
    fErrFlagIsPair = 3;
    return kFALSE;
  }
  if(!trk1->GetInnerParam()){
    fErrFlagIsPair = 4;
    return kFALSE;
  }

  AliExternalTrackParam tmptrk[]={*(trk0->GetInnerParam()), *(trk1->GetInnerParam())};

  if(fDebugLevel & 2){
    printf("\n************************ raw ESD:\n");
    AliTPCCosmicUtils::PrintTrackParam(100, &(tmptrk[0]));
    AliTPCCosmicUtils::PrintTrackParam(101, &(tmptrk[1]));
    tmptrk[0].Print();
    tmptrk[1].Print();
  }

  Double_t xyz0[3], xyz1[3];
  tmptrk[0].GetXYZ(xyz0);
  tmptrk[1].GetXYZ(xyz1);
  const TVector3 gpos0(xyz0), gpos1(xyz1);

  //============================== rotate to common angle (M. Ivanov), since it is not possible to rotate from alpha1 to alpha0 via AliExternalTrackParam::Rotate
  const Double_t meanalpha = (gpos0-gpos1).Phi();
  const Double_t alpha0 = tmptrk[0].GetAlpha();

  //track0 closer to mean alpha
  if( TMath::Abs(AliTPCCosmicUtils::AngleInRange(meanalpha-alpha0)) <TMath::PiOver2() ){
    if( !AliTPCCosmicUtils::RotateSafe(&(tmptrk[0]), meanalpha) || 
        !AliTPCCosmicUtils::RotateSafe(&(tmptrk[1]), meanalpha+TMath::Pi()) ){
      fErrFlagIsPair = 5;
      return kFALSE;
    }
  }
  //track1 closer to mean alpha
  else{
    if( !AliTPCCosmicUtils::RotateSafe(&(tmptrk[1]), meanalpha) || 
        !AliTPCCosmicUtils::RotateSafe(&(tmptrk[0]), meanalpha+TMath::Pi()) ){
      fErrFlagIsPair = 6;
      return kFALSE;
    }
  }

  if(fDebugLevel & 2){
    printf("\n************************ after rotation!!\n");
    AliTPCCosmicUtils::PrintTrackParam(300, &(tmptrk[0]));
    AliTPCCosmicUtils::PrintTrackParam(301, &(tmptrk[1]));
    tmptrk[0].Print();
    tmptrk[1].Print();
  }

  //============================== propagate from TPC inner wall to x=0 with correct dedx
  const Double_t xTogo = 0.0;
  const Double_t maxStep = 1;
  const Bool_t rotateTo = kFALSE;
  const Double_t maxSnp = 0.8;
  Double_t eloss[2]={-1, 1};
  //default [0]=upper [1]=lower
  //tmptrk[0].phi<0 ==> [0]=lower
  if(AliTPCCosmicUtils::AngleInRange(tmptrk[0].Phi())<0){
    eloss[0]= 1;
    eloss[1]=-1;
  }

  for(Int_t ii=0; ii<2; ii++){
    if(!AliTrackerBase::PropagateTrackToBxByBz(&(tmptrk[ii]), xTogo, AliTPCCosmicUtils::Mass(), maxStep, rotateTo, maxSnp, (Int_t)(eloss[ii]))){
      fErrFlagIsPair = 7;
      return kFALSE;
    }
  }
  
  if(fDebugLevel & 2){
    printf("\n************************ after dedx corr:\n");
    AliTPCCosmicUtils::PrintTrackParam(200, &(tmptrk[0]));
    AliTPCCosmicUtils::PrintTrackParam(201, &(tmptrk[1]));
    tmptrk[0].Print();
    tmptrk[1].Print();
  }

  //____________________________________________________________________________________
  //____________________________________________________________________________________

  //ESD tracks after reconstruction all have x=0 and
  //TMath::Abs(alpha0 - alpha1)~ pi ==> back-to-back with angular resolution
  //[0]: local Y-coordinate of a track (cm);  9.945702e+01 -9.961257e+01 ==> opposite
  //[1]: local Z-coordinate of a track (cm); 2.677805e+01 2.711143e+01 ==> same
  //[2]: local sine of the track momentum azimuthal angle; should be the same!! bug-fixed
  //[3]: tangent of the track momentum dip angle; 1.563563e-01 -1.542005e-01 ==> opposite
  //[4]: 1/pt (1/(GeV/c)); -1.774935e-01 1.705480e-01 ==> opposite

  fDelta[0] = tmptrk[0].GetParameter()[0] + tmptrk[1].GetParameter()[0];
  fDelta[1] = tmptrk[0].GetParameter()[1] - tmptrk[1].GetParameter()[1];
  fDelta[2] = tmptrk[0].GetParameter()[2] - tmptrk[1].GetParameter()[2];//bug-fixed!! should use "-"
  fDelta[3] = tmptrk[0].GetParameter()[3] + tmptrk[1].GetParameter()[3];
  fDelta[4] = tmptrk[0].GetParameter()[4] + tmptrk[1].GetParameter()[4];

  fPull[0] = fDelta[0]/TMath::Sqrt(tmptrk[0].GetCovariance()[0] +tmptrk[1].GetCovariance()[0]);
  fPull[1] = fDelta[1]/TMath::Sqrt(tmptrk[0].GetCovariance()[2] +tmptrk[1].GetCovariance()[2]);
  fPull[2] = fDelta[2]/TMath::Sqrt(tmptrk[0].GetCovariance()[5] +tmptrk[1].GetCovariance()[5]);
  fPull[3] = fDelta[3]/TMath::Sqrt(tmptrk[0].GetCovariance()[9] +tmptrk[1].GetCovariance()[9]);
  fPull[4] = fDelta[4]/TMath::Sqrt(tmptrk[0].GetCovariance()[14]+tmptrk[1].GetCovariance()[14]);

  if(fDebugLevel & 2){
    for(Int_t ii=0; ii<5; ii++){
      printf("test %d %e %e -- %e\n", ii, tmptrk[0].GetParameter()[ii], tmptrk[1].GetParameter()[ii], fPull[ii]);
    }
  }

  for(Int_t ii=0; ii<5; ii++){
    if( TMath::Abs(fPull[ii])  > fCutPull[ii] ){
      fErrFlagIsPair = 10+ii;
      return kFALSE;
    }
    if( TMath::Abs(fDelta[ii]) > fCutDelta[ii] ){
      fErrFlagIsPair = 20+ii;
      return kFALSE;
    }
  }

  fTrack0 = tmptrk[0];
  fTrack1 = tmptrk[1];

  return kTRUE;
}

Bool_t AliCosmicTracker::ESDtrackCut(AliESDtrack * trk, Double_t &findabler) 
{
  //
  //cut on track quality (kink, TPCrefit, findable ratio) and require TPC seed
  //

  if(fUserCut){
    if(!fUserCut(trk)){
      fErrFlagESDtrackCut = 1;
      return kFALSE;
    }
  }

  //reject kink
  if(trk->GetKinkIndex(0)>0){
    fErrFlagESDtrackCut = 2;
    return kFALSE;
  }

  //require refit
  if(!trk->IsOn(AliESDtrack::kTPCrefit)){
    fErrFlagESDtrackCut = 3;
    return kFALSE;
  }

  // due to drift velocity calibration, a track crossing Z=0 may be reconstructed as 2 ESD tracks, so two pairs are formed, each with one part of this track. Solution: cut on findable ratio (require > 0.5) to remove split tracks due to drift velocity calibration systematics on different sides
  //there is some remaining with isreuse = true, user should cut on findable ratio according to the fraction of isreuse
  findabler = -999;
  if(!trk->GetTPCNclsF()){
    fErrFlagESDtrackCut = 4;
    return kFALSE;
  }

  findabler = (Double_t)trk->GetTPCNcls()/(Double_t) trk->GetTPCNclsF();

  if(findabler < CutFindable() ){
    fErrFlagESDtrackCut = 5;
    return kFALSE;
  }

  //cut on # TPC ncls on each ESDtrack
  if(trk->GetTPCncls()<AliTPCCosmicUtils::NclsMin()){
    fErrFlagESDtrackCut = 6;
    return kFALSE;
  }

  //require ESDfriends
  if(!AliTPCCosmicUtils::GetTPCseed(trk)){
    fErrFlagESDtrackCut = 7;
    return kFALSE;
  }

  return kTRUE;
}

Int_t AliCosmicTracker::GetErrFlag() const
{
  //
  //return the error status in process
  //
  return fErrFlagESDtrackCut + fErrFlagIsPair*100 + fErrFlagCosmicTrackfit*10000;
}

void AliCosmicTracker::WriteStreamer(Int_t ntrk, AliESDCosmicTrack *costrk)
{
  //
  //output to streamer
  //

  (*fStreamer)<<"CosmicTracker_Streamer"<<
    "ntrk="<<ntrk<<

    "costrk="<<costrk<<

    "rawvtx="<<&fRawVtx<<
    "rawdca="<<fRawDCA<<

    "dphi="<<fdPhi<<
    "dtheta="<<fdTheta<<
    "pull0="<<fPull[0]<<
    "pull1="<<fPull[1]<<
    "pull2="<<fPull[2]<<
    "pull3="<<fPull[3]<<
    "pull4="<<fPull[4]<<
    "delta0="<<fDelta[0]<<
    "delta1="<<fDelta[1]<<
    "delta2="<<fDelta[2]<<
    "delta3="<<fDelta[3]<<
    "delta4="<<fDelta[4]<<
    "\n";
}

TClonesArray *AliCosmicTracker::FindCosmic(AliESDEvent *event, const Bool_t kadd)
{
  //
  //do cosmic combined trackfit
  //

  AliCosmicTracker cosmicTracker;
  cosmicTracker.SetESDEvent(event);
  const Int_t npair = cosmicTracker.Process();
  const TClonesArray *arr = cosmicTracker.GetTrackStack();

  TClonesArray *stackCosmic = 0x0;
  if(kadd){
    for(Int_t ip=0; ip<npair; ip++){
      const AliESDCosmicTrack * esdcos = (AliESDCosmicTrack*) arr->At(ip);
      event->AddCosmicTrack(esdcos);
    }
    printf("AliCosmicTracker::FindCosmic: event %d: Number of cosmic pairs by AliCosmicTracker %d out of %d tracks, err %d\n", event->GetEventNumberInFile(), npair, event->GetNumberOfTracks(), cosmicTracker.GetErrFlag());
  }
  else{
    stackCosmic = new TClonesArray(*arr);
  }

  return stackCosmic;
}

