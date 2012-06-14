/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                Dielectron Analysis Main class                         //
//                                                                       //
/*
Framework to perform event selectoin, single track selection and track pair
selection.

Convention for the signs of the pair in fPairCandidates:
The names are available via the function PairClassName(Int_t i)

0: ev1+ ev1+  (same event like sign +)
1: ev1+ ev1-  (same event unlike sign)
2: ev1- ev1-  (same event like sign -)

3: ev1+ ev2+  (mixed event like sign +)
4: ev1- ev2+  (mixed event unlike sign -+)
6: ev1+ ev2-  (mixed event unlike sign +-)
7: ev1- ev2-  (mixed event like sign -)

5: ev2+ ev2+  (same event like sign +)
8: ev2+ ev2-  (same event unlike sign)
9: ev2- ev2-  (same event like sign -)

10: ev1+ ev1- (same event track rotation)

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TList.h>
#include <TMath.h>
#include <TObject.h>

#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliKFParticle.h>

#include <AliEventplane.h>
#include <AliVEvent.h>
#include <AliVParticle.h>
#include <AliVTrack.h>
#include "AliDielectronPair.h"
#include "AliDielectronHistos.h"
#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronTrackRotator.h"
#include "AliDielectronDebugTree.h"
#include "AliDielectronSignalMC.h"
#include "AliDielectronMixingHandler.h"

#include "AliDielectron.h"

ClassImp(AliDielectron)

const char* AliDielectron::fgkTrackClassNames[4] = {
  "ev1+",
  "ev1-",
  "ev2+",
  "ev2-"
};

const char* AliDielectron::fgkPairClassNames[11] = {
  "ev1+_ev1+",
  "ev1+_ev1-",
  "ev1-_ev1-",
  "ev1+_ev2+",
  "ev1-_ev2+",
  "ev2+_ev2+",
  "ev1+_ev2-",
  "ev1-_ev2-",
  "ev2+_ev2-",
  "ev2-_ev2-",
  "ev1+_ev1-_TR"
};

//________________________________________________________________
AliDielectron::AliDielectron() :
  TNamed("AliDielectron","AliDielectron"),
  fEventFilter("EventFilter"),
  fTrackFilter("TrackFilter"),
  fPairPreFilter("PairPreFilter"),
  fPairPreFilterLegs("PairPreFilterLegs"),
  fPairFilter("PairFilter"),
  fEventPlanePreFilter("EventPlanePreFilter"),
  fEventPlanePOIPreFilter("EventPlanePOIPreFilter"),
  fPdgMother(443),
  fPdgLeg1(11),
  fPdgLeg2(11),
  fSignalsMC(0x0),
  fNoPairing(kFALSE),
  fHistos(0x0),
  fPairCandidates(new TObjArray(11)),
  fCfManagerPair(0x0),
  fTrackRotator(0x0),
  fDebugTree(0x0),
  fMixing(0x0),
  fPreFilterEventPlane(kFALSE),
  fLikeSignSubEvents(kFALSE),
  fPreFilterUnlikeOnly(kFALSE),
  fPreFilterAllSigns(kFALSE),
  fHasMC(kFALSE),
  fStoreRotatedPairs(kFALSE),
  fDontClearArrays(kFALSE),
  fEstimatorFilename(""),
  fTRDpidCorrectionFilename(""),
  fVZEROCalibrationFilename(""),
  fVZERORecenteringFilename("")
{
  //
  // Default constructor
  //

}

//________________________________________________________________
AliDielectron::AliDielectron(const char* name, const char* title) :
  TNamed(name,title),
  fEventFilter("EventFilter"),
  fTrackFilter("TrackFilter"),
  fPairPreFilter("PairPreFilter"),
  fPairPreFilterLegs("PairPreFilterLegs"),
  fPairFilter("PairFilter"),
  fEventPlanePreFilter("EventPlanePreFilter"),
  fEventPlanePOIPreFilter("EventPlanePOIPreFilter"),
  fPdgMother(443),
  fPdgLeg1(11),
  fPdgLeg2(11),
  fSignalsMC(0x0),
  fNoPairing(kFALSE),
  fHistos(0x0),
  fPairCandidates(new TObjArray(11)),
  fCfManagerPair(0x0),
  fTrackRotator(0x0),
  fDebugTree(0x0),
  fMixing(0x0),
  fPreFilterEventPlane(kFALSE),
  fLikeSignSubEvents(kFALSE),
  fPreFilterUnlikeOnly(kFALSE),
  fPreFilterAllSigns(kFALSE),
  fHasMC(kFALSE),
  fStoreRotatedPairs(kFALSE),
  fDontClearArrays(kFALSE),
  fEstimatorFilename(""),
  fTRDpidCorrectionFilename(""),
  fVZEROCalibrationFilename(""),
  fVZERORecenteringFilename("")
{
  //
  // Named constructor
  //
  
}

//________________________________________________________________
AliDielectron::~AliDielectron()
{
  //
  // Default destructor
  //
  if (fHistos) delete fHistos;
  if (fPairCandidates) delete fPairCandidates;
  if (fDebugTree) delete fDebugTree;
  if (fMixing) delete fMixing;
  if (fSignalsMC) delete fSignalsMC;
  if (fCfManagerPair) delete fCfManagerPair;
}

//________________________________________________________________
void AliDielectron::Init()
{
  //
  // Initialise objects
  //

  if(GetHasMC()) AliDielectronMC::Instance()->SetHasMC(GetHasMC());
   
  InitPairCandidateArrays();
   
  if (fCfManagerPair) {
    fCfManagerPair->SetSignalsMC(fSignalsMC);
    fCfManagerPair->InitialiseContainer(fPairFilter);
  }
  if (fTrackRotator)  {
    fTrackRotator->SetTrackArrays(&fTracks[0],&fTracks[1]);
    fTrackRotator->SetPdgLegs(fPdgLeg1,fPdgLeg2);
  }
  if (fDebugTree) fDebugTree->SetDielectron(this);
  if(fEstimatorFilename.Contains(".root")) AliDielectronVarManager::InitEstimatorAvg(fEstimatorFilename.Data());
  if(fTRDpidCorrectionFilename.Contains(".root")) AliDielectronVarManager::InitTRDpidEffHistograms(fTRDpidCorrectionFilename.Data());
  if(fVZEROCalibrationFilename.Contains(".root")) AliDielectronVarManager::SetVZEROCalibrationFile(fVZEROCalibrationFilename.Data());
  if(fVZERORecenteringFilename.Contains(".root")) AliDielectronVarManager::SetVZERORecenteringFile(fVZERORecenteringFilename.Data());
  
  if (fMixing) fMixing->Init(this);
} 

//________________________________________________________________
void AliDielectron::Process(AliVEvent *ev1, AliVEvent *ev2)
{
  //
  // Process the events
  //

  //at least first event is needed!
  if (!ev1){
    AliError("At least first event must be set!");
    return;
  }
  AliDielectronVarManager::SetEvent(ev1);
  if (fMixing){
    //set mixing bin to event data
    Int_t bin=fMixing->FindBin(AliDielectronVarManager::GetData());
    AliDielectronVarManager::SetValue(AliDielectronVarManager::kMixingBin,bin);
  }

  //in case we have MC load the MC event and process the MC particles
  if (AliDielectronMC::Instance()->HasMC()) {
    if (!AliDielectronMC::Instance()->ConnectMCEvent()){
      AliError("Could not properly connect the MC event, skipping this event!");
      return;
    }
    ProcessMC(ev1);
  }
  
  //if candidate array doesn't exist, create it
  if (!fPairCandidates->UncheckedAt(0)) {
    InitPairCandidateArrays();
  } else {
    ClearArrays();
  }

  //mask used to require that all cuts are fulfilled
  UInt_t selectedMask=(1<<fEventFilter.GetCuts()->GetEntries())-1;

  //apply event cuts
    if ((ev1&&fEventFilter.IsSelected(ev1)!=selectedMask) ||
        (ev2&&fEventFilter.IsSelected(ev2)!=selectedMask)) return;
  
//   AliDielectronVarManager::SetEvent(ev1); // why a second time???

  //fill track arrays for the first event
  if (ev1){
    FillTrackArrays(ev1);
    if (((fPreFilterAllSigns)||(fPreFilterUnlikeOnly)) && ( fPairPreFilter.GetCuts()->GetEntries()>0 )) PairPreFilter(0, 1, fTracks[0], fTracks[1]);
  }


  //fill track arrays for the second event
  if (ev2) {
    FillTrackArrays(ev2,1);
    if (((fPreFilterAllSigns)||(fPreFilterUnlikeOnly)) && ( fPairPreFilter.GetCuts()->GetEntries()>0 )) PairPreFilter(2, 3, fTracks[2], fTracks[3]);
  }

  // TPC event plane correction
  AliEventplane *cevplane = new AliEventplane();
  if (ev1 && fPreFilterEventPlane && ( fEventPlanePreFilter.GetCuts()->GetEntries()>0 || fEventPlanePOIPreFilter.GetCuts()->GetEntries()>0)) 
    EventPlanePreFilter(0, 1, fTracks[0], fTracks[1], ev1, cevplane);
  
  if (!fNoPairing){
    // create pairs and fill pair candidate arrays
    for (Int_t itrackArr1=0; itrackArr1<4; ++itrackArr1){
      for (Int_t itrackArr2=itrackArr1; itrackArr2<4; ++itrackArr2){
        FillPairArrays(itrackArr1, itrackArr2);
      }
    }

    //track rotation
    if (fTrackRotator) {
      fTrackRotator->SetEvent(ev1);
      FillPairArrayTR();
    }
  }

  //fill debug tree if a manager is attached
  if (fDebugTree) FillDebugTree();

  //process event mixing
  if (fMixing) {
    fMixing->Fill(ev1,this);
//     FillHistograms(0x0,kTRUE);
  }

  //in case there is a histogram manager, fill the QA histograms
  if (fHistos) FillHistograms(ev1);

  // clear arrays
  if (!fDontClearArrays) ClearArrays();
  AliDielectronVarManager::SetTPCEventPlane(0x0);
  delete cevplane;
}

//________________________________________________________________
void AliDielectron::ProcessMC(AliVEvent *ev1)
{
  //
  // Process the MC data
  //

  AliDielectronMC *dieMC=AliDielectronMC::Instance();

  if (fHistos) FillHistogramsMC(dieMC->GetMCEvent(), ev1);

  if(!fSignalsMC) return;
  //loop over all MC data and Fill the CF container if it exist
  if (!fCfManagerPair) return;
  fCfManagerPair->SetPdgMother(fPdgMother);
  if(!fCfManagerPair->GetStepForMCtruth()) return;

  // signals to be studied
  Int_t nSignals = fSignalsMC->GetEntries();

  // initialize 2D arrays of labels for particles from each MC signal
  Int_t** labels1;      // labels for particles satisfying branch 1
  Int_t** labels2;      // labels for particles satisfying branch 2
  Int_t** labels12;     // labels for particles satisfying both branches
  labels1 = new Int_t*[nSignals];
  labels2 = new Int_t*[nSignals];
  labels12 = new Int_t*[nSignals];
  Int_t* indexes1=new Int_t[nSignals];
  Int_t* indexes2=new Int_t[nSignals];
  Int_t* indexes12=new Int_t[nSignals];
  for(Int_t isig=0;isig<nSignals;++isig) {
    *(labels1+isig) = new Int_t[dieMC->GetNMCTracks()];
    *(labels2+isig) = new Int_t[dieMC->GetNMCTracks()];
    *(labels12+isig) = new Int_t[dieMC->GetNMCTracks()];
    for(Int_t ip=0; ip<dieMC->GetNMCTracks();++ip) {
      labels1[isig][ip] = -1;
      labels2[isig][ip] = -1;
      labels12[isig][ip] = -1;
    }
    indexes1[isig]=0;
    indexes2[isig]=0;
    indexes12[isig]=0;
  }

  Bool_t truth1=kFALSE;
  Bool_t truth2=kFALSE;
  // loop over the MC tracks
  for(Int_t ipart=0; ipart<dieMC->GetNMCTracks(); ++ipart) {
    for(Int_t isig=0; isig<nSignals; ++isig) {       // loop over signals
      // Proceed only if this signal is required in the pure MC step
      // NOTE: Some signals can be satisfied by many particles and this leads to high
      //       computation times (e.g. secondary electrons from the GEANT transport). Be aware of this!!
      if(!((AliDielectronSignalMC*)fSignalsMC->At(isig))->GetFillPureMCStep()) continue;

      truth1 = dieMC->IsMCTruth(ipart, (AliDielectronSignalMC*)fSignalsMC->At(isig), 1);
      truth2 = dieMC->IsMCTruth(ipart, (AliDielectronSignalMC*)fSignalsMC->At(isig), 2);

      // particles satisfying both branches are treated separately to avoid double counting during pairing
      if(truth1 && truth2) {
	labels12[isig][indexes12[isig]] = ipart;
	++indexes12[isig];
      }
      else {
	if(truth1) {
	  labels1[isig][indexes1[isig]] = ipart;
	  ++indexes1[isig];
	}
	if(truth2) {
	  labels2[isig][indexes2[isig]] = ipart;
	  ++indexes2[isig];
	}
      }
    }
  }  // end loop over MC particles

  // Do the pairing and fill the CF container with pure MC info
  for(Int_t isig=0; isig<nSignals; ++isig) {
    // mix the particles which satisfy only one of the signal branches
    for(Int_t i1=0;i1<indexes1[isig];++i1) {
      for(Int_t i2=0;i2<indexes2[isig];++i2) {
	fCfManagerPair->FillMC(labels1[isig][i1], labels2[isig][i2], isig);
      }
    }
    // mix the particles which satisfy both branches
    for(Int_t i1=0;i1<indexes12[isig];++i1) {
      for(Int_t i2=0; i2<i1; ++i2) {
	fCfManagerPair->FillMC(labels12[isig][i1], labels12[isig][i2], isig);
      }
    }
  }    // end loop over signals

  // release the memory
  for(Int_t isig=0;isig<nSignals;++isig) {
    delete [] *(labels1+isig);
    delete [] *(labels2+isig);
    delete [] *(labels12+isig);
  }
  delete [] labels1;
  delete [] labels2;
  delete [] labels12;
  delete [] indexes1;
  delete [] indexes2;
  delete [] indexes12;
}

//________________________________________________________________
void AliDielectron::FillHistogramsTracks(TObjArray **tracks)
{
  //
  // Fill Histogram information for tracks after prefilter
  // ignore mixed events - for prefilter, only single tracks +/- are relevant 
  //
  
  TString  className,className2;
  Double_t values[AliDielectronVarManager::kNMaxValues];
  
  //Fill track information, separately for the track array candidates
  for (Int_t i=0; i<2; ++i){
    className.Form("Pre_%s",fgkTrackClassNames[i]);
    if (!fHistos->GetHistogramList()->FindObject(className.Data())) continue;
    Int_t ntracks=tracks[i]->GetEntriesFast();
    for (Int_t itrack=0; itrack<ntracks; ++itrack){
      AliDielectronVarManager::Fill(tracks[i]->UncheckedAt(itrack), values);
      fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
    }
  }
}


//________________________________________________________________
void AliDielectron::FillHistogramsMC(const AliMCEvent *ev, AliVEvent *ev1)
{
  //
  // Fill Histogram information for MCEvents
  //

  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  // Fill event information
  AliDielectronVarManager::Fill(ev1, values);    // ESD/AOD information
  AliDielectronVarManager::Fill(ev, values);     // MC truth info
  if (fHistos->GetHistogramList()->FindObject("MCEvent"))
    fHistos->FillClass("MCEvent", AliDielectronVarManager::kNMaxValues, values);
}


//________________________________________________________________
void AliDielectron::FillHistograms(const AliVEvent *ev, Bool_t pairInfoOnly)
{
  //
  // Fill Histogram information for tracks and pairs
  //
  
  TString  className,className2;
  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  //Fill event information
  if (ev){  //TODO: Why not use GetData() ??? See below event plane stuff!!!
    AliDielectronVarManager::Fill(ev, values); //data should already be stored in AliDielectronVarManager from SetEvent, does EV plane correction rely on this???
      if (fMixing){
    //set mixing bin to event data
    Int_t bin=fMixing->FindBin(values);
    values[AliDielectronVarManager::kMixingBin]=bin;
  }

    if (fHistos->GetHistogramList()->FindObject("Event"))
//       fHistos->FillClass("Event", AliDielectronVarManager::kNMaxValues, AliDielectronVarManager::GetData());
      fHistos->FillClass("Event", AliDielectronVarManager::kNMaxValues, values); //check event plane stuff and replace with above...
  }
  
  //Fill track information, separately for the track array candidates
  if (!pairInfoOnly){
    for (Int_t i=0; i<4; ++i){
      className.Form("Track_%s",fgkTrackClassNames[i]);
      if (!fHistos->GetHistogramList()->FindObject(className.Data())) continue;
      Int_t ntracks=fTracks[i].GetEntriesFast();
      for (Int_t itrack=0; itrack<ntracks; ++itrack){
        AliDielectronVarManager::Fill(fTracks[i].UncheckedAt(itrack), values);
        fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
      }
    }
  }

  //Fill Pair information, separately for all pair candidate arrays and the legs
  TObjArray arrLegs(100);
  for (Int_t i=0; i<10; ++i){
    className.Form("Pair_%s",fgkPairClassNames[i]);
    className2.Form("Track_Legs_%s",fgkPairClassNames[i]);
    Bool_t pairClass=fHistos->GetHistogramList()->FindObject(className.Data())!=0x0;
    Bool_t legClass=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;
    if (!pairClass&&!legClass) continue;
    Int_t ntracks=PairArray(i)->GetEntriesFast();
    for (Int_t ipair=0; ipair<ntracks; ++ipair){
      AliDielectronPair *pair=static_cast<AliDielectronPair*>(PairArray(i)->UncheckedAt(ipair));
      
      //fill pair information
      if (pairClass){
        AliDielectronVarManager::Fill(pair, values);
        fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
      }

      //fill leg information, don't fill the information twice
      if (legClass){
        AliVParticle *d1=pair->GetFirstDaughter();
        AliVParticle *d2=pair->GetSecondDaughter();
        if (!arrLegs.FindObject(d1)){
          AliDielectronVarManager::Fill(d1, values);
          fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
          arrLegs.Add(d1);
        }
        if (!arrLegs.FindObject(d2)){
          AliDielectronVarManager::Fill(d2, values);
          fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
          arrLegs.Add(d2);
        }
      }
    }
    if (legClass) arrLegs.Clear();
  }
  
}
//________________________________________________________________
void AliDielectron::FillHistogramsPair(AliDielectronPair *pair,Bool_t fromPreFilter/*=kFALSE*/)
{
  //
  // Fill Histogram information for pairs and the track in the pair
  // NOTE: in this funtion the leg information may be filled multiple
  //       times. This funtion is used in the track rotation pairing
  //       and those legs are not saved!
  //
  TString  className,className2;
  Double_t values[AliDielectronVarManager::kNMaxValues];
  
  //Fill Pair information, separately for all pair candidate arrays and the legs
  TObjArray arrLegs(100);
  const Int_t type=pair->GetType();
  if (fromPreFilter) {
    className.Form("RejPair_%s",fgkPairClassNames[type]);
    className2.Form("RejTrack_%s",fgkPairClassNames[type]);
  } else {
    className.Form("Pair_%s",fgkPairClassNames[type]);
    className2.Form("Track_Legs_%s",fgkPairClassNames[type]);
  }
  
  Bool_t pairClass=fHistos->GetHistogramList()->FindObject(className.Data())!=0x0;
  Bool_t legClass=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;
  
  //fill pair information
  if (pairClass){
    AliDielectronVarManager::Fill(pair, values);
    fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
  }

  if (legClass){
    AliVParticle *d1=pair->GetFirstDaughter();
    AliDielectronVarManager::Fill(d1, values);
    fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
    
    AliVParticle *d2=pair->GetSecondDaughter();
    AliDielectronVarManager::Fill(d2, values);
    fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
  }
}

//________________________________________________________________
void AliDielectron::FillTrackArrays(AliVEvent * const ev, Int_t eventNr)
{
  //
  // select tracks and fill track candidate arrays
  // eventNr = 0: First  event, use track arrays 0 and 1
  // eventNr = 1: Second event, use track arrays 2 and 3
  //
  
  Int_t ntracks=ev->GetNumberOfTracks();
  UInt_t selectedMask=(1<<fTrackFilter.GetCuts()->GetEntries())-1;
  for (Int_t itrack=0; itrack<ntracks; ++itrack){
    //get particle
    AliVParticle *particle=ev->GetTrack(itrack);
    //TODO: temporary solution, perhaps think about a better implementation
    //      This is needed to use AliESDpidCuts, which relies on the ESD event
    //      is set as a AliESDtrack attribute... somehow ugly!
    if (ev->IsA()==AliESDEvent::Class()){
      AliESDtrack *track=static_cast<AliESDtrack*>(particle);
      track->SetESDEvent(static_cast<AliESDEvent*>(ev)); //only in trunk...
    }
    //apply track cuts
    if (fTrackFilter.IsSelected(particle)!=selectedMask) continue;

    //fill selected particle into the corresponding track arrays
    Short_t charge=particle->Charge();
    if (charge>0)      fTracks[eventNr*2].Add(particle);
    else if (charge<0) fTracks[eventNr*2+1].Add(particle);
  }
}

//________________________________________________________________
void AliDielectron::EventPlanePreFilter(Int_t arr1, Int_t arr2, TObjArray arrTracks1, TObjArray arrTracks2, const AliVEvent *ev, AliEventplane *cevplane)
{
  //
  // Prefilter tracks and tracks from pairs
  // Needed for rejection in the Q-Vector of the event plane
  // remove contribution of all tracks to the Q-vector that are in invariant mass window 
  //
  AliEventplane *evplane = const_cast<AliVEvent *>(ev)->GetEventplane();
//   AliEventplane *evplane = ev->GetEventplane();
  if(!evplane) return;
  
  // do not change these vectors qref
  TVector2 * const qref  = evplane->GetQVector();
  if(!qref) return;
  // random subevents
  TVector2 *qrsub1 = evplane->GetQsub1();
  TVector2 *qrsub2 = evplane->GetQsub2();

  // copy references
  TVector2 *qcorr  = new TVector2(*qref);
  TVector2 *qcsub1 = 0x0;
  TVector2 *qcsub2 = 0x0;
  //  printf("qrsub1 %p %f \n",qrsub1,qrsub1->X());


  // eta gap ?
  Bool_t etagap = kFALSE;
  for (Int_t iCut=0; iCut<fEventPlanePreFilter.GetCuts()->GetEntries();++iCut) {
    TString cutName=fEventPlanePreFilter.GetCuts()->At(iCut)->GetName();
    if(cutName.Contains("eta") || cutName.Contains("Eta"))  etagap=kTRUE;
  }

  // subevent separation
  if(fLikeSignSubEvents || etagap) {
    qcsub1 = new TVector2(*qcorr);
    qcsub2 = new TVector2(*qcorr);

    Int_t ntracks=ev->GetNumberOfTracks();
    
    // track removals
    for (Int_t itrack=0; itrack<ntracks; ++itrack){
      AliVParticle *particle=ev->GetTrack(itrack);
      AliVTrack *track= static_cast<AliVTrack*>(particle);
      if (!track) continue;

      Double_t cQX     = evplane->GetQContributionX(track);
      Double_t cQY     = evplane->GetQContributionY(track);
      
      // by charge sub1+ sub2-
      if(fLikeSignSubEvents) {
	Short_t charge=track->Charge();
	if (charge<0) qcsub1->Set(qcsub1->X()-cQX, qcsub1->Y()-cQY);
	if (charge>0) qcsub2->Set(qcsub2->X()-cQX, qcsub2->Y()-cQY);
      }
      // by eta sub1+ sub2-
      if(etagap) {
	Double_t eta=track->Eta();
	if (eta<0.0) qcsub1->Set(qcsub1->X()-cQX, qcsub1->Y()-cQY);
	if (eta>0.0) qcsub2->Set(qcsub2->X()-cQX, qcsub2->Y()-cQY);
      }
    }
  }
  else {
    // by a random
    qcsub1 = new TVector2(*qrsub1);
    qcsub2 = new TVector2(*qrsub2);
  }
  
  // apply cuts, e.g. etagap 
  if(fEventPlanePreFilter.GetCuts()->GetEntries()) {
    UInt_t selectedMask=(1<<fEventPlanePreFilter.GetCuts()->GetEntries())-1;
    Int_t ntracks=ev->GetNumberOfTracks();
    for (Int_t itrack=0; itrack<ntracks; ++itrack){
      AliVParticle *particle=ev->GetTrack(itrack);
      AliVTrack *track= static_cast<AliVTrack*>(particle);
      if (!track) continue;
      
      //event plane cuts
      UInt_t cutMask=fEventPlanePreFilter.IsSelected(track);
      //apply cut
      if (cutMask==selectedMask) continue;

      Double_t cQX     = 0.0;
      Double_t cQY     = 0.0;
      if(!etagap) {
	cQX = evplane->GetQContributionX(track);
	cQY = evplane->GetQContributionY(track);
      }
      Double_t cQXsub1 = evplane->GetQContributionXsub1(track);
      Double_t cQYsub1 = evplane->GetQContributionYsub1(track);
      Double_t cQXsub2 = evplane->GetQContributionXsub2(track);
      Double_t cQYsub2 = evplane->GetQContributionYsub2(track);      

      // update Q vectors
      qcorr->Set(qcorr->X()-cQX, qcorr->Y()-cQY);
      qcsub1->Set(qcsub1->X()-cQXsub1, qcsub1->Y()-cQYsub1);
      qcsub2->Set(qcsub2->X()-cQXsub2, qcsub2->Y()-cQYsub2);
    }
  }

  // POI (particle of interest) rejection
  Int_t pairIndex=GetPairIndex(arr1,arr2);
  
  Int_t ntrack1=arrTracks1.GetEntriesFast();
  Int_t ntrack2=arrTracks2.GetEntriesFast();
  AliDielectronPair candidate;
  
  UInt_t selectedMask=(1<<fEventPlanePOIPreFilter.GetCuts()->GetEntries())-1;
  for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1){
    Int_t end=ntrack2;
    if (arr1==arr2) end=itrack1;
    Bool_t accepted=kFALSE;
    for (Int_t itrack2=0; itrack2<end; ++itrack2){
      TObject *track1=arrTracks1.UncheckedAt(itrack1);
      TObject *track2=arrTracks2.UncheckedAt(itrack2);
      if (!track1 || !track2) continue;
      //create the pair
      candidate.SetTracks(static_cast<AliVTrack*>(track1), fPdgLeg1,
                          static_cast<AliVTrack*>(track2), fPdgLeg2);
      
      candidate.SetType(pairIndex);
      candidate.SetLabel(AliDielectronMC::Instance()->GetLabelMotherWithPdg(&candidate,fPdgMother));
      
      //event plane cuts
      UInt_t cutMask=fEventPlanePOIPreFilter.IsSelected(&candidate);
      //apply cut
      if (cutMask==selectedMask) continue;

      accepted=kTRUE;
      //remove the tracks from the Track arrays
      arrTracks2.AddAt(0x0,itrack2);
    }
      if ( accepted ) arrTracks1.AddAt(0x0,itrack1);
  }
  //compress the track arrays
  arrTracks1.Compress();
  arrTracks2.Compress();
  
  
  //Modify the components: subtract the tracks
  ntrack1=arrTracks1.GetEntriesFast();
  ntrack2=arrTracks2.GetEntriesFast();
  
  // remove leg1 contribution
  for (Int_t itrack=0; itrack<ntrack1; ++itrack){
    AliVTrack *track= static_cast<AliVTrack*>(arrTracks1.UncheckedAt(itrack));
    if (!track) continue;
    
    Double_t cQX     = evplane->GetQContributionX(track);
    Double_t cQY     = evplane->GetQContributionY(track);
    Double_t cQXsub1 = evplane->GetQContributionXsub1(track);
    Double_t cQYsub1 = evplane->GetQContributionYsub1(track);
    Double_t cQXsub2 = evplane->GetQContributionXsub2(track);
    Double_t cQYsub2 = evplane->GetQContributionYsub2(track);
    
    // update Q vectors
    qcorr->Set(qcorr->X()-cQX, qcorr->Y()-cQY);
    qcsub1->Set(qcsub1->X()-cQXsub1, qcsub1->Y()-cQYsub1);
    qcsub2->Set(qcsub2->X()-cQXsub2, qcsub2->Y()-cQYsub2);
  }
  // remove leg2 contribution
  for (Int_t itrack=0; itrack<ntrack2; ++itrack){
    AliVTrack *track= static_cast<AliVTrack*>(arrTracks2.UncheckedAt(itrack));
    if (!track) continue;
    
    Double_t cQX     = evplane->GetQContributionX(track);
    Double_t cQY     = evplane->GetQContributionY(track);
    Double_t cQXsub1 = evplane->GetQContributionXsub1(track);
    Double_t cQYsub1 = evplane->GetQContributionYsub1(track);
    Double_t cQXsub2 = evplane->GetQContributionXsub2(track);
    Double_t cQYsub2 = evplane->GetQContributionYsub2(track);
    
    // update Q vectors
    qcorr->Set(qcorr->X()-cQX, qcorr->Y()-cQY);
    qcsub1->Set(qcsub1->X()-cQXsub1, qcsub1->Y()-cQYsub1);
    qcsub2->Set(qcsub2->X()-cQXsub2, qcsub2->Y()-cQYsub2);
  }

  //  printf("qrsub1 %p %f \t qcsub1 %p %f \n",qrsub1,qrsub1->X(),qcsub1,qcsub1->X());
  // set AliEventplane with corrected values
  cevplane->SetQVector(qcorr);
  cevplane->SetQsub(qcsub1, qcsub2);
  AliDielectronVarManager::SetTPCEventPlane(cevplane);
}

//________________________________________________________________
void AliDielectron::PairPreFilter(Int_t arr1, Int_t arr2, TObjArray &arrTracks1, TObjArray &arrTracks2)
{
  //
  // Prefilter tracks from pairs
  // Needed for datlitz rejections
  // remove all tracks from the Single track arrays that pass the cuts in this filter
  //

  Int_t ntrack1=arrTracks1.GetEntriesFast();
  Int_t ntrack2=arrTracks2.GetEntriesFast();
  AliDielectronPair candidate;

  // flag arrays for track removal
  Bool_t *bTracks1 = new Bool_t[ntrack1];
  for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1) bTracks1[itrack1]=kFALSE;
  Bool_t *bTracks2 = new Bool_t[ntrack2];
  for (Int_t itrack2=0; itrack2<ntrack2; ++itrack2) bTracks2[itrack2]=kFALSE;

  UInt_t selectedMask=(1<<fPairPreFilter.GetCuts()->GetEntries())-1;
  UInt_t selectedMaskPair=(1<<fPairFilter.GetCuts()->GetEntries())-1;

  Int_t nRejPasses = 1; //for fPreFilterUnlikeOnly and no set flag 
  if (fPreFilterAllSigns) nRejPasses = 3;

  for (Int_t iRP=0; iRP < nRejPasses; ++iRP) {
	Int_t arr1RP=arr1, arr2RP=arr2;
	TObjArray *arrTracks1RP=&arrTracks1;
	TObjArray *arrTracks2RP=&arrTracks2;
	Bool_t *bTracks1RP = bTracks1;
	Bool_t *bTracks2RP = bTracks2;
	switch (iRP) {
		case 1: arr1RP=arr1;arr2RP=arr1;
				arrTracks1RP=&arrTracks1;
				arrTracks2RP=&arrTracks1;
				bTracks1RP = bTracks1;
				bTracks2RP = bTracks1;
				break;
		case 2: arr1RP=arr2;arr2RP=arr2;
				arrTracks1RP=&arrTracks2;
				arrTracks2RP=&arrTracks2;
				bTracks1RP = bTracks2;
				bTracks2RP = bTracks2;
				break;
		default: ;//nothing to do
	}
	Int_t ntrack1RP=(*arrTracks1RP).GetEntriesFast();
	Int_t ntrack2RP=(*arrTracks2RP).GetEntriesFast();

	Int_t pairIndex=GetPairIndex(arr1RP,arr2RP);

	for (Int_t itrack1=0; itrack1<ntrack1RP; ++itrack1){
	  Int_t end=ntrack2RP;
	  if (arr1RP==arr2RP) end=itrack1;
	  for (Int_t itrack2=0; itrack2<end; ++itrack2){
		TObject *track1=(*arrTracks1RP).UncheckedAt(itrack1);
		TObject *track2=(*arrTracks2RP).UncheckedAt(itrack2);
		if (!track1 || !track2) continue;
		//create the pair
		candidate.SetTracks(static_cast<AliVTrack*>(track1), fPdgLeg1,
			static_cast<AliVTrack*>(track2), fPdgLeg2);

		candidate.SetType(pairIndex);
		candidate.SetLabel(AliDielectronMC::Instance()->GetLabelMotherWithPdg(&candidate,fPdgMother));
		//relate to the production vertex
		//       if (AliDielectronVarManager::GetKFVertex()) candidate.SetProductionVertex(*AliDielectronVarManager::GetKFVertex());

		//pair cuts
		UInt_t cutMask=fPairPreFilter.IsSelected(&candidate);

		//apply cut
		if (cutMask!=selectedMask) continue;
		if (fCfManagerPair) fCfManagerPair->Fill(selectedMaskPair+1 ,&candidate);
		if (fHistos) FillHistogramsPair(&candidate,kTRUE);
		//set flags for track removal
		bTracks1RP[itrack1]=kTRUE;
		bTracks2RP[itrack2]=kTRUE;
	  }
	}
  }

  //remove the tracks from the Track arrays
  for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1){
    if(bTracks1[itrack1]) arrTracks1.AddAt(0x0, itrack1);
  }
  for (Int_t itrack2=0; itrack2<ntrack2; ++itrack2){
    if(bTracks2[itrack2]) arrTracks2.AddAt(0x0, itrack2);
  }

  // clean up
  delete [] bTracks1;
  delete [] bTracks2;

  //compress the track arrays
  arrTracks1.Compress();
  arrTracks2.Compress();
  
  //apply leg cuts after the pre filter
  if ( fPairPreFilterLegs.GetCuts()->GetEntries()>0 ) {
    selectedMask=(1<<fPairPreFilterLegs.GetCuts()->GetEntries())-1;
    //loop over tracks from array 1
    for (Int_t itrack=0; itrack<arrTracks1.GetEntriesFast();++itrack){
      //test cuts
      UInt_t cutMask=fPairPreFilterLegs.IsSelected(arrTracks1.UncheckedAt(itrack));
      
      //apply cut
      if (cutMask!=selectedMask) arrTracks1.AddAt(0x0,itrack);;
    }
    arrTracks1.Compress();
    
    //in case of like sign don't loop over second array
    if (arr1==arr2) {
      arrTracks2=arrTracks1;
    } else {
      
      //loop over tracks from array 2
      for (Int_t itrack=0; itrack<arrTracks2.GetEntriesFast();++itrack){
      //test cuts
        UInt_t cutMask=fPairPreFilterLegs.IsSelected(arrTracks2.UncheckedAt(itrack));
      //apply cut
        if (cutMask!=selectedMask) arrTracks2.AddAt(0x0,itrack);
      }
      arrTracks2.Compress();
      
    }
  }
  //For unlike-sign monitor track-cuts:
  if (arr1!=arr2&&fHistos) {
    TObjArray *unlikesignArray[2] = {&arrTracks1,&arrTracks2};
    FillHistogramsTracks(unlikesignArray);
  }
}

//________________________________________________________________
void AliDielectron::FillPairArrays(Int_t arr1, Int_t arr2)
{
  //
  // select pairs and fill pair candidate arrays
  //

  TObjArray arrTracks1=fTracks[arr1];
  TObjArray arrTracks2=fTracks[arr2];

  //process pre filter if set
  if ((!fPreFilterAllSigns) && (!fPreFilterUnlikeOnly) && ( fPairPreFilter.GetCuts()->GetEntries()>0 ))  PairPreFilter(arr1, arr2, arrTracks1, arrTracks2);
  
  Int_t pairIndex=GetPairIndex(arr1,arr2);

  Int_t ntrack1=arrTracks1.GetEntriesFast();
  Int_t ntrack2=arrTracks2.GetEntriesFast();

  AliDielectronPair *candidate=new AliDielectronPair;

  UInt_t selectedMask=(1<<fPairFilter.GetCuts()->GetEntries())-1;
  
  for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1){
    Int_t end=ntrack2;
    if (arr1==arr2) end=itrack1;
    for (Int_t itrack2=0; itrack2<end; ++itrack2){
      //create the pair
      candidate->SetTracks(static_cast<AliVTrack*>(arrTracks1.UncheckedAt(itrack1)), fPdgLeg1,
                             static_cast<AliVTrack*>(arrTracks2.UncheckedAt(itrack2)), fPdgLeg2);
      candidate->SetType(pairIndex);
      Int_t label=AliDielectronMC::Instance()->GetLabelMotherWithPdg(candidate,fPdgMother);
      candidate->SetLabel(label);
      if (label>-1) candidate->SetPdgCode(fPdgMother);

      //pair cuts
      UInt_t cutMask=fPairFilter.IsSelected(candidate);
      
      //CF manager for the pair
      if (fCfManagerPair) fCfManagerPair->Fill(cutMask,candidate);

      //apply cut
      if (cutMask!=selectedMask) continue;

      //add the candidate to the candidate array 
      PairArray(pairIndex)->Add(candidate);
      //get a new candidate
      candidate=new AliDielectronPair;
    }
  }
  //delete the surplus candidate
  delete candidate;
}

//________________________________________________________________
void AliDielectron::FillPairArrayTR()
{
  //
  // select pairs and fill pair candidate arrays
  //
  UInt_t selectedMask=(1<<fPairFilter.GetCuts()->GetEntries())-1;
  
  while ( fTrackRotator->NextCombination() ){
    AliDielectronPair candidate;
    candidate.SetTracks(&fTrackRotator->GetKFTrackP(), &fTrackRotator->GetKFTrackN(),
                        fTrackRotator->GetVTrackP(),fTrackRotator->GetVTrackN());
    candidate.SetType(kEv1PMRot);
    
    //pair cuts
    UInt_t cutMask=fPairFilter.IsSelected(&candidate);
    
    //CF manager for the pair
    if (fCfManagerPair) fCfManagerPair->Fill(cutMask,&candidate);
    
    //apply cut
    if (cutMask==selectedMask) {
     if(fHistos) FillHistogramsPair(&candidate);
     if(fStoreRotatedPairs) PairArray(kEv1PMRot)->Add(new AliDielectronPair(candidate));
    } 
  }
}

//________________________________________________________________
void AliDielectron::FillDebugTree()
{
  //
  // Fill Histogram information for tracks and pairs
  //
  
  //Fill Debug tree
  for (Int_t i=0; i<10; ++i){
    Int_t ntracks=PairArray(i)->GetEntriesFast();
    for (Int_t ipair=0; ipair<ntracks; ++ipair){
      fDebugTree->Fill(static_cast<AliDielectronPair*>(PairArray(i)->UncheckedAt(ipair)));
    }
  }
}

//________________________________________________________________
void AliDielectron::SaveDebugTree()
{
  //
  // delete the debug tree, this will also write the tree
  //
  if (fDebugTree) fDebugTree->DeleteStreamer();
}


//__________________________________________________________________
void AliDielectron::AddSignalMC(AliDielectronSignalMC* signal) {
  //
  //  Add an MC signal to the signals list
  //
  if(!fSignalsMC) {
    fSignalsMC = new TObjArray();
    fSignalsMC->SetOwner();
  }
  fSignalsMC->Add(signal);
}
