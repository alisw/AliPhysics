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

#include <AliESDEvent.h>
#include <AliESDtrack.h>

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
  fPdgMother(443),
  fPdgLeg1(11),
  fPdgLeg2(11),
  fNoPairing(kFALSE),
  fHistos(0x0),
  fPairCandidates(new TObjArray(10)),
  fCfManagerPair(0x0),
  fTrackRotator(0x0),
  fDebugTree(0x0),
  fPreFilterUnlikeOnly(kFALSE),
  fHasMC(kFALSE)
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
  fPdgMother(443),
  fPdgLeg1(11),
  fPdgLeg2(11),
  fNoPairing(kFALSE),
  fHistos(0x0),
  fPairCandidates(new TObjArray(10)),
  fCfManagerPair(0x0),
  fTrackRotator(0x0),
  fDebugTree(0x0),
  fPreFilterUnlikeOnly(kFALSE),
  fHasMC(kFALSE)
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
}

//________________________________________________________________
void AliDielectron::Init()
{
  //
  // Initialise objects
  //

  if(GetHasMC()) AliDielectronMC::Instance()->SetHasMC(GetHasMC());
  
  if (fCfManagerPair) fCfManagerPair->InitialiseContainer(fPairFilter);
  if (fTrackRotator)  {
    fTrackRotator->SetTrackArrays(&fTracks[0],&fTracks[1]);
    fTrackRotator->SetPdgLegs(fPdgLeg1,fPdgLeg2);
  }
  if (fDebugTree) fDebugTree->SetDielectron(this);
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
   
  //in case we have MC load the MC event and process the MC particles
  if (AliDielectronMC::Instance()->HasMC()) {
    if (!AliDielectronMC::Instance()->ConnectMCEvent()){
      AliError("Could not properly connect the MC event, skipping this event!");
      return;
    }
    ProcessMC();
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
  
  AliDielectronVarManager::SetEvent(ev1);
  
  //fill track arrays for the first event
  if (ev1){
    FillTrackArrays(ev1);
    if ((fPreFilterUnlikeOnly) && ( fPairPreFilter.GetCuts()->GetEntries()>0 )) PairPreFilter(0, 1, fTracks[0], fTracks[1]);
  }


  //fill track arrays for the second event
  if (ev2) {
    FillTrackArrays(ev2,1);
    if ((fPreFilterUnlikeOnly) && ( fPairPreFilter.GetCuts()->GetEntries()>0 )) PairPreFilter(2, 3, fTracks[2], fTracks[3]);
  }

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
  
  //in case there is a histogram manager, fill the QA histograms
  if (fHistos) FillHistograms(ev1);

  //fill debug tree if a manager is attached
  if (fDebugTree) FillDebugTree();
}

//________________________________________________________________
void AliDielectron::ProcessMC()
{
  //
  // Process the MC data
  //

  //loop over all MC data and Fill the CF container if it exist
  if (!fCfManagerPair) return;
  fCfManagerPair->SetPdgMother(fPdgMother);
  AliDielectronMC *dieMC=AliDielectronMC::Instance();
  for (Int_t ipart=0; ipart<dieMC->GetNMCTracks();++ipart){
    //TODO: MC truth cut properly!!!
    AliVParticle *mcPart=dieMC->GetMCTrackFromMCEvent(ipart);
    if (!dieMC->IsMCMotherToEE(mcPart, fPdgMother)) continue;
    fCfManagerPair->FillMC(mcPart);
  }
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
void AliDielectron::FillHistograms(const AliVEvent *ev)
{
  //
  // Fill Histogram information for tracks and pairs
  //
  
  TString  className,className2;
  Double_t values[AliDielectronVarManager::kNMaxValues];
  //Fill event information
  AliDielectronVarManager::Fill(ev, values);
  if (fHistos->GetHistogramList()->FindObject("Event"))
    fHistos->FillClass("Event", AliDielectronVarManager::kNMaxValues, values);
  
  //Fill track information, separately for the track array candidates
  for (Int_t i=0; i<4; ++i){
    className.Form("Track_%s",fgkTrackClassNames[i]);
    if (!fHistos->GetHistogramList()->FindObject(className.Data())) continue;
    Int_t ntracks=fTracks[i].GetEntriesFast();
    for (Int_t itrack=0; itrack<ntracks; ++itrack){
      AliDielectronVarManager::Fill(fTracks[i].UncheckedAt(itrack), values);
      fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
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
void AliDielectron::PairPreFilter(Int_t arr1, Int_t arr2, TObjArray &arrTracks1, TObjArray &arrTracks2)
{
  //
  // Prefilter tracks from pairs
  // Needed for datlitz rejections
  // remove all tracks from the Single track arrays that pass the cuts in this filter
  //
  Int_t pairIndex=GetPairIndex(arr1,arr2);
  
  Int_t ntrack1=arrTracks1.GetEntriesFast();
  Int_t ntrack2=arrTracks2.GetEntriesFast();
  
  AliDielectronPair candidate;
  
  UInt_t selectedMask=(1<<fPairPreFilter.GetCuts()->GetEntries())-1;
  UInt_t selectedMaskPair=(1<<fPairFilter.GetCuts()->GetEntries())-1;
  
  for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1){
    Int_t end=ntrack2;
    if (arr1==arr2) end=itrack1;
    Bool_t accepted=kFALSE;
    for (Int_t itrack2=0; itrack2<end; ++itrack2){
      AliVTrack *track1=static_cast<AliVTrack*>(arrTracks1.UncheckedAt(itrack1));
      AliVTrack *track2=static_cast<AliVTrack*>(arrTracks2.UncheckedAt(itrack2));
      if (!track1 || !track2) continue;
      //create the pair
      candidate.SetTracks(track1, fPdgLeg1,
                          track2, fPdgLeg2);
      candidate.SetType(pairIndex);
      candidate.SetLabel(AliDielectronMC::Instance()->GetLabelMotherWithPdg(&candidate,fPdgMother));
      
      //pair cuts
      UInt_t cutMask=fPairPreFilter.IsSelected(&candidate);
      
      //apply cut
      if (cutMask!=selectedMask) continue;
      if (fCfManagerPair) fCfManagerPair->Fill(selectedMaskPair+1 ,&candidate);
      accepted=kTRUE;
      FillHistogramsPair(&candidate,kTRUE);
      //remove the tracks from the Track arrays
      arrTracks2.AddAt(0x0,itrack2);
      //in case of like sign remove the track from both arrays!
      if (arr1==arr2) arrTracks1.AddAt(0x0, itrack2);
    }
    if ( accepted ) arrTracks1.AddAt(0x0,itrack1);
  }
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
  if (arr1!=arr2) {
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
  if ((!fPreFilterUnlikeOnly) && ( fPairPreFilter.GetCuts()->GetEntries()>0 ))  PairPreFilter(arr1, arr2, arrTracks1, arrTracks2);
  
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
      candidate->SetLabel(AliDielectronMC::Instance()->GetLabelMotherWithPdg(candidate,fPdgMother));

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
    if (cutMask==selectedMask&&fHistos) FillHistogramsPair(&candidate);
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

