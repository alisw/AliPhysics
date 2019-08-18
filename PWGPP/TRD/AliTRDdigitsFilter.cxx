/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author:             *
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
// The task:
// write out raw
//
//
// Author:
//
//


#include "AliTRDdigitsFilter.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDv0KineCuts.h"
#include "AliESDv0.h"
#include "AliCentrality.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDarrayADC.h"

#include "TChain.h"
#include "TFile.h"
#include "TRandom.h"

class TCanvas;
class TAxis;
class TFile;
class TStyle;
class TString;
class TH1F;
class TH2D;
class THnSparse;
class TLegend;
class TVirtualFitter;
class AliESDtrackCuts;
class AliStack;
class AliMCParticle;


using namespace std;

ClassImp(AliTRDdigitsFilter)



//________________________________________________________________________
AliTRDdigitsFilter::AliTRDdigitsFilter(const char *name)
    : AliAnalysisTaskSE(name), fV0cuts(0x0), fPidTags(0x0),
    fESDEvent(0), fOutputContainer(0), fESDtrackCuts(0),
    fESDtrackCutsV0(0), fListQA(0x0),
    fDigitsInputFile(0), fDigitsOutputFile(0),
    fEventNoInFile(-1), fDigMan(0)
{
  //
  // Constructor
  //

  fDigMan = new AliTRDdigitsManager;
  fDigMan->CreateArrays();

  // V0 Kine cuts
  fV0cuts = new AliESDv0KineCuts();

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

}


//_________________________________________________
AliTRDdigitsFilter::~AliTRDdigitsFilter()
{
  //
  // Destructor
  //

  delete fDigitsInputFile;
  delete fDigitsOutputFile;
  delete fDigMan;
  delete fV0cuts;
}

void AliTRDdigitsFilter::AcceptParticles(TString label, EPID_t pid,
    Float_t minPt, Float_t maxPt, Float_t fraction)
{

  AcceptCrit crit;

  crit.fLabel     = label;
  crit.fPid       = pid;
  crit.fMinPt     = minPt;
  crit.fMaxPt     = maxPt;
  crit.fFraction  = fraction;

  fAcceptCriteria.push_back(crit);
}

//________________________________________________________________________
void AliTRDdigitsFilter::UserCreateOutputObjects()
{
  //
  // Definition of user output ntuple and histogram file
  //

  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!inputHandler)
    printf("Inputhandler not available \n");


  OpenFile(1);
  fListQA=new TList;
  fListQA->SetOwner();

  // V0 QA histograms
  fhArmenteros  = new TH2F( "Armenteros","Armenteros plot",
                            200,-1.,1.,200,0.,0.4);

  // statistics of accepted events
  fhEventCuts = new TH1F("EventCuts","statistics of event cuts",10,0.,10.);

  // pT histograms
  fhPtTag  = new TH1F("fhPtTag",  "pT of PID-tagged tracks", 100,0.,20.);
  fhPtGood = new TH1F("fhPtGood", "pT after quality cuts",   100,0.,20.);
  fhPtAcc  = new TH1F("fhPtAcc",  "pT of accepted tracks",   100,0.,20.);

  // add everything to the list
  fListQA->Add(fhArmenteros);
  fListQA->Add(fhEventCuts);
  fListQA->Add(fhPtTag);
  fListQA->Add(fhPtGood);
  fListQA->Add(fhPtAcc);

  PostData(1,fListQA);


}

//_____________________________________________________________________________
Bool_t AliTRDdigitsFilter::UserNotify()
{
  delete fDigitsInputFile;
  delete fDigitsOutputFile;

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  TString ofname = esdH->GetTree()->GetCurrentFile()->GetName();
  TString ifname = ofname;

  ifname.ReplaceAll("AliESDs.root", "TRD.Digits.root");
  ofname.ReplaceAll("AliESDs.root", "TRD.FltDigits.root");

  fDigitsInputFile  = new TFile(ifname);
  fDigitsOutputFile = new TFile(ofname,"RECREATE");

  fEventNoInFile = 0;

  return kTRUE;
}

//_____________________________________________________________________________
void AliTRDdigitsFilter::UserExec(Option_t *)
{
    //
    //calls the Process function
    //

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      printf("ERROR: Could not get ESDInputHandler \n");
    }
    else fESDEvent = (AliESDEvent *) esdH->GetEvent();

    // allocate space for the PID tags
    fPidTags.resize(fESDEvent->GetNumberOfTracks());

    FillV0PIDlist();
    Process(fESDEvent);

    PostData(1,fListQA);

    // increment the event counter for this file
    fEventNoInFile++;

}



//________________________________________________________________________
void AliTRDdigitsFilter::Process(AliESDEvent *const esdEvent)
{
  //
  //called for each event
  //

  fhEventCuts->Fill("event_in",1);

  if (!esdEvent) {
    Printf("ERROR: esdEvent not available");
    return;
  }

  fhEventCuts->Fill("event_esd",1);

  // check for a valid event vertex
  const AliESDVertex* fESDEventvertex = esdEvent->GetPrimaryVertexTracks();
  if (!fESDEventvertex)
    return;

  Int_t ncontr = fESDEventvertex->GetNContributors();
  if (ncontr <= 0) return;

  fhEventCuts->Fill("event_vtx",1);

  Bool_t keepEvent = kFALSE;

  for (int iTrack=0; iTrack<fPidTags.size(); iTrack++) {

    Bool_t keepTrack = kFALSE;

    if (fPidTags[iTrack] == kPidUndef) continue;
    if (fPidTags[iTrack] == kPidError) continue;

    AliESDtrack* track = esdEvent->GetTrack(iTrack);

    fhPtTag->Fill(track->Pt());

    // general track cuts
    if ( ! PassTrackCuts(track,3) ) continue;

    fhPtGood->Fill(track->Pt());

    for (std::list<AcceptCrit>::iterator iCrit = fAcceptCriteria.begin();
         iCrit != fAcceptCriteria.end(); iCrit++) {

        // check track criteria
        if ( iCrit->fPid != fPidTags[iTrack] ) continue;
        if ( iCrit->fMinPt > track->Pt()) continue;
        if ( iCrit->fMaxPt < track->Pt()) continue;

        // accept given fraction of tracks
        if ( gRandom->Uniform() > iCrit->fFraction ) continue;

        keepEvent = kTRUE;
        keepTrack = kTRUE;
     }

     if (keepTrack) {
       fhPtAcc->Fill(track->Pt());
     }

  }

  if (keepEvent) {

    fhEventCuts->Fill("event_acc",1);

    // load the digits from TRD.Digits.root
    ReadDigits();

    // store the digits in TRD.FltDigits.root
    WriteDigits();
  }

  PostData(1,fListQA);
}


//________________________________________________________________________
void AliTRDdigitsFilter::ReadDigits()
{
  TTree* tr = (TTree*)fDigitsInputFile->Get(Form("Event%d/TreeD",
						 fEventNoInFile));
  for (Int_t det=0; det<540; det++) {
    fDigMan->ClearArrays(det);
    fDigMan->ClearIndexes(det);
  }

  fDigMan->ReadDigits(tr);
  delete tr;
}

//________________________________________________________________________
void AliTRDdigitsFilter::WriteDigits()
{
  TDirectory* evdir =
    fDigitsOutputFile->mkdir(Form("Event%d", fEventNoInFile),
			     Form("Event%d", fEventNoInFile));

  evdir->Write();
  evdir->cd();

  TTree* tr = new TTree("TreeD", "TreeD");
  fDigMan->MakeBranch(tr);
  fDigMan->WriteDigits();
  delete tr;
}


//________________________________________________________________________
void AliTRDdigitsFilter::Terminate(const Option_t *)
{
    //
    // Terminate function
    //
}


//______________________________________________________________________________
void AliTRDdigitsFilter::FillV0PIDlist(){

  //
  // Fill the PID object arrays holding the pointers to identified particle tracks
  //

  // Dynamic cast to ESD events (DO NOTHING for AOD events)
  AliESDEvent *event = dynamic_cast<AliESDEvent *>(InputEvent());
  if ( !event )  return;


  // V0 selection
  // set event
  fV0cuts->SetEvent(event);


  // loop over V0 particles
  for(Int_t iv0=0; iv0<event->GetNumberOfV0s();iv0++){

    AliESDv0 *v0 = (AliESDv0 *) event->GetV0(iv0);

    if(!v0) continue;
    if(v0->GetOnFlyStatus()) continue;

    // Get the particle selection
    Bool_t foundV0 = kFALSE;
    Int_t pdgV0, pdgP, pdgN;
    foundV0 = fV0cuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if(!foundV0) continue;
    Int_t iTrackP = v0->GetPindex();  // positive track
    Int_t iTrackN = v0->GetNindex();  // negative track

    // v0 Armenteros plot (QA)
    Float_t armVar[2] = {0.0,0.0};
    fV0cuts->Armenteros(v0, armVar);
    // if ( !(TMath::Power(armVar[0]/0.95,2)+TMath::Power(armVar[1]/0.05,2) < 1) ) continue;

    if(fListQA&&fhArmenteros) fhArmenteros->Fill(armVar[0],armVar[1]);

    // fill the tags

    if( pdgP ==   -11 ) { fPidTags[iTrackP] = kPidV0Electron; }
    if( pdgN ==    11 ) { fPidTags[iTrackN] = kPidV0Electron; }

    if( pdgP ==   211 ) { fPidTags[iTrackP] = kPidV0Pion; }
    if( pdgN ==  -211 ) { fPidTags[iTrackN] = kPidV0Pion; }

    if( pdgP ==  2212 ) { fPidTags[iTrackP] = kPidV0Proton; }
    if( pdgN == -2212 ) { fPidTags[iTrackN] = kPidV0Proton; }

  }
}


//________________________________________________________________________
Bool_t AliTRDdigitsFilter::PassTrackPIDCuts(AliESDtrack *fESDTrack)
{
    //
    // check if tracks pass minimum quality critieria
    //
    if(!fESDTrack) return kFALSE;
    if(fESDTrack->GetTPCsignal()<85) return kFALSE;
    if(fESDTrack->GetTPCsignal()>115) return kFALSE;
    return kTRUE;
}

//________________________________________________________________________
Bool_t AliTRDdigitsFilter::PassTrackCuts(AliESDtrack *fESDTrack, Int_t threshold)
{
    //
    // check if tracks pass minimum quality critieria
    //

    if(!fESDTrack) return kFALSE;

    // DCA to PV
    Float_t dca[2];
    fESDTrack->GetImpactParameters(dca[0],dca[1]);
    if(dca[0]>5||dca[1]>10) return kFALSE;

    // eta cut
    if((TMath::Abs(fESDTrack->Eta()))>0.9) return kFALSE;

    //TRD out
    if((fESDTrack->GetStatus()&AliVTrack::kTRDout)==0)return kFALSE;

    // TPC refit
    if((fESDTrack->GetStatus()&AliVTrack::kTPCrefit)==0)return kFALSE;
    // remove kinks
    if(fESDTrack->GetKinkIndex(0)>0) return kFALSE;

    Float_t tpcchi2=99;
    Int_t tpcnclusF=fESDTrack->GetTPCNclsF();
    if(tpcnclusF!=0) tpcchi2=(Float_t)fESDTrack->GetTPCchi2()/tpcnclusF;
    else tpcchi2=1000;
    if(tpcchi2 > 4) return kFALSE;


    Int_t ntrackletstracking=fESDTrack->GetTRDntracklets();
    if(ntrackletstracking<threshold) return kFALSE;

     // QA #TRD PID tracklets
    if(fESDTrack->GetTRDntrackletsPID()<threshold) return kFALSE;

    Int_t ntrl=0;
    for(Int_t jPl=0;jPl<6;jPl++){
      Double_t signal=0;
      for(int isl= 0; isl<= 8;isl++){
        Double_t sigsl=fESDTrack->GetTRDslice(jPl,isl);
        if(sigsl>0)signal+=sigsl;
      }
      // if signal is missing, stop counting
      if(signal<=0||fESDTrack->GetTRDmomentum(jPl)<=0)break;
      ntrl++;
    }
    if(ntrl<threshold) return kFALSE;

    return kTRUE;
}
