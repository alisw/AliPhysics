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
#include "AliMultSelection.h"
#include "AliESDv0.h"
#include "AliCentrality.h"
#include "AliLog.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDarrayADC.h"

#include "TChain.h"
#include "TFile.h"
#include "TRandom.h"
#include "THnSparse.h"


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
    : AliTRDdigitsTask(name)
{

  // V0 Kine cuts
  fV0cuts = new AliESDv0KineCuts();

  SetDigitsInputFilename("TRD.Digits.root");
  SetDigitsOutputFilename("TRD.FltDigits.root");

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

}


//_________________________________________________
AliTRDdigitsFilter::~AliTRDdigitsFilter()
{
  //
  // Destructor
  //

  delete fV0cuts;
}

//________________________________________________________________________
void AliTRDdigitsFilter::AcceptTracks(TString label, EPID_t pid,
    Float_t minPt, Float_t maxPt, Float_t fraction)
{

  TrackCrit crit;

  crit.fLabel     = label;
  crit.fPid       = pid;
  crit.fMinPt     = minPt;
  crit.fMaxPt     = maxPt;
  crit.fFraction  = fraction;

  fTrackCriteria.push_back(crit);
  AliInfoF("nTrackCrit: %lu", fTrackCriteria.size());
}

//________________________________________________________________________
void AliTRDdigitsFilter::AcceptEvents(TString label,
    Float_t minCent, Float_t maxCent, Float_t fraction)
{

  EventCrit crit;

  crit.fLabel     = label;
  crit.fMinCent     = minCent;
  crit.fMaxCent     = maxCent;
  crit.fFraction  = fraction;

  fEventCriteria.push_back(crit);
  AliInfoF("nEventCrit: %lu", fEventCriteria.size());
}

//________________________________________________________________________
void AliTRDdigitsFilter::PrintSettings()
{

  AliInfo("Accept these particle classes:");
  for (std::vector<TrackCrit>::iterator iCrit = fTrackCriteria.begin();
       iCrit != fTrackCriteria.end(); iCrit++) {

         AliInfoF("%10s: pidclass=%d, %f < pT(GeV/c) < %f",
         iCrit->fLabel.Data(), iCrit->fPid, iCrit->fMinPt, iCrit->fMaxPt
       );
   }


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
  fOutputList=new TList;
  fOutputList->SetOwner();

  CreateV0Plots();

  // statistics of accepted events
  int nClasses = fEventCriteria.size() + fTrackCriteria.size();
  fhEventCuts = new TH1F("EventCuts","statistics of event cuts",
    nClasses+4, 0., nClasses+4.0);

  fhEventCuts->GetXaxis()->SetBinLabel(1,"event_in");
  fhEventCuts->GetXaxis()->SetBinLabel(2,"event_esd");
  fhEventCuts->GetXaxis()->SetBinLabel(3,"event_vtx");
  fhEventCuts->GetXaxis()->SetBinLabel(4,"event_acc");

  for (Int_t i = 0; i<fEventCriteria.size(); i++) {
    fhEventCuts->GetXaxis()->SetBinLabel(i+5, "event_acc_"+fEventCriteria[i].fLabel);
  }

  for (Int_t i = 0; i<fTrackCriteria.size(); i++) {
    fhEventCuts->GetXaxis()->SetBinLabel( i+5+fEventCriteria.size(),
    "event_acc_"+fTrackCriteria[i].fLabel);
  }


  // THnSparse for accepted tracks
  const Int_t ntc = 2<<fTrackCriteria.size();
  Int_t nbins[]   = { ntc,    64,   10,  10 };
  Double_t xmin[] = { 0.0,   0.0, -1.0, 0.0 };
  Double_t xmax[] = { float(ntc),  64.0,  1.0, 2*TMath::Pi() };

  Double_t pTedges[] = {
    1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
    2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
    3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
    5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8,
    7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
    10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
    20., 25., 30., 35., 40., 45., 50., 55., 60.
  };

  fhAcc = new THnSparseF("fhAcc", "accepted tracks", 4, nbins, xmin, xmax);
  fhAcc->SetBinEdges(1,pTedges);

  //fhAcc->GetAxis(0)->SetBinLabel(1,"foo");
  fhAcc->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  fhAcc->GetAxis(2)->SetTitle("#eta");
  fhAcc->GetAxis(3)->SetTitle("#phi (rad)");

  // pT histograms
  fhPtTag  = new TH1F("fhPtTag",  "pT of PID-tagged tracks", 100,0.,20.);
  fhPtGood = new TH1F("fhPtGood", "pT after quality cuts",   100,0.,20.);
  fhPtAcc  = new TH1F("fhPtAcc",  "pT of accepted tracks",   100,0.,20.);

  // multiplicity histograms
  fhCent = new TH1F("fhCentralityAll", "Centrality of Events", 105, 0., 105.);
  fhCentAcc = new TH1F("fhCentralityAccepted", "Centrality of Accepted Events",
                       105, 0., 105.);

  // add everything to the list
  fOutputList->Add(fhAcc);
  fOutputList->Add(fhEventCuts);
  fOutputList->Add(fhPtTag);
  fOutputList->Add(fhPtGood);
  fOutputList->Add(fhPtAcc);
  fOutputList->Add(fhCent);
  fOutputList->Add(fhCentAcc);

  CreateTriggerHistos();

  PostData(1,fOutputList);


}

//_____________________________________________________________________________
void AliTRDdigitsFilter::UserExec(Option_t *)
{
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
  // IMPORTANT: call NextEvent() for book-keeping
  // -----------------------------------------------------------------
  NextEvent();
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------


  //
  //calls the Process function
  //

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) {
    printf("ERROR: Could not get ESDInputHandler \n");
    fESDevent = NULL;
  } else {
    fESDevent = (AliESDEvent *) esdH->GetEvent();
  }

  // allocate space for the PID tags
  fPidTags.resize(fESDevent->GetNumberOfTracks());

  FillV0PIDlist();
  Process();

  PostData(1,fOutputList);

}



//________________________________________________________________________
void AliTRDdigitsFilter::Process()
{
  //
  //called for each event
  //

  fhEventCuts->Fill("event_in",1);
  FillTriggerHisto(fhTrgAll);

  if (!fESDevent) {
    Printf("ERROR: fESDevent not available");
    return;
  }

  fhEventCuts->Fill("event_esd",1);

  // check for a valid event vertex
  const AliESDVertex* fESDEventvertex = fESDevent->GetPrimaryVertexTracks();
  if (!fESDEventvertex) return;

  Int_t ncontr = fESDEventvertex->GetNContributors();
  if (ncontr <= 0) return;

  fhEventCuts->Fill("event_vtx",1);


  //-------------------------------------------------------------------
  Int_t keepEvent = 0;
  Int_t keepNTracks = 0;

  Int_t nEventCrit = fEventCriteria.size();
  Int_t nTrackCrit = fTrackCriteria.size();

  //-------------------------------------------------------------------

  AliMultSelection *multSelection =
  static_cast<AliMultSelection*>(fESDevent->FindListObject("MultSelection"));

  if(multSelection) {

    Float_t centrality = multSelection->GetMultiplicityPercentile("V0M");

    fhCent->Fill(centrality);

    for (Int_t i = 0; i<nEventCrit; i++) {

      // check event criteria
      if ( fEventCriteria[i].fMinCent > centrality) continue;
      if ( fEventCriteria[i].fMaxCent <= centrality) continue;

      // accept given fraction of events
      if ( gRandom->Uniform() > fEventCriteria[i].fFraction ) continue;

      keepEvent |= 1<<i;

   }

   if (keepEvent) {
     fhCentAcc->Fill(centrality);
   }

 }

 //-------------------------------------------------------------------

 for (int iTrack=0; iTrack<fPidTags.size(); iTrack++) {

    Int_t keepTrack = 0;

    if (fPidTags[iTrack] == kPidUndef) continue;
    if (fPidTags[iTrack] == kPidError) continue;

    AliESDtrack* track = fESDevent->GetTrack(iTrack);

    fhPtTag->Fill(track->Pt());

    // general track cuts
    if ( ! PassTrackCuts(track,3) ) continue;

    fhPtGood->Fill(track->Pt());

    for (Int_t i = 0; i<nTrackCrit; i++) {

      // check track criteria
      if ( fTrackCriteria[i].fPid != fPidTags[iTrack] ) continue;
      if ( fTrackCriteria[i].fMinPt > track->Pt()) continue;
      if ( fTrackCriteria[i].fMaxPt < track->Pt()) continue;

      // accept given fraction of tracks
      if ( gRandom->Uniform() > fTrackCriteria[i].fFraction ) continue;

      keepEvent |= 1<<(nEventCrit + i);
      keepTrack |= 1<<i;
      keepNTracks++;
    }

    if (keepTrack) {
      Double_t data[4];
      data[0] = keepTrack;
      data[1] = track->Pt();
      data[2] = track->Eta();
      data[3] = track->Phi();

      fhAcc->Fill(data);
      fhPtAcc->Fill(track->Pt());
    }

  }

  if (keepEvent) {

    fhEventCuts->Fill("event_acc",1);
    FillTriggerHisto(fhTrgAcc);

    for (Int_t i = 0; i<nEventCrit; i++) {
      if (keepEvent & ( 1 << i ) ) {
        fhEventCuts->Fill("event_acc_"+fEventCriteria[i].fLabel,1);
      }
    }

    for (Int_t i = 0; i<nTrackCrit; i++) {
      if (keepEvent & ( 1 << (nEventCrit+i) ) ) {
        fhEventCuts->Fill("event_acc_"+fTrackCriteria[i].fLabel,1);
      }
    }

    if (ReadDigits()) {
      WriteDigits();
    }

  }

  PostData(1,fOutputList);
}

// //________________________________________________________________________
// Bool_t AliTRDdigitsFilter::ReadDigits()
// {
//   if (!fDigitsInputFile) {
//     AliError("Digits file not open");
//     return kFALSE;
//   }
//
//   TTree* tr = (TTree*)fDigitsInputFile->Get(Form("Event%d/TreeD",
// 						 fEventNoInFile));
//
//   if (!tr) {
//     AliErrorF("Digits tree for event %d not found", fEventNoInFile);
//     return kFALSE;
//   }
//
//   for (Int_t det=0; det<540; det++) {
//     fDigMan->ClearArrays(det);
//     fDigMan->ClearIndexes(det);
//   }
//
//   fDigMan->ReadDigits(tr);
//   delete tr;
//   return kTRUE;
// }
//
// //________________________________________________________________________
// void AliTRDdigitsFilter::WriteDigits()
// {
//   if (!fDigitsOutputFile) {
//     AliError("Filtered digits file not open");
//     return;
//   }
//
//   TDirectory* evdir =
//     fDigitsOutputFile->mkdir(Form("Event%d", fEventNoInFile),
// 			     Form("Event%d", fEventNoInFile));
//
//   evdir->Write();
//   evdir->cd();
//
//   TTree* tr = new TTree("TreeD", "TreeD");
//   fDigMan->MakeBranch(tr);
//   fDigMan->WriteDigits();
//   delete tr;
// }


// //________________________________________________________________________
// void AliTRDdigitsFilter::Terminate(const Option_t *)
// {
//     //
//     // Terminate function
//     //
// }


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
