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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskCreateNUE.h"
#include "AliMultSelection.h"
#include <vector>

class AliAnalysisTaskCreateNUE; // your analysis class

using namespace std; // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskCreateNUE) // classimp: necessary for root

    AliAnalysisTaskCreateNUE::AliAnalysisTaskCreateNUE() : AliAnalysisTaskSE(),
                                                           fAOD(0), fOutputList(0), fIsMC(false),
  fPeriod("LHC17"),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fMinPt(0.2), fMaxPt(3.0),
  fEtaCut(0.8),
  hEventCount(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCreateNUE::AliAnalysisTaskCreateNUE(const char *name) : AliAnalysisTaskSE(name),
                                                                       fAOD(0), fOutputList(0), fIsMC(false), fPeriod("LHC17"),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fMinPt(0.2), fMaxPt(3.0),
  fEtaCut(0.8),
  hEventCount(0)
{
  // constructor
  DefineInput(0, TChain::Class()); // define the input of the analysis: in this case we take a 'chain' of events
                                   // this chain is created by the analysis manager, so no need to worry about it,
                                   // it does its work automatically
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this case it's a list of histograms
                                   // you can add more output objects by calling DefineOutput(2, classname::Class())
                                   // if you add more output objects, make sure to call PostData for all of them, and to
                                   // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskCreateNUE::~AliAnalysisTaskCreateNUE()
{
  // destructor
  if (fOutputList)
  {
    delete fOutputList; // at the end of your task, it is deleted from memory by calling this function
  }
  if (fGFWSelection)
    delete fGFWSelection;
  if (fGFWSelection15o)
    delete fGFWSelection15o;
}
//_____________________________________________________________________________
void AliAnalysisTaskCreateNUE::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output file
  //
  fOutputList = new TList(); // this is a list which will contain all of your histograms
                             // at the end of the analysis, the contents of this list are written
                             // to the output file
  fOutputList->SetName("WeightList");
  fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all objects it contains and will delete them
                                // if requested (dont worry about this now)
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1 and LHC17n
    fGFWSelection15o = new AliGFWNFCuts();
    fGFWSelection15o->PrintSetup();
  } else {
    fGFWSelection = new AliGFWMCuts();
    fGFWSelection->PrintSetup();
  }

  //avoid changing the index of weight
  hEventCount = new TH1D("hEventCount", "; centrality;;", 5, 0, 5);
	fOutputList->Add(hEventCount);

  double ptBins[] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,
     1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
     2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0};
  hPtMCGen = new TH1D("hMCGen", "histogram of MC generation", 36, ptBins);
  hPtMCGen->Sumw2();
	fOutputList->Add(hPtMCGen);

  hPtMCRec = new TH1D("hMCRec", "histogram of MC reconstruction",36, ptBins);
  hPtMCRec->Sumw2();
	fOutputList->Add(hPtMCRec);

  // don't forget to add it to the list! the list will be written to file, so if you want
  // your histogram in the output file, add it to the list!

  PostData(1, fOutputList); // postdata will notify the analysis manager of changes / updates to the
                            // fOutputList object. the manager will in the end take care of writing your output to file
                            // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskCreateNUE::UserExec(Option_t *)
{
  // user exec
  hEventCount->GetXaxis()->SetBinLabel(1,"Loop Number");
  hEventCount->Fill(0.5);

  fAOD = dynamic_cast<AliAODEvent *>(InputEvent()); // get an event (called fAOD) from the input file
  if (!fAOD) return;
  hEventCount->GetXaxis()->SetBinLabel(2,"AOD OK");
  hEventCount->Fill(1.5);
  //Standard AliEvent Cuts
  if(fUseHM){
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
  }
  else{
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kINT7, true);
  }

  //Standard AliEventCuts for events
  if (!fEventCuts.AcceptEvent(fAOD))
  { // automatic event selection for Run2
    PostData(1, fOutputList);
    return;
  }
  hEventCount->GetXaxis()->SetBinLabel(3,"After fEventCuts");
  hEventCount->Fill(2.5);

  //AliGFWCuts for Sysmatics
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1
    fGFWSelection15o->ResetCuts();
    fGFWSelection15o->SetupCuts(fCurrSystFlag);
    if (!fGFWSelection15o->AcceptVertex(fAOD))
    {
      PostData(1, fOutputList);
      return;
    }
  } else {
    fGFWSelection->ResetCuts();
    fGFWSelection->SetupCuts(fCurrSystFlag);
    if (!fGFWSelection->AcceptVertex(fAOD))
      {
        PostData(1, fOutputList);
        return;
    }
  }
  hEventCount->GetXaxis()->SetBinLabel(4,"After AliGFWCuts");
  hEventCount->Fill(3.5);


  Int_t iTracks(fAOD->GetNumberOfTracks()); // see how many tracks there are in the event

  //..for DCA
  Double_t pos[3], vz, vx, vy;
  vz = fAOD->GetPrimaryVertex()->GetZ();
  vx = fAOD->GetPrimaryVertex()->GetX();
  vy = fAOD->GetPrimaryVertex()->GetY();
  double vtxp[3] = {vx, vy, vz};

  // Start looping over the Reco tracks
  //..LOOP OVER TRACKS........
  //........................................
  const int nAODTracks = fAOD->GetNumberOfTracks();
  for(Int_t nt = 0; nt < nAODTracks; nt++) {

    AliAODTrack *aodTrk = (AliAODTrack*) fInputEvent->GetTrack(nt);

    if (!aodTrk) {
      continue;
    }

    aodTrk->GetXYZ(pos);
    double dcaX = pos[0] - vtxp[0]; 
    double dcaY = pos[1] - vtxp[1];
    double dcaZ = abs(pos[2] - vtxp[2]);
    double dcaXY = TMath::Sqrt(dcaX*dcaX+dcaY*dcaY);

    if (!AcceptAODTrack(aodTrk, pos, vtxp)) continue;
    hPtMCRec->Fill(aodTrk->Pt());
 }




  // Start looping over the Truth tracks
  //..LOOP OVER TRACKS........
  //........................................
  TClonesArray* farray = (TClonesArray*)fAOD->FindListObject("mcparticles");
  const int nAODTracksMC = farray->GetEntries();
  for(Int_t nt = 0; nt < nAODTracksMC; nt++) {

    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(nt, fMCEvent)) continue;
    AliAODMCParticle *track = (AliAODMCParticle*) farray->At(TMath::Abs(nt));

    if (!track) {
      continue;
    }

    // track->GetXYZ(pos);
    if (!AcceptMCTruthTrack(track)) continue;
    hPtMCGen->Fill(track->Pt());
  }

  hEventCount->GetXaxis()->SetBinLabel(5,"Final pass");
  hEventCount->Fill(4.5);
  PostData(1, fOutputList); // stream the results the analysis of this event to
                            // the output manager which will take care of writing
                            // it to a file
}
//_____________________________________________________________________________
void AliAnalysisTaskCreateNUE::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________

Bool_t AliAnalysisTaskCreateNUE::AcceptMCTruthTrack(AliAODMCParticle *mtrk) {
  // Pt cut
  if(mtrk->Pt() < fMinPt) return kFALSE;
  if(mtrk->Pt() > fMaxPt) return kFALSE;

  if(TMath::Abs(mtrk->Eta()) > fEtaCut) return kFALSE;

  if (!(mtrk->IsPhysicalPrimary())) return kFALSE;
  if (mtrk->Charge() == 0) return kFALSE;
  return kTRUE;
}


Bool_t AliAnalysisTaskCreateNUE::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp)
{
  // Pt cut
  if (mtr->Pt() < fMinPt)
    return kFALSE;
  if (mtr->Pt() > fMaxPt)
    return kFALSE;

  // DCA cut
  if (ltrackXYZ && vtxp)
  {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0] - vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1] - vtxp[1];
    ltrackXYZ[2] = abs(ltrackXYZ[2] - vtxp[2]);
  }
  else
    return kFALSE; //DCA cut is a must for now

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1
    return fGFWSelection15o->AcceptTrack(mtr, ltrackXYZ, 0, kFALSE);
  } else {
    return fGFWSelection->AcceptTrack(mtr, ltrackXYZ, 0, kFALSE);
  }
}

