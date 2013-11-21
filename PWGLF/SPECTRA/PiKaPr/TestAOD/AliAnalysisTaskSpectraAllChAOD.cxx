/**************************************************************************
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

//-----------------------------------------------------------------
//         AliAnalysisTaskSpectraAllChAOD class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskSpectraAllChAOD.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODEventCuts.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include <TMCProcess.h>

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskSpectraAllChAOD)

//________________________________________________________________________
AliAnalysisTaskSpectraAllChAOD::AliAnalysisTaskSpectraAllChAOD(const char *name) : AliAnalysisTaskSE(name), fAOD(0), fTrackCuts(0), fEventCuts(0), fIsMC(0), fOutput(0)
{
  // Default constructor
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliSpectraAODEventCuts::Class());
  DefineOutput(3, AliSpectraAODTrackCuts::Class());
  
}
//________________________________________________________________________
//________________________________________________________________________
void AliAnalysisTaskSpectraAllChAOD::UserCreateOutputObjects()
{
  // create output objects
  fOutput=new TList();
  fOutput->SetOwner();
  fOutput->SetName("fOutput");
  
  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");
  
  //dimensions of standard THnSparse
  const Int_t nvar=4;
  const Double_t ptBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
  Int_t nptBins=52;
  
  Int_t    binsHistReal[nvar] = {   100,   100,      nptBins,     20};
  Double_t xminHistReal[nvar] = {    0.,     0,            0.,      0.};
  Double_t xmaxHistReal[nvar] = {  10.,   10.,            5.,   100.};
  
  THnSparseF* NSparseHist = new THnSparseF("NSparseHist","NSparseHist",nvar,binsHistReal,xminHistReal,xmaxHistReal);
  NSparseHist->GetAxis(0)->SetTitle("Q vec VZERO-C");
  NSparseHist->GetAxis(1)->SetTitle("Q vec VZERO-A");
  NSparseHist->GetAxis(2)->SetTitle("#it{p}_{T}");
  NSparseHist->SetBinEdges(2,ptBins);
  NSparseHist->GetAxis(3)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  fOutput->Add(NSparseHist);
  
  TH2F* hGen = new TH2F("hGen","hGen",nptBins,ptBins,20,0.,100.);
  hGen->GetXaxis()->SetTitle("#it{p}_{T,Gen}");
  hGen->GetYaxis()->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  fOutput->Add(hGen);
  
  TH2F* hRec = new TH2F("hRec","hRec",nptBins,ptBins,20,0.,100.);
  hRec->GetXaxis()->SetTitle("#it{p}_{T,Rec}");
  hRec->GetYaxis()->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  fOutput->Add(hRec);
  
  PostData(1, fOutput  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  
}
//________________________________________________________________________
void AliAnalysisTaskSpectraAllChAOD::UserExec(Option_t *)
{
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (!fAOD) {
    AliWarning("ERROR: AliAODEvent not available \n");
    return;
  }
  
  if (strcmp(fAOD->ClassName(), "AliAODEvent"))
    {
      AliFatal("Not processing AODs");
    }
  
  if(!fEventCuts->IsSelected(fAOD,fTrackCuts))return;//event selection
  
  // First do MC to fill up the MC particle array, such that we can use it later
  TClonesArray *arrayMC = 0;
  if (fIsMC)
    {
      arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!arrayMC) {
	AliFatal("Error: MC particles branch not found!\n");
      }
      Int_t nMC = arrayMC->GetEntries();
      for (Int_t iMC = 0; iMC < nMC; iMC++)
	{
	  AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
	  if(!partMC->Charge()) continue;//Skip neutrals
	  if(!partMC->IsPhysicalPrimary()) continue; //skip secondaries
	  if(partMC->Eta() < fTrackCuts->GetEtaMin() || partMC->Eta() > fTrackCuts->GetEtaMax())continue;//ETA CUT ON GENERATED!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
	  ((TH2F*)fOutput->FindObject("hGen"))->Fill(partMC->Pt(),fEventCuts->GetCent());
	  //Printf("a particle");
 
	}
    }
  
  //main loop on tracks
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* track = fAOD->GetTrack(iTracks);
    if (!fTrackCuts->IsSelected(track,kTRUE)) continue;
    
    ((TH2F*)fOutput->FindObject("hRec"))->Fill(track->Pt(),fEventCuts->GetCent());
    
    Double_t var[4];
    var[0]=fEventCuts->GetqV0C();
    var[1]=fEventCuts->GetqV0A();
    var[2]=track->Pt();
    var[3]=fEventCuts->GetCent();
    ((THnSparseF*)fOutput->FindObject("NSparseHist"))->Fill(var);
    
    //Printf("a track");
    
  } // end loop on tracks
  
  PostData(1, fOutput  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
}

//_________________________________________________________________
void   AliAnalysisTaskSpectraAllChAOD::Terminate(Option_t *)
{
  // Terminate
}
