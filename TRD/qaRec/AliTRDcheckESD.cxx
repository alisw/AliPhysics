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


#include <TClonesArray.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TH2I.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliTRDgeometry.h"

#include "AliTRDcheckESD.h"

ClassImp(AliTRDcheckESD)

const Float_t AliTRDcheckESD::xTPC = 290.;
const Float_t AliTRDcheckESD::xTOF = 365.;

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD():
  AliAnalysisTask("ESDchecker", "TRD checker @ ESD level")
  ,fStatus(0)
  ,fESD(0x0)
  ,fMC(0x0)
  ,fHistos(0x0)
{
  //
  // Default constructor
  //

  DefineInput(0, TChain::Class());
  DefineOutput(0, TObjArray::Class());
}

//____________________________________________________________________
AliTRDcheckESD::~AliTRDcheckESD()
{
  if(fHistos){
    //fHistos->Delete();
    delete fHistos;
  }
}

//____________________________________________________________________
void AliTRDcheckESD::ConnectInputData(Option_t *)
{
  //
  // Link the Input Data
  //
  TTree *tree = dynamic_cast<TChain*>(GetInputData(0));
  if(tree) tree->SetBranchStatus("Tracks", 1);

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fESD = esdH ? esdH->GetEvent() : 0x0;

  if(!HasMC()) return;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  fMC = mcH ? mcH->MCEvent() : 0x0;
}

//____________________________________________________________________
void AliTRDcheckESD::CreateOutputObjects()
{	
  //
  // Create Output Containers (TObjectArray containing 1D histograms)
  //
  OpenFile(0, "RECREATE");  
  fHistos = new TObjArray(5);
  //fHistos->SetOwner(kTRUE);
  
  TH1 *h = 0x0;

  // clusters per tracklet
  if(!(h = (TH1I*)gROOT->FindObject("hNCl"))){
    h = new TH1I("hNCl", "Clusters per TRD track", 100, 0., 200.);
    h->GetXaxis()->SetTitle("N_{cl}^{TRD}");
    h->GetYaxis()->SetTitle("entries");
  } else h->Reset();
  fHistos->AddAt(h, kNCl);

  // TPC out
  if(!(h = (TH2I*)gROOT->FindObject("hTRDstat"))){
    h = new TH2I("hTRDstat", "TRD status bits", 100, 0., 20., 4, -.5, 3.5);
    h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h->GetYaxis()->SetTitle("status bits");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fHistos->AddAt(h, kTRDstat);
}

//____________________________________________________________________
void AliTRDcheckESD::Exec(Option_t *){
  //
  // Run the Analysis
  //
  if(!fESD){
    AliError("ESD not found");
    return;
  }

  // Get MC information if available
  AliStack * fStack = 0x0;
  if(HasMC() && !fMC){ 
    AliWarning("Monte Carlo Event not available");
    SetMC(kFALSE);
  } else {
    if(!(fStack = fMC->Stack())){
      AliWarning("Cannot get the Monte Carlo Stack");
      SetMC(kFALSE);
    }
  }
  Bool_t TPCout(0), TRDin(0), TRDout(0), TRDpid(0);


  Int_t nTRD = 0, nTPC = 0;
  //Int_t nTracks = fESD->GetNumberOfTracks();
  AliESDtrack *esdTrack = 0x0;
  for(Int_t itrk = 0; itrk < fESD->GetNumberOfTracks(); itrk++){
    TPCout=0;TRDin=0;TRDout=0;TRDpid=0;
    esdTrack = fESD->GetTrack(itrk);
    if(esdTrack->GetNcls(1)) nTPC++;
    if(esdTrack->GetNcls(2)) nTRD++;

    // track status
    ULong_t status = esdTrack->GetStatus();

    // TRD PID
    Double_t p[AliPID::kSPECIES]; esdTrack->GetTRDpid(p);
    // pid quality
    esdTrack->GetTRDpidQuality();
    // kink index
    esdTrack->GetKinkIndex(0);
    // TPC clusters
    esdTrack->GetNcls(1);

    // look at external track param
    const AliExternalTrackParam *op = esdTrack->GetOuterParam();
    Double_t xyz[3];
    if(op){
      op->GetXYZ(xyz);
      op->Global2LocalPosition(xyz, op->GetAlpha());
      //printf("op @ X[%7.3f]\n", xyz[0]);
    }

    // read MC info
    if(!HasMC()) continue;

    Int_t fLabel = esdTrack->GetLabel();
    if(TMath::Abs(fLabel) > fStack->GetNtrack()) continue; 
    
    // read MC particle
    AliMCParticle *mcParticle = 0x0; 
    if(!(mcParticle = fMC->GetTrack(TMath::Abs(fLabel)))){
      AliWarning(Form("MC particle missing for ESD fLabel %d.", fLabel));
      continue;
    }

    AliTrackReference *ref = 0x0; 
    Int_t nRefs = mcParticle->GetNumberOfTrackReferences();
    Int_t iref = 0;
    while(iref<nRefs){
      ref = mcParticle->GetTrackReference(iref);
      if(ref->LocalX() > xTPC) break;
      ref=0x0; iref++;
    }

    // read TParticle
    TParticle *tParticle = mcParticle->Particle(); 
    Int_t fPdg = tParticle->GetPdgCode();
    //tParticle->IsPrimary();

    //printf("[%c] ref[%2d]=", tParticle->IsPrimary() ? 'P' : 'S', iref);

    TPCout=1;
    if(ref){
      if(ref->LocalX() > xTOF){ 
        //printf("  TOF   [");
        ref = mcParticle->GetTrackReference(iref-1);
      } else {
        //printf("%7.2f [", ref->LocalX());
        TRDin=1;
        if(esdTrack->GetNcls(2)) TRDout=1;
        if(esdTrack->GetTRDpidQuality()) TRDpid=1;
      }
    } else { 
      //printf("  TPC   [");
      ref = mcParticle->GetTrackReference(iref-1);
    }
    Float_t pt = ref->Pt();
    //printf("%f]\n", pt);
    
    TH2 *h = (TH2I*)fHistos->At(kTRDstat);
    if(/*status & AliESDtrack::k*/TPCout) h->Fill(pt, 0);
    if(/*status & AliESDtrack::k*/TRDin) h->Fill(pt, 1);
    if(/*status & AliESDtrack::k*/TRDout){ 
      ((TH1*)fHistos->At(kNCl))->Fill(esdTrack->GetNcls(2));
      h->Fill(pt, 2);
    }
    if(/*status & AliESDtrack::k*/TRDpid) h->Fill(pt, 3);
  }  

  PostData(0, fHistos);
}


//____________________________________________________________________
void AliTRDcheckESD::Terminate(Option_t *)
{
}
