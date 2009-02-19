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
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDkink.h"
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

  // status bits histogram
  if(!(h = (TH2I*)gROOT->FindObject("hTRDstat"))){
    h = new TH2I("hTRDstat", "TRD status bits", 100, 0., 20., kNbits, .5, kNbits+.5);
    h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h->GetYaxis()->SetTitle("status bits");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fHistos->AddAt(h, kTRDstat);

  // results array
  TObjArray *res = new TObjArray();
  res->SetName("Results");
  fHistos->AddAt(res, kResults);
}

//____________________________________________________________________
void AliTRDcheckESD::Exec(Option_t *){
  //
  // Run the Analysis
  //
  if(!fESD){
    AliError("ESD event missing.");
    return;
  }

  // Get MC information if available
  AliStack * fStack = 0x0;
  if(HasMC() && !fMC){ 
    AliWarning("MC event missing");
    SetMC(kFALSE);
  } else {
    if(!(fStack = fMC->Stack())){
      AliWarning("MC stack missing");
      SetMC(kFALSE);
    }
  }
  Bool_t TRDin(0), TRDout(0), TRDpid(0);

  AliESDtrack *esdTrack = 0x0;
  for(Int_t itrk = 0; itrk < fESD->GetNumberOfTracks(); itrk++){
    TRDin=0;TRDout=0;TRDpid=0;
    esdTrack = fESD->GetTrack(itrk);

//     if(esdTrack->GetNcls(1)) nTPC++;
//     if(esdTrack->GetNcls(2)) nTRD++;

    // track status
    ULong_t status = esdTrack->GetStatus();

    // define TPC out tracks
    if(!Bool_t(status & AliESDtrack::kTPCout)) continue;
    if(esdTrack->GetKinkIndex(0) > 0) continue;

    // TRD PID
    Double_t p[AliPID::kSPECIES]; esdTrack->GetTRDpid(p);
    // pid quality
    esdTrack->GetTRDpidQuality();

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
      AliWarning(Form("MC particle missing. Label[ %d].", fLabel));
      continue;
    }

    AliTrackReference *ref = 0x0; 
    Int_t nRefs = mcParticle->GetNumberOfTrackReferences();
    if(!nRefs){
      AliWarning(Form("Track refs missing. Label[%d].", fLabel));
      continue;
    }
    Int_t iref = 0;
    while(iref<nRefs){
      ref = mcParticle->GetTrackReference(iref);
      if(ref->LocalX() > xTPC) break;
      ref=0x0; iref++;
    }

    // read TParticle
    TParticle *tParticle = mcParticle->Particle(); 
    //Int_t fPdg = tParticle->GetPdgCode();
    // reject secondaries
    //if(!tParticle->IsPrimary()) continue;

    if(ref){ 
      if(ref->LocalX() > xTOF){ // track skipping TRD fiducial volume
        ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
      } else {
        TRDin=1;
        if(esdTrack->GetNcls(2)) TRDout=1;
        if(esdTrack->GetTRDpidQuality()) TRDpid=1;
      }
    } else { // track stopped in TPC 
      ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
    }
    // get the MC pt !!
    Float_t pt = ref->Pt();

    TH2 *h = (TH2I*)fHistos->At(kTRDstat);
    if(status & AliESDtrack::kTPCout) h->Fill(pt, kTPCout);
    if(/*status & AliESDtrack::k*/TRDin) h->Fill(pt, kTRDin);
    if(/*status & AliESDtrack::k*/TRDout){ 
      ((TH1*)fHistos->At(kNCl))->Fill(esdTrack->GetNcls(2));
      h->Fill(pt, kTRDout);
    }
    if(status & AliESDtrack::kTRDpid) h->Fill(pt, kTRDpid);
  }  
  PostData(0, fHistos);
}


//____________________________________________________________________
void AliTRDcheckESD::Terminate(Option_t *)
{
  TH2I *h2 = (TH2I*)fHistos->At(kTRDstat);
  TAxis *ax = h2->GetXaxis();
  TObjArray *res = (TObjArray*)fHistos->At(kResults);
  
  TH1 *h1[2] = {0x0, 0x0};
  TGraphErrors *g = 0x0;
  res->Expand(3);
  Int_t n1 = 0, n2 = 0, ip=0;
  Double_t eff = 0.;
  
  // geometrical efficiency
  h1[0] = h2->ProjectionX("px0", kTPCout, kTPCout);
  h1[1] = h2->ProjectionX("px1", kTRDin, kTRDin);
  res->Add(g = new TGraphErrors());
  g->SetNameTitle("geom", "TRD geometrical efficiency (TRDin/TPCout)");
  for(Int_t ib=1; ib<=ax->GetNbins(); ib++){
    if(!(n1 = (Int_t)h1[0]->GetBinContent(ib))) continue;
    n2 = (Int_t)h1[1]->GetBinContent(ib);
    eff = n2/Float_t(n1);

    ip=g->GetN();
    g->SetPoint(ip, ax->GetBinCenter(ib), eff);
    g->SetPointError(ip, 0., n2 ? eff*TMath::Sqrt(1./n1+1./n2) : 0.);
  }

  // tracking efficiency
  h1[0] = h2->ProjectionX("px0", kTRDin, kTRDin);
  h1[1] = h2->ProjectionX("px1", kTRDout, kTRDout);
  res->Add(g = new TGraphErrors());
  g->SetNameTitle("tracking", "TRD tracking efficiency (TRDout/TRDin)");
  for(Int_t ib=1; ib<=ax->GetNbins(); ib++){
    if(!(n1 = (Int_t)h1[0]->GetBinContent(ib))) continue;
    n2 = (Int_t)h1[1]->GetBinContent(ib);
    eff = n2/Float_t(n1);

    ip=g->GetN();
    g->SetPoint(ip, ax->GetBinCenter(ib), eff);
    g->SetPointError(ip, 0., n2 ? eff*TMath::Sqrt(1./n1+1./n2) : 0.);
  }

  // PID efficiency
  h1[0] = h2->ProjectionX("px0", kTRDin, kTRDin);
  h1[1] = h2->ProjectionX("px1", kTRDpid, kTRDpid);
  res->Add(g = new TGraphErrors());
  g->SetNameTitle("PID", "TRD PID efficiency (TRDpid/TRDin)");
  for(Int_t ib=1; ib<=ax->GetNbins(); ib++){
    if(!(n1 = (Int_t)h1[0]->GetBinContent(ib))) continue;
    n2 = (Int_t)h1[1]->GetBinContent(ib);
    eff = n2/Float_t(n1);

    ip=g->GetN();
    g->SetPoint(ip, ax->GetBinCenter(ib), eff);
    g->SetPointError(ip, 0., n2 ? eff*TMath::Sqrt(1./n1+1./n2) : 0.);
  }
}

// if(TRDin && !TRDout){
//   printf("[%c] x[%2d]=%7.2f pt[%7.3f] id[%d] kink[%d] label[%d] in[%d] out[%d] pdg[%d]\n", tParticle->IsPrimary() ? 'P' : 'S', iref, ref->LocalX(), pt, esdTrack->GetID(), esdTrack->GetKinkIndex(0), fLabel, TRDin, TRDout, tParticle->GetPdgCode());
//   printf("\tstatus ITS[%d %d] TPC[%d %d %d] \n"
//     ,Bool_t(status & AliESDtrack::kITSin)
//     ,Bool_t(status & AliESDtrack::kITSout)
//     ,Bool_t(status & AliESDtrack::kTPCin)
//     ,Bool_t(status & AliESDtrack::kTPCout)
//     ,Bool_t(status & AliESDtrack::kTPCrefit)
//   );
//   printf("clusters ITS[%d] TPC[%d] TRD[%d]\n", esdTrack->GetNcls(0), esdTrack->GetNcls(1), esdTrack->GetNcls(2));
//   if(esdTrack->GetKinkIndex(0) != 0){
//     AliESDkink *kink = fESD->GetKink(TMath::Abs(esdTrack->GetKinkIndex(0)));
//     if(kink) printf("index[%d %d] label[%d %d]\n", kink->GetIndex(0), kink->GetIndex(1), kink->GetLabel(0), kink->GetLabel(1));
//   }
// }
