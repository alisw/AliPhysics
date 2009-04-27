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

#include "AliTRDcheckESD.h"

ClassImp(AliTRDcheckESD)

const Int_t   AliTRDcheckESD::fgkNgraphs = 4;
const Float_t AliTRDcheckESD::fgkxTPC = 290.;
const Float_t AliTRDcheckESD::fgkxTOF = 365.;

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD():
  AliAnalysisTask("checkESD", "ESD checker for TRD info")
  ,fStatus(0)
  ,fESD(0x0)
  ,fMC(0x0)
  ,fHistos(0x0)
{
  //
  // Default constructor
  //
  SetMC(kTRUE);
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
  Histos();
}

//____________________________________________________________________
TGraphErrors* AliTRDcheckESD::GetGraph(Int_t id, Option_t*)
{
  Bool_t kBUILD = 1, // build graph if none found
         kCLEAR = 1; // clear existing graph

  const Char_t *name[] = {
    "Geo", "Trk", "Pid", "Ref"
  };
  const Char_t *title[] = {
    "TRD geometrical efficiency (TRDin/TPCout)"
    ,"TRD tracking efficiency (TRDout/TRDin)"
    ,"TRD PID efficiency (TRDpid/TRDin)"
    ,"TRD refit efficiency (TRDrefit/TRDin)"
  };
  const Int_t ngr = sizeof(name)/sizeof(Char_t*);
  if(ngr != fgkNgraphs){
    AliWarning("No of graphs defined different from definition");
    return 0x0;
  }

  TObjArray *res = 0x0;
  if(!(res = (TObjArray*)fHistos->At(kResults)) ||
      (id < 0 || id >= ngr)){
    AliWarning("Graph missing.");
    return 0x0;
  }

  TGraphErrors *g = 0x0;
  if((g = dynamic_cast<TGraphErrors*>(res->At(id)))){
    if(kCLEAR) for(Int_t ip=g->GetN(); ip--;) g->RemovePoint(ip);
  } else {
    if(kBUILD){
      g = new TGraphErrors();
      g->SetNameTitle(name[id], title[id]);
      res->AddAt(g, id);
    }
  }
  return g;
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
  if(HasMC()){
    if(!fMC){ 
      AliWarning("MC event missing");
      SetMC(kFALSE);
    } else {
      if(!(fStack = fMC->Stack())){
        AliWarning("MC stack missing");
        SetMC(kFALSE);
      }
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
    //PrintStatus(status);

    // define TPC out tracks
    if(!Bool_t(status & AliESDtrack::kTPCout)) continue;
    if(esdTrack->GetKinkIndex(0) > 0) continue;

    // TRD PID
    Double_t p[AliPID::kSPECIES]; esdTrack->GetTRDpid(p);
    // pid quality
    //esdTrack->GetTRDntrackletsPID();

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
      if(ref->LocalX() > fgkxTPC) break;
      ref=0x0; iref++;
    }

    // read TParticle
    //TParticle *tParticle = mcParticle->Particle(); 
    //Int_t fPdg = tParticle->GetPdgCode();
    // reject secondaries
    //if(!tParticle->IsPrimary()) continue;

    if(ref){ 
      if(ref->LocalX() > fgkxTOF){ // track skipping TRD fiducial volume
        ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
      } else {
        TRDin=1;
        if(esdTrack->GetNcls(2)) TRDout=1;
        if(esdTrack->GetTRDntrackletsPID()>=4) TRDpid=1;
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
    if(/*status & AliESDtrack::k*/TRDpid) h->Fill(pt, kTRDpid);
    if(status & AliESDtrack::kTRDrefit) h->Fill(pt, kTRDref);
  }  
  PostData(0, fHistos);
}

//____________________________________________________________________
TObjArray* AliTRDcheckESD::Histos()
{
  if(fHistos) return fHistos;

  fHistos = new TObjArray(kNhistos);
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
  return fHistos;
}

//____________________________________________________________________
Bool_t AliTRDcheckESD::Load(const Char_t *filename, const Char_t *name)
{
  if(!TFile::Open(filename)){
    AliWarning(Form("Couldn't open file %s.", filename));
    return kFALSE;
  }
  TObjArray *o = 0x0;
  if(!(o = (TObjArray*)gFile->Get(name ? name : GetName()))){
    AliWarning("Missing histogram container.");
    return kFALSE;
  }
  fHistos = (TObjArray*)o->Clone(GetName());
  gFile->Close();
  return kTRUE;
}

//____________________________________________________________________
void AliTRDcheckESD::Terminate(Option_t *)
{
  TObjArray *res = 0x0;
  if(!(res = (TObjArray*)fHistos->At(kResults))){
    AliWarning("Graph container missing.");
    return;
  }
  if(!res->GetEntriesFast()) res->Expand(fgkNgraphs);
  
  // geometrical efficiency
  TH2I *h2 = (TH2I*)fHistos->At(kTRDstat);
  TH1 *h1[2] = {0x0, 0x0};
  h1[0] = h2->ProjectionX("px0", kTPCout, kTPCout);
  h1[1] = h2->ProjectionX("px1", kTRDin, kTRDin);
  Process(h1, GetGraph(0));

  // tracking efficiency
  h1[0] = h2->ProjectionX("px0", kTRDin, kTRDin);
  h1[1] = h2->ProjectionX("px1", kTRDout, kTRDout);
  Process(h1, GetGraph(1));

  // PID efficiency
  //h1[0] = h2->ProjectionX("px0", kTRDin, kTRDin);
  h1[1] = h2->ProjectionX("px1", kTRDpid, kTRDpid);
  Process(h1, GetGraph(2));

  // Refit efficiency
  //h1[0] = h2->ProjectionX("px0", kTRDin, kTRDin);
  h1[1] = h2->ProjectionX("px1", kTRDref, kTRDref);
  Process(h1, GetGraph(3));
}

//____________________________________________________________________
void AliTRDcheckESD::Process(TH1 **h1, TGraphErrors *g)
{
  Int_t n1 = 0, n2 = 0, ip=0;
  Double_t eff = 0.;

  TAxis *ax = h1[0]->GetXaxis();
  for(Int_t ib=1; ib<=ax->GetNbins(); ib++){
    if(!(n1 = (Int_t)h1[0]->GetBinContent(ib))) continue;
    n2 = (Int_t)h1[1]->GetBinContent(ib);
    eff = n2/Float_t(n1);

    ip=g->GetN();
    g->SetPoint(ip, ax->GetBinCenter(ib), eff);
    g->SetPointError(ip, 0., n2 ? eff*TMath::Sqrt(1./n1+1./n2) : 0.);
  }
}  

//____________________________________________________________________
void AliTRDcheckESD::PrintStatus(ULong_t status)
{
  printf("ITS[i(%d) o(%d) r(%d)] TPC[i(%d) o(%d) r(%d) p(%d)] TRD[i(%d) o(%d) r(%d) p(%d) s(%d)] HMPID[o(%d) p(%d)]\n"
    ,Bool_t(status & AliESDtrack::kITSin)
    ,Bool_t(status & AliESDtrack::kITSout)
    ,Bool_t(status & AliESDtrack::kITSrefit)
    ,Bool_t(status & AliESDtrack::kTPCin)
    ,Bool_t(status & AliESDtrack::kTPCout)
    ,Bool_t(status & AliESDtrack::kTPCrefit)
    ,Bool_t(status & AliESDtrack::kTPCpid)
    ,Bool_t(status & AliESDtrack::kTRDin)
    ,Bool_t(status & AliESDtrack::kTRDout)
    ,Bool_t(status & AliESDtrack::kTRDrefit)
    ,Bool_t(status & AliESDtrack::kTRDpid)
    ,Bool_t(status & AliESDtrack::kTRDStop)
    ,Bool_t(status & AliESDtrack::kHMPIDout)
    ,Bool_t(status & AliESDtrack::kHMPIDpid)
  );
}

