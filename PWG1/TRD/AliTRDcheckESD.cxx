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

/////////////////////////////////////////////////////
//
// Check basic detector results at ESD level
//   - Geometrical efficiency  
//   - Tracking efficiency  
//   - PID efficiency  
//   - Refit efficiency  
//
// Author
//   Alex Bercuci <A.Bercuci@gsi.de>
//
//////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TPad.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH2I.h>
#include <TH3S.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
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

const Float_t AliTRDcheckESD::fgkxTPC = 290.;
const Float_t AliTRDcheckESD::fgkxTOF = 365.;
const UChar_t AliTRDcheckESD::fgkNgraph[AliTRDcheckESD::kNrefs] ={
8, 4, 2, 20};
FILE* AliTRDcheckESD::fgFile = NULL;

const Float_t AliTRDcheckESD::fgkEvVertexZ = 15.;
const Int_t   AliTRDcheckESD::fgkEvVertexN = 1;
const Float_t AliTRDcheckESD::fgkTrkDCAxy  = 40.;
const Float_t AliTRDcheckESD::fgkTrkDCAz   = 15.;
const Int_t   AliTRDcheckESD::fgkNclTPC    = 100;
const Float_t AliTRDcheckESD::fgkPt        = 0.2;
const Float_t AliTRDcheckESD::fgkEta       = 0.9;

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD():
  AliAnalysisTaskSE()
  ,fStatus(0)
  ,fNRefFigures(0)
  ,fESD(NULL)
  ,fMC(NULL)
  ,fHistos(NULL)
  ,fResults(NULL)
{
  //
  // Default constructor
  //
  SetNameTitle("checkESD", "Check TRD @ ESD level");
  SetMC(kTRUE);
}

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD(char* name):
  AliAnalysisTaskSE(name)
  ,fStatus(0)
  ,fNRefFigures(0)
  ,fESD(NULL)
  ,fMC(NULL)
  ,fHistos(NULL)
  ,fResults(NULL)
{
  //
  // Default constructor
  //
  SetMC(kTRUE);
  SetTitle("Check TRD @ ESD level");
  DefineOutput(1, TObjArray::Class());
}

//____________________________________________________________________
AliTRDcheckESD::~AliTRDcheckESD()
{
// Destructor
  if(fHistos){
    //fHistos->Delete();
    delete fHistos;
  }
  if(fResults){
    fResults->Delete();
    delete fResults;
  }
}

//____________________________________________________________________
void AliTRDcheckESD::UserCreateOutputObjects()
{	
  //
  // Create Output Containers (TObjectArray containing 1D histograms)
  //
  Histos();
}


//____________________________________________________________________
Bool_t AliTRDcheckESD::GetRefFigure(Int_t ifig)
{
  if(ifig>=fNRefFigures){
    AliWarning(Form("Ref plot %d not available. Valid only up to %d", ifig, fNRefFigures));
    return kFALSE;
  }
  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  } else {
    gPad->SetLogx(0);gPad->SetLogy(0);
    gPad->SetMargin(0.125, 0.015, 0.1, 0.015);
  }

  const Char_t *title[20];
  TH1 *hF(NULL);
  if((hF=(TH1S*)gROOT->FindObject("hFcheckESD"))) delete hF;
  TLegend *leg(NULL);
  TList *l(NULL); TVirtualPad *pad(NULL);
  TGraphErrors *g(NULL);TGraphAsymmErrors *ga(NULL);
  TObjArray *arr(NULL);
  switch(ifig){
  case kNCl: // number of clusters/track
    if(!(arr = (TObjArray*)fResults->At(kNCl))) return kFALSE;

    leg = new TLegend(.83, .7, .99, .96);
    leg->SetHeader("Species");
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    for(Int_t ig(0); ig<fgkNgraph[kNCl]; ig++){
      if(!(g = (TGraphErrors*)arr->At(ig))) return kFALSE;
      if(!g->GetN()) continue;
      g->Draw(ig?"pc":"apc"); leg->AddEntry(g, g->GetTitle(), "pl");
      if(ig) continue;
      hF=g->GetHistogram();
      hF->SetXTitle("no of clusters");
      hF->SetYTitle("entries"); 
      hF->GetYaxis()->CenterTitle(1);
      hF->GetYaxis()->SetTitleOffset(1.2);
      hF->SetMinimum(5);
    }
    leg->Draw(); gPad->SetLogy();
    break;
  case kTRDstat: // Efficiency
    if(!(arr = (TObjArray*)fResults->At(kTRDstat))) return kFALSE;
    leg = new TLegend(.62, .77, .98, .98);
    leg->SetHeader("TRD Efficiency");
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    title[0] = "Geometrical (TRDin/TPCout)";
    title[1] = "Tracking (TRDout/TRDin)";
    title[2] = "PID (TRDpid/TRDin)";
    title[3] = "Refit (TRDrefit/TRDin)";
    hF = new TH1S("hFcheckESD", ";p [GeV/c];Efficiency", 10, 0.1, 10.);
    hF->SetMaximum(1.4);
    hF->GetXaxis()->SetMoreLogLabels();
    hF->GetYaxis()->CenterTitle(1);
    hF->Draw("p");
    for(Int_t ig(0); ig<fgkNgraph[kTRDstat]; ig++){
      if(!(g = (TGraphErrors*)arr->At(ig))) return kFALSE;
      g->Draw("pl"); leg->AddEntry(g, title[ig], "pl");
      //PutTrendValue(name[id], g->GetMean(2));
      //PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
    }
    leg->Draw(); gPad->SetLogx();
    break;
  case kTRDmom: // Energy loss
    if(!(arr = (TObjArray*)fResults->At(kTRDmom))) return kFALSE;
    leg = new TLegend(.65, .7, .95, .99);
    leg->SetHeader("Energy Loss");
    leg->SetBorderSize(1); leg->SetFillColor(0);
    title[0] = "Max & 90% quantile";
    title[1] = "Mean & 60% quantile";
    hF = new TH1S("hFcheckESD", ";layer;#Delta E", 6, -0.5, 5.5);
    hF->SetMaximum(1.3);hF->SetMinimum(-.3);
    hF->Draw("p");
    for(Int_t ig(0); ig<fgkNgraph[kTRDmom]; ig++){
      if(!(ga = (TGraphAsymmErrors*)arr->At(ig))) return kFALSE;
      ga->Draw("pl"); leg->AddEntry(ga, title[ig], "pl");
      //PutTrendValue(name[id], g->GetMean(2));
      //PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
    }
    leg->Draw();gPad->SetLogx(kFALSE);
    break;
  case kPtRes: // Pt resolution @ vertex
    if(!(arr = (TObjArray*)fResults->At(kPtRes))) return kFALSE;
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l=gPad->GetListOfPrimitives(); 
    pad = ((TVirtualPad*)l->At(0)); pad->cd(); pad->SetLogx();
    pad->SetMargin(0.1, 0.022, 0.1, 0.023);
    hF = new TH1S("hFcheckESD", "ITS+TPC+TRD;p_{t} [GeV/c];#Delta p_{t} / p_{t} [%]", 10, 0.2, 10.);
    hF->SetMaximum(10.);hF->SetMinimum(-3.);
    hF->GetXaxis()->SetMoreLogLabels();
    hF->GetXaxis()->SetTitleOffset(1.2);
    hF->GetYaxis()->CenterTitle();
    hF->Draw("p");
    //for(Int_t ig(0); ig<fgkNgraph[kPtRes]/2; ig++){
    for(Int_t ig(2); ig<6; ig++){
      if(!(g = (TGraphErrors*)arr->At(ig))) continue;
      if(!g->GetN()) continue;
      g->Draw("pl");
      //PutTrendValue(name[id], g->GetMean(2));
      //PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
    }
    pad = ((TVirtualPad*)l->At(1)); pad->cd(); pad->SetLogx();
    pad->SetMargin(0.1, 0.22, 0.1, 0.023);
    hF = (TH1*)hF->Clone("hFcheckESD1");
    hF->SetTitle("ITS+TPC");
    hF->SetMaximum(10.);hF->SetMinimum(-3.);
    hF->Draw("p");
    leg = new TLegend(.78, .1, .99, .98);
    leg->SetHeader("P_{t} @ DCA");
    leg->SetBorderSize(1); leg->SetFillColor(0);
    leg->SetTextAlign(22);
    leg->SetTextFont(12);
    leg->SetTextSize(0.03813559);
    {
      Int_t nPlots(0);
      //for(Int_t ig(fgkNgraph[kPtRes]/2); ig<fgkNgraph[kPtRes]; ig++){
      for(Int_t ig(12); ig<16; ig++){
        if(!(g = (TGraphErrors*)arr->At(ig))) continue;
        if(!g->GetN()) continue;
        nPlots++;
        g->Draw("pl"); leg->AddEntry(g, g->GetTitle(), "pl");
        //PutTrendValue(name[id], g->GetMean(2));
        //PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
      }
      if(nPlots) leg->Draw();
    }

    break;
  }
  return kTRUE;
}

//____________________________________________________________________
void AliTRDcheckESD::UserExec(Option_t *){
  //
  // Run the Analysis
  //
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fMC = MCEvent();

  if(!fESD){
    AliError("ESD event missing.");
    return;
  }
  
  // Get MC information if available
  AliStack * fStack = NULL;
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
  TH2 *h(NULL);
  
  AliESDtrack *esdTrack(NULL);
  for(Int_t itrk = 0; itrk < fESD->GetNumberOfTracks(); itrk++){
    esdTrack = fESD->GetTrack(itrk);

    // track status
    ULong_t status = esdTrack->GetStatus(); //PrintStatus(status);

    // track selection
    Bool_t selected(kTRUE);
    if(esdTrack->Pt() < fgkPt){ 
      AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] Pt[%5.2f]", fESD->GetEventNumberInFile(), itrk, esdTrack->Pt()));
      selected = kFALSE;
    }
    if(TMath::Abs(esdTrack->Eta()) > fgkEta){
      AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] Eta[%5.2f]", fESD->GetEventNumberInFile(), itrk, TMath::Abs(esdTrack->Eta())));
      selected = kFALSE;
    }
    if(!Bool_t(status & AliESDtrack::kTPCout)){
      AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] !TPCout", fESD->GetEventNumberInFile(), itrk));
      selected = kFALSE;
    }
    if(esdTrack->GetKinkIndex(0) > 0){
      AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] Kink", fESD->GetEventNumberInFile(), itrk));
      selected = kFALSE;
    }
    if(esdTrack->GetTPCNcls() < fgkNclTPC){ 
      AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] NclTPC[%d]", fESD->GetEventNumberInFile(), itrk, esdTrack->GetTPCNcls()));
      selected = kFALSE;
    }
    Float_t par[2], cov[3];
    esdTrack->GetImpactParameters(par, cov);
    if(IsCollision()){ // cuts on DCA
      if(TMath::Abs(par[0]) > fgkTrkDCAxy){ 
        AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] DCAxy[%f]", fESD->GetEventNumberInFile(), itrk, TMath::Abs(par[0])));
        selected = kFALSE;
      }
      if(TMath::Abs(par[1]) > fgkTrkDCAz){ 
        AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] DCAz[%f]", fESD->GetEventNumberInFile(), itrk, TMath::Abs(par[1])));
        selected = kFALSE;
      }
    }
    if(!selected) continue;

    //Int_t nTPC(esdTrack->GetNcls(1));
    Int_t nTRD(esdTrack->GetNcls(2));
    Double_t pt(esdTrack->Pt());
    //Double_t eta(esdTrack->Eta());
    //Double_t phi(esdTrack->Phi());
    Double_t p[AliPID::kSPECIES]; esdTrack->GetTRDpid(p);
    // pid quality
    //esdTrack->GetTRDntrackletsPID();
    Bool_t kBarrel = Bool_t(status & AliESDtrack::kTRDin);

    // look at external track param
    const AliExternalTrackParam *op = esdTrack->GetOuterParam();
    const AliExternalTrackParam *ip = esdTrack->GetInnerParam();

    Double_t pt0(0.), eta0(0.), phi0(0.), ptTRD(0.); 
    // read MC info if available
    Bool_t kFOUND(kFALSE), kPhysPrim(kFALSE);
    AliMCParticle *mcParticle(NULL);
    if(HasMC()){
      AliTrackReference *ref(NULL); 
      Int_t fLabel(esdTrack->GetLabel());
      Int_t fIdx(TMath::Abs(fLabel));
      if(fIdx > fStack->GetNtrack()) continue; 
      
      // read MC particle 
      if(!(mcParticle = (AliMCParticle*) fMC->GetTrack(fIdx))) {
        AliWarning(Form("MC particle missing. Label[ %d].", fLabel));
        continue;
      }
      pt0  = mcParticle->Pt();
      eta0 = mcParticle->Eta();
      phi0 = mcParticle->Phi();
      kPhysPrim = fMC->IsPhysicalPrimary(fIdx);

      // read track references
      Int_t nRefs = mcParticle->GetNumberOfTrackReferences();
      if(!nRefs){
        AliWarning(Form("No TR found for track @ Label[%d].", fLabel));
        continue;
      }
      Int_t iref = 0;
      while(iref<nRefs){
        ref = mcParticle->GetTrackReference(iref);
        if(ref->LocalX() > fgkxTPC) break;
        ref=NULL; iref++;
      }
      if(ref){ 
        if(ref->LocalX() > fgkxTOF){ // track skipping TRD fiducial volume
          ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
        }
      } else { // track stopped in TPC 
        ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
      }
      ptTRD = ref->Pt();kFOUND=kTRUE;
    } else { // use reconstructed values
      if(op){
        Double_t x(op->GetX());
        if(x<fgkxTOF && x>fgkxTPC){
          ptTRD=op->Pt();
          kFOUND=kTRUE;
        }
      }

      if(!kFOUND && ip){
        ptTRD=ip->Pt();
        kFOUND=kTRUE;
      }
    }

    if(kFOUND){
      h = (TH2I*)fHistos->At(kTRDstat);
      if(status & AliESDtrack::kTPCout) h->Fill(ptTRD, kTPCout);
      if(status & AliESDtrack::kTRDin) h->Fill(ptTRD, kTRDin);
      if(kBarrel && (status & AliESDtrack::kTRDout)) h->Fill(ptTRD, kTRDout);
      if(kBarrel && (status & AliESDtrack::kTRDpid)) h->Fill(ptTRD, kTRDpid);
      if(kBarrel && (status & AliESDtrack::kTRDrefit)) h->Fill(ptTRD, kTRDref);
    }
    Int_t idx(HasMC() ? Pdg2Idx(TMath::Abs(mcParticle->PdgCode())): 0)
         ,sgn(esdTrack->Charge()<0?0:1);
    if(kBarrel && kPhysPrim) {
      TH3 *h3 = (TH3S*)fHistos->At(kPtRes);
      Int_t offset = (status & AliESDtrack::kTRDrefit) ? 0 : 10; 
      h3->Fill(pt0, 1.e2*(pt/pt0-1.), 
        offset + 2*idx + sgn);
    }
    ((TH1*)fHistos->At(kNCl))->Fill(nTRD, 2*idx + sgn);
    if(ip){
      h = (TH2I*)fHistos->At(kTRDmom);
      Float_t pTRD(0.);
      for(Int_t ily=6; ily--;){
        if((pTRD=esdTrack->GetTRDmomentum(ily))<0.) continue;
        h->Fill(ip->GetP()-pTRD, ily);
      }
    }
  }  
  PostData(1, fHistos);
}

//____________________________________________________________________
TObjArray* AliTRDcheckESD::Histos()
{
// Retrieve histograms array if already build or build it

  if(fHistos) return fHistos;

  fHistos = new TObjArray(kNhistos);
  //fHistos->SetOwner(kTRUE);

  TH1 *h = NULL;

  // clusters per track
  const Int_t kNpt(30);
  Float_t Pt(0.2);
  Float_t binsPt[kNpt+1];
  for(Int_t i=0;i<kNpt+1; i++,Pt+=(TMath::Exp(i*i*.001)-1.)) binsPt[i]=Pt;
  if(!(h = (TH2S*)gROOT->FindObject("hNCl"))){
    h = new TH2S("hNCl", "Clusters per TRD track;N_{cl}^{TRD};SPECIES;entries", 60, 0., 180., 10, -0.5, 9.5);
    TAxis *ay(h->GetYaxis());
    ay->SetLabelOffset(0.015);
    for(Int_t i(0); i<AliPID::kSPECIES; i++){
      ay->SetBinLabel(2*i+1, Form("%s^{-}", AliPID::ParticleLatexName(i)));
      ay->SetBinLabel(2*i+2, Form("%s^{+}", AliPID::ParticleLatexName(i)));
    }
  } else h->Reset();
  fHistos->AddAt(h, kNCl); fNRefFigures++;

  // status bits histogram
  const Int_t kNbits(5);
  Float_t Bits(.5);
  Float_t binsBits[kNbits+1];
  for(Int_t i=0; i<kNbits+1; i++,Bits+=1.) binsBits[i]=Bits;
  if(!(h = (TH2I*)gROOT->FindObject("hTRDstat"))){
    h = new TH2I("hTRDstat", "TRD status bits;p_{t} @ TRD [GeV/c];status;entries", kNpt, binsPt, kNbits, binsBits);
    TAxis *ay(h->GetYaxis());
    ay->SetBinLabel(1, "kTPCout");
    ay->SetBinLabel(2, "kTRDin");
    ay->SetBinLabel(3, "kTRDout");
    ay->SetBinLabel(4, "kTRDpid");
    ay->SetBinLabel(5, "kTRDrefit");
  } else h->Reset();
  fHistos->AddAt(h, kTRDstat);

  // energy loss
  if(!(h = (TH2I*)gROOT->FindObject("hTRDmom"))){
    h = new TH2I("hTRDmom", "TRD energy loss;p_{inner} - p_{ly} [GeV/c];ly;entries", 100, -1., 2., 6, -0.5, 5.5);
  } else h->Reset();
  fHistos->AddAt(h, kTRDmom);
  if(!HasMC()) return fHistos;

  // pt resolution
  const Int_t kNdpt(100), kNspec(4*AliPID::kSPECIES);
  Float_t DPt(-3.), Spec(-0.5);
  Float_t binsDPt[kNdpt+1], binsSpec[kNspec+1];
  for(Int_t i=0; i<kNdpt+1; i++,DPt+=6.e-2) binsDPt[i]=DPt;
  for(Int_t i=0; i<kNspec+1; i++,Spec+=1.) binsSpec[i]=Spec;
  if(!(h = (TH3S*)gROOT->FindObject("hPtRes"))){
    h = new TH3S("hPtRes", "P_{t} resolution @ DCA;p_{t}^{MC} [GeV/c];#Delta p_{t}/p_{t}^{MC} [%];SPECIES", kNpt, binsPt, kNdpt, binsDPt, kNspec, binsSpec);
    TAxis *az(h->GetZaxis());
    az->SetLabelOffset(0.015);
    for(Int_t i(0); i<AliPID::kSPECIES; i++){
      az->SetBinLabel(2*i+1, Form("%s^{-}", AliPID::ParticleLatexName(i)));
      az->SetBinLabel(2*i+2, Form("%s^{+}", AliPID::ParticleLatexName(i)));
      az->SetBinLabel(10+2*i+1, Form("%s^{-}", AliPID::ParticleLatexName(i)));
      az->SetBinLabel(10+2*i+2, Form("%s^{+}", AliPID::ParticleLatexName(i)));
    }
  } else h->Reset();
  fHistos->AddAt(h, kPtRes);

  return fHistos;
}

//____________________________________________________________________
Bool_t AliTRDcheckESD::Load(const Char_t *file, const Char_t *dir, const Char_t *name)
{
// Load data from performance file

  if(!TFile::Open(file)){
    AliWarning(Form("Couldn't open file %s.", file));
    return kFALSE;
  }
  if(dir){
    if(!gFile->cd(dir)){
      AliWarning(Form("Couldn't cd to %s in %s.", dir, file));
      return kFALSE;
    }
  }
  TObjArray *o(NULL);
  const Char_t *tn=(name ? name : GetName());
  if(!(o = (TObjArray*)gDirectory->Get(tn))){
    AliWarning(Form("Missing histogram container %s.", tn));
    return kFALSE;
  }
  fHistos = (TObjArray*)o->Clone(GetName());
  gFile->Close();
  return kTRUE;
}

//_______________________________________________________
Bool_t AliTRDcheckESD::PutTrendValue(const Char_t *name, Double_t val)
{
// Dump trending value to default file

  if(!fgFile){
    fgFile = fopen("TRD.Performance.txt", "at");
  }
  fprintf(fgFile, "%s_%s %f\n", GetName(), name, val);
  return kTRUE;
}

//____________________________________________________________________
void AliTRDcheckESD::Terminate(Option_t *)
{
// Steer post-processing 
  if(!fHistos){
    fHistos = dynamic_cast<TObjArray *>(GetOutputData(1));
    if(!fHistos){
      AliError("Histogram container not found in output");
      return;
    }
  }

  const Char_t *name[kNrefs] = {
    "Ncl", "Eff", "Eloss", "PtResDCA"
  };
  TObjArray *arr(NULL); TGraph *g(NULL);
  if(!fResults){
    fResults = new TObjArray(kNrefs);
    fResults->SetOwner();
    fResults->SetName("results");
    for(Int_t iref(0); iref<kNrefs; iref++){
      fResults->AddAt(arr = new TObjArray(fgkNgraph[iref]), iref);
      arr->SetName(name[iref]);  arr->SetOwner();
      switch(iref){
      case kNCl:
        for(Int_t ig(0); ig<fgkNgraph[iref]; ig++){
          arr->AddAt(g = new TGraphErrors(), ig);
          g->SetLineColor(ig+1); 
          g->SetMarkerColor(ig+1); 
          g->SetMarkerStyle(ig+20); 
          g->SetName(Form("s%d", ig));
          switch(ig){
          case 0: g->SetTitle("ALL"); break;
          case 1: g->SetTitle("NEG"); break;
          case 2: g->SetTitle("POS"); break;
          default: g->SetTitle(AliPID::ParticleLatexName(ig-3)); break;
          };
        }
        break;
      case kTRDmom:
        for(Int_t ig(0); ig<fgkNgraph[iref]; ig++){
          arr->AddAt(g = new TGraphAsymmErrors(), ig);
          g->SetLineColor(ig+1); 
          g->SetMarkerColor(ig+1); 
          g->SetMarkerStyle(ig+20); 
        }
        break;
      case kPtRes:
        for(Int_t idx(0); idx<AliPID::kSPECIES; idx++){
          Int_t ig(2*idx);
          arr->AddAt(g = new TGraphErrors(), ig);
          g->SetLineColor(kRed-idx); 
          g->SetMarkerColor(kRed-idx); 
          g->SetMarkerStyle(20+idx); 
          g->SetNameTitle(Form("s%d", ig), Form("res %s", AliPID::ParticleLatexName(idx)));
          arr->AddAt(g = new TGraphErrors(), ig+1);
          g->SetLineColor(kBlue-idx); 
          g->SetMarkerColor(kBlue-idx); 
          g->SetMarkerStyle(20+idx); 
          g->SetNameTitle(Form("m%d", ig+1), Form("sys %s", AliPID::ParticleLatexName(idx)));

          ig+=10;
          arr->AddAt(g = new TGraphErrors(), ig);
          g->SetLineColor(kRed-idx); 
          g->SetMarkerColor(kRed-idx); 
          g->SetMarkerStyle(20+idx); 
          g->SetNameTitle(Form("s%d", ig), Form("sigma %s", AliPID::ParticleLatexName(idx)));
          arr->AddAt(g = new TGraphErrors(), ig+1);
          g->SetLineColor(kBlue-idx); 
          g->SetMarkerColor(kBlue-idx); 
          g->SetMarkerStyle(20+idx); 
          g->SetNameTitle(Form("m%d", ig+1), Form("mean %s", AliPID::ParticleLatexName(idx)));
        }
        break;
      default:
        for(Int_t ig(0); ig<fgkNgraph[iref]; ig++){
          arr->AddAt(g = new TGraphErrors(), ig);
          g->SetLineColor(ig+1); 
          g->SetMarkerColor(ig+1); 
          g->SetMarkerStyle(ig+20); 
        }
        break;
      }
    }
  }
  TH1 *h1[2] = {NULL, NULL};
  TH2I *h2(NULL);
  TAxis *ax(NULL);

  // No of clusters
  if(!(h2 = (TH2I*)fHistos->At(kNCl))) return;
  ax = h2->GetXaxis();
  arr = (TObjArray*)fResults->At(kNCl);
  // All tracks
  h1[0] = h2->ProjectionX("Ncl_px");
  TGraphErrors *ge=(TGraphErrors*)arr->At(0);
  for(Int_t ib=2; ib<=ax->GetNbins(); ib++){
    ge->SetPoint(ib-2, ax->GetBinCenter(ib), h1[0]->GetBinContent(ib));
  }
  // All charged tracks
  TH1 *hNclCh[2] = {(TH1D*)h1[0]->Clone("NEG"), (TH1D*)h1[0]->Clone("POS")};
  hNclCh[0]->Reset();hNclCh[1]->Reset();
  for(Int_t is(1); is<=AliPID::kSPECIES; is++){
    hNclCh[0]->Add(h2->ProjectionX("Ncl_px", 2*is-1, 2*is-1)); // neg
    hNclCh[1]->Add(h2->ProjectionX("Ncl_px", 2*is, 2*is));     // pos
  }
  if(Int_t(hNclCh[0]->GetEntries())){
    ge=(TGraphErrors*)arr->At(1);
    for(Int_t ib=2; ib<=ax->GetNbins(); ib++){
      ge->SetPoint(ib-2, ax->GetBinCenter(ib), hNclCh[0]->GetBinContent(ib));
    }
  }
  if(Int_t(hNclCh[1]->GetEntries())){
    ge=(TGraphErrors*)arr->At(2);
    for(Int_t ib=2; ib<=ax->GetNbins(); ib++){
      ge->SetPoint(ib-2, ax->GetBinCenter(ib), hNclCh[1]->GetBinContent(ib));
    }
  }
  // Species wise
  for(Int_t is(1); is<=AliPID::kSPECIES; is++){
    h1[0] = h2->ProjectionX("Ncl_px", 2*is-1, 2*is);
    if(!Int_t(h1[0]->GetEntries())) continue;
    ge=(TGraphErrors*)arr->At(2+is);
    for(Int_t ib=2; ib<=ax->GetNbins(); ib++){
      ge->SetPoint(ib-2, ax->GetBinCenter(ib), h1[0]->GetBinContent(ib));
    }
  }
  fNRefFigures = 1;

  // EFFICIENCY
  // geometrical efficiency
  if(!(h2 = (TH2I*)fHistos->At(kTRDstat))) return;
  arr = (TObjArray*)fResults->At(kTRDstat);
  h1[0] = h2->ProjectionX("checkESDx0", kTPCout, kTPCout);
  h1[1] = h2->ProjectionX("checkESDx1", kTRDin, kTRDin);
  Process(h1, (TGraphErrors*)arr->At(0));
  delete h1[0];delete h1[1];
  // tracking efficiency
  h1[0] = h2->ProjectionX("checkESDx0", kTRDin, kTRDin);
  h1[1] = h2->ProjectionX("checkESDx1", kTRDout, kTRDout);
  Process(h1, (TGraphErrors*)arr->At(1));
  delete h1[1];
  // PID efficiency
  h1[1] = h2->ProjectionX("checkESDx1", kTRDpid, kTRDpid);
  Process(h1, (TGraphErrors*)arr->At(2));
  delete h1[1];
  // Refit efficiency
  h1[1] = h2->ProjectionX("checkESDx1", kTRDref, kTRDref);
  Process(h1, (TGraphErrors*)arr->At(3));
  delete h1[1];
  fNRefFigures++;

  // ENERGY LOSS
  if(!(h2 = dynamic_cast<TH2I*>(fHistos->At(kTRDmom)))) return;
  arr = (TObjArray*)fResults->At(kTRDmom);
  TGraphAsymmErrors *g06 = (TGraphAsymmErrors*)arr->At(0), *g09 = (TGraphAsymmErrors*)arr->At(1);
  ax=h2->GetXaxis();
  const Int_t nq(4);
  const Double_t xq[nq] = {0.05, 0.2, 0.8, 0.95};
  Double_t yq[nq];
  for(Int_t ily=6; ily--;){
    h1[0] = h2->ProjectionX("checkESDp0", ily+1, ily+1);
    h1[0]->GetQuantiles(nq,yq,xq);
    g06->SetPoint(ily, Float_t(ily), ax->GetBinCenter(h1[0]->GetMaximumBin()));
    g06->SetPointError(ily, 0., 0., TMath::Abs(yq[0]), yq[3]);
    g09->SetPoint(ily, Float_t(ily), h1[0]->GetMean());
    g09->SetPointError(ily, 0., 0., TMath::Abs(yq[1]), yq[2]);

    //printf(" max[%f] mean[%f] q[%f %f %f %f]\n", ax->GetBinCenter(h1[0]->GetMaximumBin()), h1[0]->GetMean(), yq[0], yq[1], yq[2], yq[3]);
    delete h1[0];
  }
  fNRefFigures++;
  if(!HasMC()) return;

  // Pt RESOLUTION @ DCA
  TH3S* h3(NULL); TGraphErrors *gg[2] = {NULL,NULL};
  if(!(h3 = dynamic_cast<TH3S*>(fHistos->At(kPtRes)))) return;
  arr = (TObjArray*)fResults->At(kPtRes);
  TAxis *az(h3->GetZaxis());
  for(Int_t i(0); i<AliPID::kSPECIES; i++){
    Int_t idx(2*i);
    az->SetRange(idx+1, idx+2); 
    gg[1] = (TGraphErrors*)arr->At(idx);
    gg[0] = (TGraphErrors*)arr->At(idx+1);
    Process2D((TH2*)h3->Project3D("yx"), gg);

    idx+=10;
    az->SetRange(idx+1, idx+2); 
    gg[1] = (TGraphErrors*)arr->At(idx);
    gg[0] = (TGraphErrors*)arr->At(idx+1);
    Process2D((TH2*)h3->Project3D("yx"), gg);
  }
  fNRefFigures++;
}

//____________________________________________________________________
Int_t AliTRDcheckESD::Pdg2Idx(Int_t pdg)
{
  switch(pdg){
  case kElectron: return AliPID::kElectron;  
  case kMuonMinus: return AliPID::kMuon;  
  case kPiPlus: return AliPID::kPion;  
  case kKPlus: return AliPID::kKaon;
  case kProton: return AliPID::kProton;
  } 
  return -1;
}

//____________________________________________________________________
void AliTRDcheckESD::Process(TH1 **h1, TGraphErrors *g)
{
// Generic function to process one reference plot

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
//________________________________________________________
void AliTRDcheckESD::Process2D(TH2 * const h2, TGraphErrors **g)
{
  //
  // Do the processing
  //

  Int_t n = 0;
  if((n=g[0]->GetN())) for(;n--;) g[0]->RemovePoint(n);
  if((n=g[1]->GetN())) for(;n--;) g[1]->RemovePoint(n);
  TF1 f("fg", "gaus", -3.,3.);
  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t x = h2->GetXaxis()->GetBinCenter(ibin);
    TH1D *h = h2->ProjectionY("py", ibin, ibin);
    if(h->GetEntries()<100) continue;
    //AdjustF1(h, f);

    h->Fit(&f, "QN");
    Int_t ip = g[0]->GetN();
    g[0]->SetPoint(ip, x, f.GetParameter(1));
    g[0]->SetPointError(ip, 0., f.GetParError(1));
    g[1]->SetPoint(ip, x, f.GetParameter(2));
    g[1]->SetPointError(ip, 0., f.GetParError(2));
  }
  return;
}
//____________________________________________________________________
void AliTRDcheckESD::PrintStatus(ULong_t status)
{
// Dump track status to stdout

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

