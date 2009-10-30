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
#include <TObject.h>
#include <TH2I.h>
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
FILE* AliTRDcheckESD::fgFile = 0x0;

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD():
  AliAnalysisTask("checkESD", "ESD checker for TRD info")
  ,fStatus(0)
  ,fESD(0x0)
  ,fMC(0x0)
  ,fHistos(0x0)
  ,fResults(0x0)
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
Bool_t AliTRDcheckESD::GetRefFigure(Int_t ifig)
{
  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  }
  TGraph *g(0x0); TH1 *h(0x0); TLegend *leg(0x0);
  switch(ifig){
  case kNCl:
    if(!fHistos || !(h=(TH1I*)fHistos->At(kNCl))) break;
    h->Draw("c");
    return kTRUE;
  case kTRDstat:
    leg=new TLegend(.5, .75, .98, .98);
    leg->SetBorderSize(1); leg->SetFillColor(0);
    leg->SetHeader("Tracking Efficiency");
    for(Int_t ieff=0; ieff<4; ieff++){
      if(!(g=GetGraph(ieff, ""))) break;
      if(!ieff){
        if((h=(TH1S*)gROOT->FindObject("frame"))) delete h;
        h=new TH1S("h","",100, 0., 15.);
        h->SetLineColor(1);h->SetLineWidth(1);
        h->SetMinimum(0.2);h->SetMaximum(1.2);
        h->SetXTitle("p_{T} [GeV/c]");
        h->SetYTitle("Efficiency");
        h->Draw();
      }
      g->Draw("pc"); leg->AddEntry(g, g->GetTitle(), "pl");
    }
    leg->Draw();
    return kTRUE;
  case kTRDmom:
    leg=new TLegend(.6, .75, .98, .98);
    leg->SetBorderSize(1); leg->SetFillColor(0);
    leg->SetHeader("Energy loss");
    for(Int_t ieff=0; ieff<2; ieff++){
      if(!(g=GetGraph(4+ieff, ""))) break;
      if(!ieff){
        if((h=(TH1S*)gROOT->FindObject("frame"))) delete h;
        h=new TH1S("h","",100, -0.5, 5.5);
        h->SetLineColor(0);h->SetLineWidth(1);
        h->SetMinimum(-0.5);h->SetMaximum(1.2);
        h->SetXTitle("layer");
        h->SetYTitle("p_{T}^{Inner TPC} - p_{T}^{layer} [GeV/c]");
        h->Draw();
      }
      g->Draw("p"); leg->AddEntry(g, g->GetTitle(), "pl");
    }
    leg->Draw();
    return kTRUE;
  default:
    break;
  }
  AliInfo(Form("Reference plot [%d] missing result", ifig));
  return kFALSE;
}

//____________________________________________________________________
TGraph* AliTRDcheckESD::GetGraph(Int_t id, Option_t *opt)
{
// Retrieve graph with "id"
// Possible options are :
//   "b" - build graph if none found
//   "c" - clear existing graph

  Bool_t kBUILD = strstr(opt, "b"), // build graph if none found
         kCLEAR = strstr(opt, "c"); // clear existing graph

  const Char_t *name[] = {
    "Geo", "Trk", "Pid", "Ref", "Max06", "Mean09"
  };
  const Char_t *title[] = {
    "TRD geometrical efficiency (TRDin/TPCout)"
    ,"TRD tracking efficiency (TRDout/TRDin)"
    ,"TRD PID efficiency (TRDpid/TRDin)"
    ,"TRD refit efficiency (TRDrefit/TRDin)"
    ,"TRD Eloss (Max/90% quantile)"
    ,"TRD Eloss (Mean/60% quantile)"
  };
  const Int_t ngr = sizeof(name)/sizeof(Char_t*);
  if(ngr != kNgraphs){
    AliWarning("No of graphs defined different from definition");
    return 0x0;
  }

  if(!fResults){
    fResults = new TObjArray(kNgraphs);
    fResults->SetOwner();
    fResults->SetName("results");
  }

  TGraph *g = 0x0;
  if((g = dynamic_cast<TGraph*>(fResults->At(id)))){
    if(kCLEAR){ 
      for(Int_t ip=g->GetN(); ip--;) g->RemovePoint(ip);
    } else {
      PutTrendValue(name[id], g->GetMean(2));
      PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
    }
  } else {
    if(kBUILD){
      switch(id){
      case 0:
        g = new TGraphErrors();
        g->SetMarkerStyle(7);g->SetMarkerColor(kBlack);
        g->SetLineColor(kBlack);
        break;
      case 1:
        g = new TGraphErrors();
        g->SetMarkerStyle(7);g->SetMarkerColor(kRed);
        g->SetLineColor(kRed);
        break;
      case 2:
        g = new TGraphErrors();
        g->SetMarkerStyle(7);g->SetMarkerColor(kBlue);
        g->SetLineColor(kBlue);
        break;
      case 3:
        g = new TGraphErrors();
        g->SetMarkerStyle(7);g->SetMarkerColor(kGreen);
        g->SetLineColor(kGreen);
        break;
      case 4:
        g = new TGraphAsymmErrors(6);
        g->SetMarkerStyle(22);g->SetMarkerColor(kRed);
        g->SetLineColor(kBlack);g->SetLineWidth(2);
        break;
      case 5:
        g = new TGraphAsymmErrors(6);
        g->SetMarkerStyle(21);
        g->SetLineColor(kRed);g->SetLineWidth(2);
        break;
      default:
        AliWarning(Form("Graph index[%d] missing/not defined.", id));
        return 0x0;
      }
      g->SetNameTitle(name[id], title[id]);
      fResults->AddAt(g, id);
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
  Bool_t bTRDin(0), bTRDout(0), bTRDpid(0);

  AliESDtrack *esdTrack = 0x0;
  for(Int_t itrk = 0; itrk < fESD->GetNumberOfTracks(); itrk++){
    bTRDin=0;bTRDout=0;bTRDpid=0;
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
    const AliExternalTrackParam *ip = esdTrack->GetInnerParam();
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
    if(!(mcParticle = (AliMCParticle*) fMC->GetTrack(TMath::Abs(fLabel)))){
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
        bTRDin=1;
        if(esdTrack->GetNcls(2)) bTRDout=1;
        if(esdTrack->GetTRDntrackletsPID()>=4) bTRDpid=1;
      }
    } else { // track stopped in TPC 
      ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
    }
    // get the MC pt !!
    Float_t pt = ref->Pt();

    TH2 *h = (TH2I*)fHistos->At(kTRDstat);
    if(status & AliESDtrack::kTPCout) h->Fill(pt, kTPCout);
    if(/*status & AliESDtrack::k*/bTRDin) h->Fill(pt, kTRDin);
    if(/*status & AliESDtrack::k*/bTRDout){ 
      ((TH1*)fHistos->At(kNCl))->Fill(esdTrack->GetNcls(2));
      h->Fill(pt, kTRDout);
    }
    if(/*status & AliESDtrack::k*/bTRDpid) h->Fill(pt, kTRDpid);
    if(status & AliESDtrack::kTRDrefit) h->Fill(pt, kTRDref);

    if(ip){
      h = (TH2I*)fHistos->At(kTRDmom);
      Float_t pp(0.);
      for(Int_t ily=6; ily--;){
        if((pp=esdTrack->GetTRDmomentum(ily))<0.) continue;
        h->Fill(ip->GetP()-pp, ily);
      }
    }
  }  
  PostData(0, fHistos);
}

//____________________________________________________________________
TObjArray* AliTRDcheckESD::Histos()
{
// Retrieve histograms array if already build or build it

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

  // energy loss
  if(!(h = (TH2I*)gROOT->FindObject("hTRDmom"))){
    h = new TH2I("hTRDmom", "TRD energy loss", 100, -1., 2., 6, -0.5, 5.5);
    h->GetXaxis()->SetTitle("p_{inner} - p_{ly} [GeV/c]");
    h->GetYaxis()->SetTitle("layer");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fHistos->AddAt(h, kTRDmom);

  return fHistos;
}

//____________________________________________________________________
Bool_t AliTRDcheckESD::Load(const Char_t *filename, const Char_t *name)
{
// Load data from performance file

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


  // geometrical efficiency
  TH2I *h2 = (TH2I*)fHistos->At(kTRDstat);
  TH1 *h1[2] = {0x0, 0x0};
  h1[0] = h2->ProjectionX("checkESDx0", kTPCout, kTPCout);
  h1[1] = h2->ProjectionX("checkESDx1", kTRDin, kTRDin);
  Process(h1, (TGraphErrors*)GetGraph(0));
  delete h1[0];delete h1[1];

  // tracking efficiency
  h1[0] = h2->ProjectionX("checkESDx0", kTRDin, kTRDin);
  h1[1] = h2->ProjectionX("checkESDx1", kTRDout, kTRDout);
  Process(h1, (TGraphErrors*)GetGraph(1));
  delete h1[1];

  // PID efficiency
  h1[1] = h2->ProjectionX("checkESDx1", kTRDpid, kTRDpid);
  Process(h1, (TGraphErrors*)GetGraph(2));
  delete h1[1];

  // Refit efficiency
  h1[1] = h2->ProjectionX("checkESDx1", kTRDref, kTRDref);
  Process(h1, (TGraphErrors*)GetGraph(3));
  delete h1[1];
  if(!(h2 = dynamic_cast<TH2I*>(fHistos->At(kTRDmom)))) return;
 
  TGraphAsymmErrors *g06 = (TGraphAsymmErrors*)GetGraph(4), *g09 = (TGraphAsymmErrors*)GetGraph(5);
  TAxis *ax=h2->GetXaxis();
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

