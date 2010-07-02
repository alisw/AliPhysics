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
#include <TCanvas.h>
#include <TObjArray.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TF1.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TH3S.h>
#include <TH3F.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TParticle.h>
#include <TTimeStamp.h>

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
const Float_t AliTRDcheckESD::fgkQs        = 0.002;

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
void AliTRDcheckESD::MakeSummary(){
  TCanvas *cOut = new TCanvas(Form("summary%s1", GetName()), Form("Summary 1 for task %s", GetName()), 1024, 768);
  cOut->cd();
  GetRefFigure(4);
  cOut->SaveAs(Form("TRDsummary%s1.gif", GetName()));

  cOut = new TCanvas(Form("summary%s2", GetName()), Form("Summary 2 for task %s", GetName()), 1024, 768);
  cOut->cd();
  GetRefFigure(5);
  cOut->SaveAs(Form("TRDsummary%s2.gif", GetName()));
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
  Float_t nada(0.0);
  TH1 *hF(NULL);
  TH1 *hFeffP(NULL); TH1 *hFeffN(NULL);
  TH2 *h2F(NULL); TH2 *h2Feff(NULL);
  TH2 *h2FtpcP(NULL); TH2 *h2FtpcN(NULL);
  TH2 *h2FtrdP(NULL); TH2 *h2FtrdN(NULL);
  TH3 *h3F(NULL);
  if((hF=(TH1S*)gROOT->FindObject("hFcheckESD"))) delete hF;
  TLegend *leg(NULL);
  TList *l(NULL); TVirtualPad *pad(NULL);
  TGraphErrors *g(NULL);TGraphAsymmErrors *ga(NULL);
  TObjArray *arr(NULL);
  TProfile2D *hProf2D(NULL);
  TProfile *hProf(NULL);
  TLatex *lat=new TLatex();
  lat->SetTextSize(0.07);
  lat->SetTextColor(2);
  TLine line;
  TTimeStamp now;
  TF1* fitFunc(NULL);
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
  case 4:            // plot a 3x3 canvas with tracking related histograms
    gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
    gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
    gPad->Divide(3,3,0.,0.);
    l=gPad->GetListOfPrimitives();
    // eta-phi distr. for positive TPC tracks
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE); 
    h3F = dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksPos));
    h2FtpcP = (TH2F*)Proj3D((TH3F*)h3F, 0x0, 1, h3F->GetZaxis()->GetNbins(), nada)->Clone();
    h2FtpcP->GetXaxis()->SetTitle("#eta");
    h2FtpcP->GetXaxis()->CenterTitle();
    h2FtpcP->GetXaxis()->SetTitleSize(0.07);
    h2FtpcP->GetXaxis()->SetTitleOffset(0.8);
    h2FtpcP->GetXaxis()->SetLabelSize(0.05);
    h2FtpcP->GetYaxis()->SetTitle("detector #varphi");
    h2FtpcP->GetYaxis()->CenterTitle();
    h2FtpcP->GetYaxis()->SetTitleSize(0.07);
    h2FtpcP->GetYaxis()->SetTitleOffset(0.8);
    h2FtpcP->GetYaxis()->SetLabelSize(0.05);
    h2FtpcP->SetTitle("");
    h2FtpcP->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TPC positive ref. tracks");
    //-----------------
    // eta-phi distr. for negative TPC tracks
    pad = ((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    h3F = dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksNeg));
    h2FtpcN = (TH2F*)Proj3D((TH3F*)h3F, 0x0, 1, h3F->GetZaxis()->GetNbins(), nada)->Clone();
    h2FtpcN->GetXaxis()->SetTitle("#eta");
    h2FtpcN->GetXaxis()->CenterTitle();
    h2FtpcN->GetXaxis()->SetTitleSize(0.07);
    h2FtpcN->GetXaxis()->SetTitleOffset(0.8);
    h2FtpcN->GetXaxis()->SetLabelSize(0.05);
    h2FtpcN->GetYaxis()->SetTitle("detector #varphi");
    h2FtpcN->GetYaxis()->CenterTitle();
    h2FtpcN->GetYaxis()->SetTitleSize(0.07);
    h2FtpcN->GetYaxis()->SetTitleOffset(0.8);
    h2FtpcN->GetYaxis()->SetLabelSize(0.05);
    h2FtpcN->SetTitle("");
    h2FtpcN->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TPC negative ref. tracks");
    // eta-phi distr. for positive TRD tracks
    pad = ((TVirtualPad*)l->At(3)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    h3F = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos));
    h2FtrdP = (TH2F*)Proj3D((TH3F*)h3F, 0x0, 1, h3F->GetZaxis()->GetNbins(), nada)->Clone();
    h2FtrdP->GetXaxis()->SetTitle("#eta");
    h2FtrdP->GetXaxis()->CenterTitle();
    h2FtrdP->GetXaxis()->SetTitleSize(0.07);
    h2FtrdP->GetXaxis()->SetTitleOffset(0.8);
    h2FtrdP->GetXaxis()->SetLabelSize(0.05);
    h2FtrdP->GetYaxis()->SetTitle("detector #varphi");
    h2FtrdP->GetYaxis()->CenterTitle();
    h2FtrdP->GetYaxis()->SetTitleSize(0.07);
    h2FtrdP->GetYaxis()->SetTitleOffset(0.8);
    h2FtrdP->GetYaxis()->SetLabelSize(0.05);
    h2FtrdP->SetMaximum(h2FtpcP->GetMaximum());
    h2FtrdP->SetTitle("");
    h2FtrdP->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TRD positive ref. tracks");
    //-----------------
    // eta-phi distr. for negative TRD tracks
    pad = ((TVirtualPad*)l->At(4)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    h3F = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg));
    h2FtrdN = (TH2F*)Proj3D((TH3F*)h3F, 0x0, 1, h3F->GetZaxis()->GetNbins(), nada)->Clone();
    h2FtrdN->GetXaxis()->SetTitle("#eta");
    h2FtrdN->GetXaxis()->CenterTitle();
    h2FtrdN->GetXaxis()->SetTitleSize(0.07);
    h2FtrdN->GetXaxis()->SetTitleOffset(0.8);
    h2FtrdN->GetXaxis()->SetLabelSize(0.05);
    h2FtrdN->GetYaxis()->SetTitle("detector #varphi");
    h2FtrdN->GetYaxis()->CenterTitle();
    h2FtrdN->GetYaxis()->SetTitleSize(0.07);
    h2FtrdN->GetYaxis()->SetTitleOffset(0.8);
    h2FtrdN->GetYaxis()->SetLabelSize(0.05);
    h2FtrdN->SetMaximum(h2FtpcN->GetMaximum());
    h2FtrdN->SetTitle("");
    h2FtrdN->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TRD negative ref. tracks");
    // eta-phi efficiency for positive TRD tracks
    pad = ((TVirtualPad*)l->At(6)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    h2Feff = (TH2F*)h2FtrdP->Clone();
    h2Feff->Reset();
    h2Feff->Divide(h2FtrdP, h2FtpcP);
    h2Feff->GetXaxis()->SetTitle("#eta");
    h2Feff->GetXaxis()->CenterTitle();
    h2Feff->GetXaxis()->SetTitleSize(0.07);
    h2Feff->GetXaxis()->SetTitleOffset(0.8);
    h2Feff->GetXaxis()->SetLabelSize(0.05);
    h2Feff->GetYaxis()->SetTitle("detector #varphi");
    h2Feff->GetYaxis()->CenterTitle();
    h2Feff->GetYaxis()->SetTitleSize(0.07);
    h2Feff->GetYaxis()->SetTitleOffset(0.8);
    h2Feff->GetYaxis()->SetLabelSize(0.05);
    h2Feff->SetMaximum(1.0);
    h2Feff->SetTitle("");
    h2Feff->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "Efficiency positive tracks");
    // eta-phi efficiency for negative TRD tracks
    pad = ((TVirtualPad*)l->At(7)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    h2Feff = (TH2F*)h2FtrdN->Clone();
    h2Feff->Reset();
    h2Feff->Divide(h2FtrdN, h2FtpcN);
    h2Feff->GetXaxis()->SetTitle("#eta");
    h2Feff->GetXaxis()->CenterTitle();
    h2Feff->GetXaxis()->SetTitleSize(0.07);
    h2Feff->GetXaxis()->SetTitleOffset(0.8);
    h2Feff->GetXaxis()->SetLabelSize(0.05);
    h2Feff->GetYaxis()->SetTitle("detector #varphi");
    h2Feff->GetYaxis()->CenterTitle();
    h2Feff->GetYaxis()->SetTitleSize(0.07);
    h2Feff->GetYaxis()->SetTitleOffset(0.8);
    h2Feff->GetYaxis()->SetLabelSize(0.05);
    h2Feff->SetMaximum(1.0);
    h2Feff->SetTitle("");
    h2Feff->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "Efficiency negative tracks");
    
    // <ntracklets> vs (phi,eta)
    pad = ((TVirtualPad*)l->At(2)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    hProf2D = dynamic_cast<TProfile2D*>(fHistos->At(kTRDEtaPhiAvNtrkl));
    hProf2D->SetStats(kFALSE);
    hProf2D->SetTitle("");
    hProf2D->GetXaxis()->SetTitle("#eta");
    hProf2D->GetXaxis()->SetTitleOffset(0.8); 
    hProf2D->GetXaxis()->SetTitleSize(0.07);
    hProf2D->GetXaxis()->CenterTitle();
    hProf2D->GetXaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->SetTitle("detector #varphi");
    hProf2D->GetYaxis()->SetTitleOffset(0.8); 
    hProf2D->GetYaxis()->SetTitleSize(0.07);
    hProf2D->GetYaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->CenterTitle();
    hProf2D->SetMinimum(0.);
    hProf2D->SetMaximum(6.);
    hProf2D->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TRD <N_{tracklets}>");
    // TPC-TRD matching efficiency vs pt
    pad = ((TVirtualPad*)l->At(5)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
    pad->SetTopMargin(0.02); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    hFeffP = TRDEfficiency(+1);
    hFeffN = TRDEfficiency(-1);
    h2F=new TH2F("rangeEffPt", "",10,0.,10.,10,0.,1.1);
    h2F->SetStats(kFALSE);
    h2F->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h2F->GetXaxis()->SetTitleOffset(0.8); 
    h2F->GetXaxis()->SetTitleSize(0.07);
    h2F->GetXaxis()->CenterTitle();
    h2F->GetXaxis()->SetLabelSize(0.05);
    h2F->GetYaxis()->SetTitle("TRD-TPC matching efficiency");
    h2F->GetYaxis()->SetTitleOffset(0.8); 
    h2F->GetYaxis()->SetTitleSize(0.07);
    h2F->GetYaxis()->SetLabelSize(0.05);
    h2F->GetYaxis()->CenterTitle();
    h2F->Draw();
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.DrawLine(h2F->GetXaxis()->GetXmin(), 0.7, h2F->GetXaxis()->GetXmax(), 0.7);
    line.DrawLine(h2F->GetXaxis()->GetXmin(), 0.9, h2F->GetXaxis()->GetXmax(), 0.9);
    hFeffP->SetMarkerStyle(20);
    hFeffP->SetMarkerColor(2);
    hFeffN->SetMarkerStyle(22);
    hFeffN->SetMarkerColor(4);
    hFeffP->Draw("same");
    hFeffN->Draw("same");
    leg=new TLegend(0.65, 0.2, 0.95, 0.4);
    leg->SetFillColor(0);
    leg->AddEntry(hFeffP, "positives", "p");
    leg->AddEntry(hFeffN, "negatives", "p");
    leg->Draw();
    // create trending values for the TPC-TRD matching efficiency
    // fit the efficiency histos with a constant in the range [1.0,1.5] GeV/c
    fitFunc = new TF1("constantFunc","[0]",1.0,1.5);
    hFeffP->Fit(fitFunc,"Q0","",1.0,1.5);
    PutTrendValue("TrackingEffPos1GeV", fitFunc->GetParameter(0));
    PutTrendValue("TrackingEffPos1GeVErr", fitFunc->GetParError(0));
    hFeffN->Fit(fitFunc,"Q0","",1.0,1.5);
    PutTrendValue("TrackingEffNeg1GeV", fitFunc->GetParameter(0));
    PutTrendValue("TrackingEffNeg1GeVErr", fitFunc->GetParError(0));
    
    // Nclusters per TRD track
    pad = ((TVirtualPad*)l->At(8)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.12);
    pad->SetTopMargin(0.02); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    pad->SetLogz();
    h2F = dynamic_cast<TH2F*>(fHistos->At(kNClsTrackTRD));
    h2F->SetStats(kFALSE);
    h2F->SetTitle("");
    h2F->GetXaxis()->SetTitle("p [GeV/c]");
    h2F->GetXaxis()->SetTitleOffset(0.8); 
    h2F->GetXaxis()->SetTitleSize(0.07);
    h2F->GetXaxis()->CenterTitle();
    h2F->GetXaxis()->SetLabelSize(0.05);
    h2F->GetYaxis()->SetTitle("#clusters per TRD track");
    h2F->GetYaxis()->SetTitleOffset(0.8); 
    h2F->GetYaxis()->SetTitleSize(0.07);
    h2F->GetYaxis()->CenterTitle();
    h2F->GetYaxis()->SetLabelSize(0.05);
    h2F->Draw("colz");
    break;
  case 5:            // plot a 3x3 canvas with PID related histograms
    gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
    gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
    gPad->Divide(3,3,0.,0.);
    l=gPad->GetListOfPrimitives();
    // eta-phi distr. for <Qtot> in layer 0
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    hProf2D = dynamic_cast<TProfile2D*>(fHistos->At(kTRDEtaPhiAvQtot+0));
    hProf2D->SetStats(kFALSE);
    hProf2D->SetTitle("");
    hProf2D->GetXaxis()->SetTitle("#eta");
    hProf2D->GetXaxis()->SetTitleOffset(0.8); 
    hProf2D->GetXaxis()->SetTitleSize(0.07);
    hProf2D->GetXaxis()->CenterTitle();
    hProf2D->GetXaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->SetTitle("detector #varphi");
    hProf2D->GetYaxis()->SetTitleOffset(0.8); 
    hProf2D->GetYaxis()->SetTitleSize(0.07);
    hProf2D->GetYaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->CenterTitle();
    hProf2D->SetMinimum(0.);
    hProf2D->SetMaximum(25.);
    hProf2D->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TRD <Q_{tot}> Layer 0");
    // eta-phi distr. for <Qtot> in layer 1
    pad = ((TVirtualPad*)l->At(3)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    hProf2D = dynamic_cast<TProfile2D*>(fHistos->At(kTRDEtaPhiAvQtot+1));
    hProf2D->SetStats(kFALSE);
    hProf2D->SetTitle("");
    hProf2D->GetXaxis()->SetTitle("#eta");
    hProf2D->GetXaxis()->SetTitleOffset(0.8); 
    hProf2D->GetXaxis()->SetTitleSize(0.07);
    hProf2D->GetXaxis()->CenterTitle();
    hProf2D->GetXaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->SetTitle("detector #varphi");
    hProf2D->GetYaxis()->SetTitleOffset(0.8); 
    hProf2D->GetYaxis()->SetTitleSize(0.07);
    hProf2D->GetYaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->CenterTitle();
    hProf2D->SetMinimum(0.);
    hProf2D->SetMaximum(25.);
    hProf2D->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TRD <Q_{tot}> Layer 1");
    // eta-phi distr. for <Qtot> in layer 2
    pad = ((TVirtualPad*)l->At(6)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    hProf2D = dynamic_cast<TProfile2D*>(fHistos->At(kTRDEtaPhiAvQtot+2));
    hProf2D->SetStats(kFALSE);
    hProf2D->SetTitle("");
    hProf2D->GetXaxis()->SetTitle("#eta");
    hProf2D->GetXaxis()->SetTitleOffset(0.8); 
    hProf2D->GetXaxis()->SetTitleSize(0.07);
    hProf2D->GetXaxis()->CenterTitle();
    hProf2D->GetXaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->SetTitle("detector #varphi");
    hProf2D->GetYaxis()->SetTitleOffset(0.8); 
    hProf2D->GetYaxis()->SetTitleSize(0.07);
    hProf2D->GetYaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->CenterTitle();
    hProf2D->SetMinimum(0.);
    hProf2D->SetMaximum(25.);
    hProf2D->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TRD <Q_{tot}> Layer 2");
    // eta-phi distr. for <Qtot> in layer 3
    pad = ((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    hProf2D = dynamic_cast<TProfile2D*>(fHistos->At(kTRDEtaPhiAvQtot+3));
    hProf2D->SetStats(kFALSE);
    hProf2D->SetTitle("");
    hProf2D->GetXaxis()->SetTitle("#eta");
    hProf2D->GetXaxis()->SetTitleOffset(0.8); 
    hProf2D->GetXaxis()->SetTitleSize(0.07);
    hProf2D->GetXaxis()->CenterTitle();
    hProf2D->GetXaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->SetTitle("detector #varphi");
    hProf2D->GetYaxis()->SetTitleOffset(0.8); 
    hProf2D->GetYaxis()->SetTitleSize(0.07);
    hProf2D->GetYaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->CenterTitle();
    hProf2D->SetMinimum(0.);
    hProf2D->SetMaximum(25.);
    hProf2D->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TRD <Q_{tot}> Layer 3");
    // eta-phi distr. for <Qtot> in layer 4
    pad = ((TVirtualPad*)l->At(4)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    hProf2D = dynamic_cast<TProfile2D*>(fHistos->At(kTRDEtaPhiAvQtot+4));
    hProf2D->SetStats(kFALSE);
    hProf2D->SetTitle("");
    hProf2D->GetXaxis()->SetTitle("#eta");
    hProf2D->GetXaxis()->SetTitleOffset(0.8); 
    hProf2D->GetXaxis()->SetTitleSize(0.07);
    hProf2D->GetXaxis()->CenterTitle();
    hProf2D->GetXaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->SetTitle("detector #varphi");
    hProf2D->GetYaxis()->SetTitleOffset(0.8); 
    hProf2D->GetYaxis()->SetTitleSize(0.07);
    hProf2D->GetYaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->CenterTitle();
    hProf2D->SetMinimum(0.);
    hProf2D->SetMaximum(25.);
    hProf2D->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TRD <Q_{tot}> Layer 4");
    // eta-phi distr. for <Qtot> in layer 5
    pad = ((TVirtualPad*)l->At(7)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    hProf2D = dynamic_cast<TProfile2D*>(fHistos->At(kTRDEtaPhiAvQtot+5));
    hProf2D->SetStats(kFALSE);
    hProf2D->SetTitle("");
    hProf2D->GetXaxis()->SetTitle("#eta");
    hProf2D->GetXaxis()->SetTitleOffset(0.8); 
    hProf2D->GetXaxis()->SetTitleSize(0.07);
    hProf2D->GetXaxis()->CenterTitle();
    hProf2D->GetXaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->SetTitle("detector #varphi");
    hProf2D->GetYaxis()->SetTitleOffset(0.8); 
    hProf2D->GetYaxis()->SetTitleSize(0.07);
    hProf2D->GetYaxis()->SetLabelSize(0.05);
    hProf2D->GetYaxis()->CenterTitle();
    hProf2D->SetMinimum(0.);
    hProf2D->SetMaximum(25.);
    hProf2D->Draw("colz");
    lat->DrawLatex(-0.9, 3.6, "TRD <Q_{tot}> Layer 5");
    // PH versus slice number
    pad = ((TVirtualPad*)l->At(2)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    h2F = dynamic_cast<TH2F*>(fHistos->At(kPHSlice));
    h2F->SetStats(kFALSE);
    h2F->SetTitle("");
    h2F->GetXaxis()->SetTitle("slice");
    h2F->GetXaxis()->SetTitleOffset(0.8); 
    h2F->GetXaxis()->SetTitleSize(0.07);
    h2F->GetXaxis()->CenterTitle();
    h2F->GetXaxis()->SetLabelSize(0.05);
    h2F->GetYaxis()->SetTitle("PH");
    h2F->GetYaxis()->SetTitleOffset(0.8); 
    h2F->GetYaxis()->SetTitleSize(0.07);
    h2F->GetYaxis()->SetLabelSize(0.05);
    h2F->GetYaxis()->CenterTitle();
    h2F->Draw("colz");
    //hProf = h2F->ProfileX("profileX");
    //hProf->SetLineWidth(2);
    //hProf->Draw("same");
    // Qtot vs P
    pad = ((TVirtualPad*)l->At(5)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    pad->SetLogz();
    h2F = dynamic_cast<TH2F*>(fHistos->At(kQtotP));
    h2F->SetStats(kFALSE);
    h2F->SetTitle("");
    h2F->GetXaxis()->SetTitle("P [GeV/c]");
    h2F->GetXaxis()->SetTitleOffset(0.8); 
    h2F->GetXaxis()->SetTitleSize(0.07);
    h2F->GetXaxis()->CenterTitle();
    h2F->GetXaxis()->SetLabelSize(0.05);
    h2F->GetYaxis()->SetRangeUser(0.0,100.0);
    h2F->GetYaxis()->SetTitle("Q_{tot}");
    h2F->GetYaxis()->SetTitleOffset(0.8); 
    h2F->GetYaxis()->SetTitleSize(0.07);
    h2F->GetYaxis()->SetLabelSize(0.05);
    h2F->GetYaxis()->CenterTitle();
    h2F->Draw("colz");
    // create trending value for the average Qtot at 1 GeV/c
    hProf = h2F->ProfileX("profileQtot",1,h2F->GetYaxis()->FindBin(40.));
    PutTrendValue("AvQtot1GeV", hProf->GetBinContent(hProf->GetXaxis()->FindBin(1.)));
    PutTrendValue("AvQtot1GeVErr", hProf->GetBinError(hProf->GetXaxis()->FindBin(1.)));
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
  TH1 *h(NULL);
  
  // fill event vertex histos
  h = (TH1F*)fHistos->At(kTPCVertex);
  if(fESD->GetPrimaryVertexTPC()) h->Fill(fESD->GetPrimaryVertexTPC()->GetZv());
  h = (TH1F*)fHistos->At(kEventVertex);
  if(fESD->GetPrimaryVertex()) h->Fill(fESD->GetPrimaryVertex()->GetZv());
  // fill the uncutted number of tracks
  h = (TH1I*)fHistos->At(kNTracksAll);
  h->Fill(fESD->GetNumberOfTracks());
  
  // counters for number of tracks in acceptance&DCA and for those with a minimum of TPC clusters
  Int_t nTracksAcc=0;
  Int_t nTracksTPC=0;
  
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
    
    Float_t par[2], cov[3];
    esdTrack->GetImpactParameters(par, cov);
    if(selected && esdTrack->GetTPCNcls()>=10) {
      // fill DCA histograms
      h = (TH1F*)fHistos->At(kDCAxy); h->Fill(par[0]);
      h = (TH1F*)fHistos->At(kDCAz); h->Fill(par[1]);
      // fill pt distribution at this stage
      h = (TH1F*)fHistos->At(kPt1); h->Fill(esdTrack->Pt());
    }
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
    Float_t theta=esdTrack->Theta();
    Float_t phi=esdTrack->Phi();
    Int_t nClustersTPC = esdTrack->GetTPCNcls();
    Float_t eta=-TMath::Log(TMath::Tan(theta/2.));
    if(selected) {
      nTracksAcc++;   // number of tracks in acceptance and DCA cut
      // fill pt distribution at this stage
      h = (TH1F*)fHistos->At(kPt2); h->Fill(esdTrack->Pt());
      // TPC nclusters distribution
      h = (TH1I*)fHistos->At(kNTPCCl); h->Fill(nClustersTPC);
      if(esdTrack->Pt()>1.0) {
        h = (TH1I*)fHistos->At(kNTPCCl2); h->Fill(nClustersTPC);
      }
      // (eta,nclustersTPC) distrib of TPC ref. tracks
      h = (TH2F*)fHistos->At(kEtaNclsTPC); h->Fill(eta, nClustersTPC);
      // (phi,nclustersTPC) distrib of TPC ref. tracks
      h = (TH2F*)fHistos->At(kPhiNclsTPC); h->Fill(phi, nClustersTPC);
      
    }
      
    if(nClustersTPC < fgkNclTPC){ 
      AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] NclTPC[%d]", fESD->GetEventNumberInFile(), itrk, nClustersTPC));
      selected = kFALSE;
    }
    if(!selected) continue;
    
    // number of TPC reference tracks
    nTracksTPC++;
    
    Int_t nTRD(esdTrack->GetNcls(2));
    Double_t pt(esdTrack->Pt());
    Double_t p[AliPID::kSPECIES]; esdTrack->GetTRDpid(p);
    // pid quality
    Bool_t kBarrel = Bool_t(status & AliESDtrack::kTRDin);

    TH3F *hhh;
    // find position and momentum of the track at entrance in TRD
    Double_t localCoord[3] = {0., 0., 0.};
    Bool_t localCoordGood = esdTrack->GetXYZAt(298., fESD->GetMagneticField(), localCoord);
    if(localCoordGood) {
      hhh = (TH3F*)fHistos->At(kPropagXYvsP); hhh->Fill(localCoord[0], localCoord[1], esdTrack->GetP());
      hhh = (TH3F*)fHistos->At(kPropagRZvsP); hhh->Fill(localCoord[2], TMath::Sqrt(localCoord[0]*localCoord[0]+localCoord[1]*localCoord[1]), esdTrack->GetP());
    }
    Double_t localMom[3] = {0., 0., 0.};
    Bool_t localMomGood = esdTrack->GetPxPyPzAt(298., fESD->GetMagneticField(), localMom);
    Double_t localPhi = (localMomGood ? TMath::ATan2(localMom[1], localMom[0]) : 0.0);
    Double_t localSagitaPhi = (localCoordGood ? TMath::ATan2(localCoord[1], localCoord[0]) : 0.0);

    // fill pt distribution at this stage
    if(esdTrack->Charge()>0) {
      h = (TH1F*)fHistos->At(kPt3pos); h->Fill(pt);
      // fill eta-phi map of TPC positive ref. tracks
      if(localCoordGood && localMomGood) {
        hhh = (TH3F*)fHistos->At(kTPCRefTracksPos); hhh->Fill(eta, localSagitaPhi, pt);
      }
    }
    if(esdTrack->Charge()<0) {
      h = (TH1F*)fHistos->At(kPt3neg); h->Fill(pt);
      // fill eta-phi map of TPC negative ref. tracks
      if(localCoordGood && localMomGood) {
        hhh = (TH3F*)fHistos->At(kTPCRefTracksNeg); hhh->Fill(eta, localSagitaPhi, pt);
      }
    }
    // TPC dE/dx vs P
    h = (TH2F*)fHistos->At(kTPCDedx); h->Fill(esdTrack->GetP(), esdTrack->GetTPCsignal());
    // (eta,phi) distrib of TPC ref. tracks
    h = (TH2F*)fHistos->At(kEtaPhi); h->Fill(eta, phi);
        
    Int_t nTRDtrkl = esdTrack->GetTRDntracklets();
    // TRD reference tracks
    if(nTRDtrkl>=1) {
      // fill pt distribution at this stage
      if(esdTrack->Charge()>0) {
        h = (TH1F*)fHistos->At(kPt4pos); h->Fill(pt);
	// fill eta-phi map of TRD positive ref. tracks
	if(localCoordGood && localMomGood) {
          hhh = (TH3F*)fHistos->At(kTRDRefTracksPos); hhh->Fill(eta, localSagitaPhi, pt);
        }
      }
      if(esdTrack->Charge()<0) {
        h = (TH1F*)fHistos->At(kPt4neg); h->Fill(pt);
	// fill eta-phi map of TRD negative ref. tracks
	if(localCoordGood && localMomGood) {
          hhh = (TH3F*)fHistos->At(kTRDRefTracksNeg); hhh->Fill(eta, localSagitaPhi, pt);
        }
      }
      TProfile2D *h2d;
      // fill eta-phi map of TRD negative ref. tracks
      if(localCoordGood && localMomGood) {
        h2d = (TProfile2D*)fHistos->At(kTRDEtaPhiAvNtrkl); h2d->Fill(eta, localSagitaPhi, (Float_t)nTRDtrkl);
	h2d = (TProfile2D*)fHistos->At(kTRDEtaDeltaPhiAvNtrkl); h2d->Fill(eta, localPhi-localSagitaPhi, (Float_t)nTRDtrkl);
      }
      // ntracklets/track vs P
      h = (TH2F*)fHistos->At(kNTrackletsTRD); h->Fill(esdTrack->GetP(), nTRDtrkl);
      // ntracklets/track vs P
      h = (TH2F*)fHistos->At(kNClsTrackTRD); h->Fill(esdTrack->GetP(), esdTrack->GetTRDncls());
      // (slicePH,sliceNo) distribution and Qtot from slices
      for(Int_t iPlane=0; iPlane<6; iPlane++) {
        Float_t Qtot=0;
        for(Int_t iSlice=0; iSlice<8; iSlice++) {
	  if(esdTrack->GetTRDslice(iPlane, iSlice)>20.) {
	    h = (TH2F*)fHistos->At(kPHSlice); h->Fill(iSlice, esdTrack->GetTRDslice(iPlane, iSlice));
	    Qtot += esdTrack->GetTRDslice(iPlane, iSlice);
	  }
        }
        // Qtot>100 to avoid noise
        if(Qtot>100.) {
          h = (TH2F*)fHistos->At(kQtotP); h->Fill(esdTrack->GetTRDmomentum(iPlane), fgkQs*Qtot);
	}
	// Qtot>100 to avoid noise
	// fgkQs*Qtot<40. so that the average will give a value close to the peak
	if(localCoordGood && localMomGood && Qtot>100. && fgkQs*Qtot<40.) {
	  h2d = (TProfile2D*)fHistos->At(kTRDEtaPhiAvQtot+iPlane);
	  h2d->Fill(eta, localSagitaPhi, fgkQs*Qtot);
	}
      }
      // theta distribution
      h = (TH1F*)fHistos->At(kTheta); h->Fill(theta);
      h = (TH1F*)fHistos->At(kPhi); h->Fill(phi);
    }  // end if nTRDtrkl>=1
    
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
    }     // end if(HasMC())

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
  }  // end loop over tracks
  
  // fill the number of tracks histograms
  h = (TH1I*)fHistos->At(kNTracksAcc);
  h->Fill(nTracksAcc);
  h = (TH1I*)fHistos->At(kNTracksTPC);
  h->Fill(nTracksTPC);
  
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
  //if(!HasMC()) return fHistos;

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

  // TPC event vertex distribution
  if(!(h = (TH1F*)gROOT->FindObject("hTPCVertex"))){
    h = new TH1F("hTPCVertex", "Event vertex Z coord. from TPC tracks", 100, -25., 25.);
  } else h->Reset();
  fHistos->AddAt(h, kTPCVertex);
  
  // Event vertex
  if(!(h = (TH1F*)gROOT->FindObject("hEventVertex"))){
    h = new TH1F("hEventVertex", "Event vertex Z coord.", 100, -25., 25.);
  } else h->Reset();
  fHistos->AddAt(h, kEventVertex);
  
  // Number of all tracks
  if(!(h = (TH1I*)gROOT->FindObject("hNTracksAll"))){
    h = new TH1I("hNTracksAll", "Number of tracks per event, event vertex cuts", 5000, 0, 5000);
  } else h->Reset();
  fHistos->AddAt(h, kNTracksAll);
  
  // Number of tracks in acceptance and DCA cut
  if(!(h = (TH1I*)gROOT->FindObject("hNTracksAcc"))){
    h = new TH1I("hNTracksAcc", Form("Number of tracks per event, |#eta|<%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
				     fgkEta, fgkTrkDCAxy, fgkTrkDCAz), 5000, 0, 5000);
  } else h->Reset();
  fHistos->AddAt(h, kNTracksAcc);
  
  // Number of tracks in TPC (Ncls>10)
  if(!(h = (TH1I*)gROOT->FindObject("hNTracksTPC"))){
    h = new TH1I("hNTracksTPC", Form("Number of tracks per event, |#eta|<%.1f, pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
				     fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 5000, 0, 5000);
  } else h->Reset();
  fHistos->AddAt(h, kNTracksTPC);
  
  // Distribution of DCA-xy
  if(!(h = (TH1F*)gROOT->FindObject("hDCAxy"))){
    h = new TH1F("hDCAxy", "Distribution of transverse DCA", 100, -100., 100.);
  } else h->Reset();
  fHistos->AddAt(h, kDCAxy);
  
  // Distribution of DCA-z
  if(!(h = (TH1F*)gROOT->FindObject("hDCAz"))){
    h = new TH1F("hDCAz", "Distribution of longitudinal DCA", 100, -100., 100.);
  } else h->Reset();
  fHistos->AddAt(h, kDCAz);
  
  Float_t binPtLimits[33] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
		 	     1.0, 1.1, 1.2, 1.3, 1.4, 
			     1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
			     3.4, 3.8, 4.2, 4.6, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  // Pt distributions
  if(!(h = (TH1F*)gROOT->FindObject("hPt1"))){
    h = new TH1F("hPt1", Form("dN/dpt, |#eta|<%.1f and pt>%.1f", fgkEta, fgkPt), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt1);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt2"))){
    h = new TH1F("hPt2", Form("dN/dpt, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
			      fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt2);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt3pos"))){
    h = new TH1F("hPt3pos", Form("dN/dpt (positives), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt3pos);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt3neg"))){
    h = new TH1F("hPt3neg", Form("dN/dpt (negatives), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt3neg);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt4pos"))){
    h = new TH1F("hPt4pos", Form("dN/dpt (positives), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq 1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt4pos);
  
  if(!(h = (TH1F*)gROOT->FindObject("hPt4neg"))){
    h = new TH1F("hPt4pos", Form("dN/dpt (negatives), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq 1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kPt4neg);
  
  // theta distribution of TRD tracks
  if(!(h = (TH1F*)gROOT->FindObject("hTheta"))){
    h = new TH1F("hTheta", Form("dN/d#theta, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq 1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 220,.5,2.7);
  } else h->Reset();
  fHistos->AddAt(h, kTheta);
  
  // phi distribution of TRD tracks
  if(!(h = (TH1F*)gROOT->FindObject("hPhi"))){
    h = new TH1F("hPhi", Form("dN/d#varphi, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq 1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 157,0,6.28);
  } else h->Reset();
  fHistos->AddAt(h, kPhi);
  
  // TPC cluster distribution
  if(!(h = (TH1F*)gROOT->FindObject("hNTPCCl"))){
    h = new TH1I("hNTPCCl", Form("Number of TPC clusters/track, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 160, 0, 160);
  } else h->Reset();
  fHistos->AddAt(h, kNTPCCl);
  
  if(!(h = (TH1I*)gROOT->FindObject("hNTPCCl2"))){
    h = new TH1F("hNTPCCl2", Form("Number of TPC clusters/track, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, pt>1.0 GeV/c",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 160, 0, 160);
  } else h->Reset();
  fHistos->AddAt(h, kNTPCCl2);
  
  // dE/dx vs P for TPC reference tracks
  if(!(h = (TH2F*)gROOT->FindObject("hTPCDedx"))){
    h = new TH2F("hTPCDedx", Form("TPC dE/dx vs P, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, 0.1,10.1, 120, 0,600.);
  } else h->Reset();
  fHistos->AddAt(h, kTPCDedx);
  
  // eta,phi distribution of TPC reference tracks
  if(!(h = (TH2F*)gROOT->FindObject("hEtaPhi"))){
    h = new TH2F("hEtaPhi", Form("TPC (#eta,#varphi), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 50, -1, 1, 157, 0, 6.28);
  } else h->Reset();
  fHistos->AddAt(h, kEtaPhi);
  
  // Nclusters vs eta distribution for TPC tracks
  if(!(h = (TH2F*)gROOT->FindObject("hEtaNclsTPC"))){
    h = new TH2F("hEtaNclsTPC", Form("TPC Nclusters vs. #eta, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 50, -1, 1, 160, 0, 160.);
  } else h->Reset();
  fHistos->AddAt(h, kEtaNclsTPC);
  
  // Nclusters vs phi distribution for TPC reference tracks
  if(!(h = (TH2F*)gROOT->FindObject("hPhiNclsTPC"))){
    h = new TH2F("hPhiNclsTPC", Form("TPC Nclusters vs. #varphi, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz), 157, 0, 6.28, 160, 0, 160.);
  } else h->Reset();
  fHistos->AddAt(h, kPhiNclsTPC);
  
  // Ntracklets/track vs P for TRD reference tracks
  Double_t binsP[19] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0,
			2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0};
  if(!(h = (TH2F*)gROOT->FindObject("hNTrackletsTRD"))){
    h = new TH2F("hNTrackletsTRD", Form("TRD Ntracklets/track vs. P, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 18, binsP, 7, -0.5, 6.5);
  } else h->Reset();
  fHistos->AddAt(h, kNTrackletsTRD);
  
  // Nclusters/track vs P for TRD reference tracks
  if(!(h = (TH2F*)gROOT->FindObject("hNClsTrackTRD"))){
    h = new TH2F("hNClsTrackTRD", Form("TRD Nclusters/track vs. P, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 18, binsP, 180, 0., 180.);
  } else h->Reset();
  fHistos->AddAt(h, kNClsTrackTRD);
  
  // <PH> vs slice number for TRD reference tracklets
  if(!(h = (TH2F*)gROOT->FindObject("hPHSlice"))){
    h = new TH2F("hPHSlice", Form("<PH> vs sliceNo, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 8, -0.5, 7.5, 2000, 0., 2000.);
  } else h->Reset();
  fHistos->AddAt(h, kPHSlice);
  
  // Qtot vs P for TRD reference tracklets
  if(!(h = (TH2F*)gROOT->FindObject("hQtotP"))){
    h = new TH2F("hQtotP", Form("Qtot(from slices) vs P, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 18, binsP, 400, 0., 200);
  } else h->Reset();
  fHistos->AddAt(h, kQtotP);
  
  // (X,Y,momentum) distribution after AliESDtrack::PropagateTo(r=300.)
  if(!(h = (TH3F*)gROOT->FindObject("hPropagXYvsP"))){
    h = new TH3F("hPropagXYvsP", Form("(x,y) vs P after AliESDtrack::PropagateTo(r=300.), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100,-500,500, 100,-500,500, 10, 0.,10.);
  } else h->Reset();
  fHistos->AddAt(h, kPropagXYvsP);
  
  // (R,Z,momentum) distribution after AliESDtrack::PropagateTo(r=300.)
  if(!(h = (TH3F*)gROOT->FindObject("hPropagRZvsP"))){
    h = new TH3F("hPropagRZvsP", Form("(r,z) vs P after AliESDtrack::PropagateTo(r=300.), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100,-350., 350., 100,0.,500., 10, 0.,10.);
  } else h->Reset();
  fHistos->AddAt(h, kPropagRZvsP);
  
  Float_t etaBinLimits[101];	
  for(Int_t i=0; i<101; i++) etaBinLimits[i] = -1.0 + i*2.0/100.;
  Float_t phiBinLimits[151];
  for(Int_t i=0; i<151; i++) phiBinLimits[i] = -1.1*TMath::Pi() + i*2.2*TMath::Pi()/150.;
  // (eta,detector phi,P) distribution of reference TPC positive tracks
  if(!(h = (TH3F*)gROOT->FindObject("hTPCRefTracksPos"))){
    h = new TH3F("hTPCRefTracksPos", Form("(#eta,detector #varphi,p) for TPC positive reference tracks, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kTPCRefTracksPos);
  
  // (eta,detector phi,P) distribution of reference TPC negative tracks
  if(!(h = (TH3F*)gROOT->FindObject("hTPCRefTracksNeg"))){
    h = new TH3F("hTPCRefTracksNeg", Form("(#eta,detector #varphi,p) for TPC negative reference tracks, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kTPCRefTracksNeg);
  
  // (eta,detector phi,P) distribution of reference TRD positive tracks
  if(!(h = (TH3F*)gROOT->FindObject("hTRDRefTracksPos"))){
    h = new TH3F("hTRDRefTracksPos", Form("(#eta,detector #varphi,p) for TRD positive reference tracks, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kTRDRefTracksPos);
  
  // (eta,detector phi,P) distribution of reference TRD negative tracks
  if(!(h = (TH3F*)gROOT->FindObject("hTRDRefTracksNeg"))){
    h = new TH3F("hTRDRefTracksNeg", Form("(#eta,detector #varphi,p) for TRD negative reference tracks, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, etaBinLimits, 150, phiBinLimits, 32, binPtLimits);
  } else h->Reset();
  fHistos->AddAt(h, kTRDRefTracksNeg);
  
  // (eta,detector phi) profile of average number of TRD tracklets/track
  if(!(h = (TProfile2D*)gROOT->FindObject("hTRDEtaPhiAvNtrkl"))){
    h = new TProfile2D("hTRDEtaPhiAvNtrkl", Form("<Ntracklets/track> vs (#eta,detector #varphi) for TRD reference tracks, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  } else h->Reset();
  fHistos->AddAt(h, kTRDEtaPhiAvNtrkl);

  // (eta,delta phi) profile of average number of TRD tracklets/track
  if(!(h = (TProfile2D*)gROOT->FindObject("hTRDEtaDeltaPhiAvNtrkl"))){
    h = new TProfile2D("hTRDEtaDeltaPhiAvNtrkl", Form("<Ntracklets/track> vs (#eta, #Delta#varphi) for TRD reference tracks, |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
		 fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, -1.0, 1.0, 50, -0.4*TMath::Pi(), 0.4*TMath::Pi());
  } else h->Reset();
  fHistos->AddAt(h, kTRDEtaDeltaPhiAvNtrkl);
  
  // (eta, detector phi) profile of average tracklet Qtot from slices
  for(Int_t iLayer=0;iLayer<6;iLayer++) {
    if(!(h = (TProfile2D*)gROOT->FindObject(Form("hTRDEtaPhiAvQtot_Layer%d",iLayer)))) {
      h = new TProfile2D(Form("hTRDEtaPhiAvQtot_Layer%d",iLayer),
			 Form("<Q_{tot}> vs (#eta, detector #varphi) for TRD reference tracks (layer %d), |#eta|<%.1f and pt>%.1f, |DCAxy|<%.1f, |DCAz|<%.1f, TPC nclusters>%d, nTRDtracklets#geq1",
		              iLayer, fgkEta, fgkPt, fgkTrkDCAxy, fgkTrkDCAz, fgkNclTPC), 100, -1.0, 1.0, 50, -1.1*TMath::Pi(), 1.1*TMath::Pi());
    } else h->Reset();
    fHistos->AddAt(h, kTRDEtaPhiAvQtot+iLayer);
  }
  
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
//  if(!HasMC()) return;

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
  
  // 3x3 tracking summary canvas
  fNRefFigures++;
  // 3x3 PID summary canvas
  fNRefFigures++;
}

//____________________________________________________________________
Int_t AliTRDcheckESD::Pdg2Idx(Int_t pdg)
{
  switch(pdg){
  case kElectron: 
  case kPositron: return AliPID::kElectron;  
  case kMuonPlus:
  case kMuonMinus: return AliPID::kMuon;  
  case kPiPlus: 
  case kPiMinus: return AliPID::kPion;  
  case kKPlus: 
  case kKMinus: return AliPID::kKaon;
  case kProton: 
  case kProtonBar: return AliPID::kProton;
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
//____________________________________________________________________
TH2F* AliTRDcheckESD::Proj3D(TH3F* hist, TH2F* accMap, Int_t zbinLow, Int_t zbinHigh, Float_t &entries) {
  //
  //  Project a 3D histogram to a 2D histogram in the Z axis interval [zbinLow,zbinHigh] 
  //  Return the 2D histogram and also the number of entries into this projection (entries)

  Int_t nBinsX = hist->GetXaxis()->GetNbins();   // X and Y axis bins are assumed to be all equal
  Float_t minX = hist->GetXaxis()->GetXmin();
  Float_t maxX = hist->GetXaxis()->GetXmax();
  Int_t nBinsY = hist->GetYaxis()->GetNbins();
  Float_t minY = hist->GetYaxis()->GetXmin();
  Float_t maxY = hist->GetYaxis()->GetXmax();
  Int_t nBinsZ = hist->GetZaxis()->GetNbins();  // Z axis bins (pt) might have different widths
  //Float_t minZ = hist->GetZaxis()->GetXmin();
  //Float_t maxZ = hist->GetZaxis()->GetXmax();

  TH2F* projHisto = (TH2F*)gROOT->FindObject("projHisto");
  if(projHisto) 
    projHisto->Reset();
  else
    projHisto = new TH2F("projHisto", "projection", nBinsX, minX, maxX, nBinsY, minY, maxY);

  const Double_t kMinAccFraction = 0.1;
  entries = 0.0;
  Double_t maxAcc = (accMap ? accMap->GetMaximum() : -1);
  
  for(Int_t iZ=1; iZ<=nBinsZ; iZ++) {
    if(iZ<zbinLow) continue;
    if(iZ>zbinHigh) continue;
    for(Int_t iX=1; iX<=nBinsX; iX++) {
      for(Int_t iY=1; iY<=nBinsY; iY++) {
	if(accMap && maxAcc>0) {
	  if(accMap->GetBinContent(iX,iY)/maxAcc>kMinAccFraction)
	    projHisto->SetBinContent(iX, iY, projHisto->GetBinContent(iX, iY)+hist->GetBinContent(iX,iY,iZ));
	}
	else    // no acc. cut 
	  projHisto->SetBinContent(iX, iY, projHisto->GetBinContent(iX, iY)+hist->GetBinContent(iX,iY,iZ));
	// count only the entries which are inside the acceptance map
	if(accMap && maxAcc>0) { 
	  if(accMap->GetBinContent(iX,iY)/maxAcc>kMinAccFraction)
	    entries+=hist->GetBinContent(iX,iY,iZ);
	}
	else    // no acc. cut
	  entries+=hist->GetBinContent(iX,iY,iZ);
      }
    }
  }
  return projHisto;
}
//____________________________________________________________________
TH1F* AliTRDcheckESD::TRDEfficiency(Short_t positives) {
  //
  // Calculate the TRD-TPC matching efficiency as function of pt
  //
  TH3F* tpc3D(NULL); TH3F* trd3D(NULL);
  if(positives>0) {      // get the histos for positive tracks  
    if(!(tpc3D = dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksPos)))) return 0x0;
    if(!(trd3D = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksPos)))) return 0x0;
  }
  else {      // get the histos for positive tracks
    if(!(tpc3D = dynamic_cast<TH3F*>(fHistos->At(kTPCRefTracksNeg)))) return 0x0;
    if(!(trd3D = dynamic_cast<TH3F*>(fHistos->At(kTRDRefTracksNeg)))) return 0x0;
  }
  
  Int_t nBinsZ = trd3D->GetZaxis()->GetNbins();
  // project everything on the eta-phi map to obtain an acceptance map (make sure there is enough statistics)
  Float_t nada = 0.;
  TH2F *trdAcc = (TH2F*)Proj3D(trd3D, 0x0, 1, nBinsZ, nada)->Clone();
  // get the bin limits from the Z axis of 3D histos
  Float_t *ptBinLimits = new Float_t[nBinsZ+1];
  for(Int_t i=1; i<=nBinsZ; i++) {
    ptBinLimits[i-1] = trd3D->GetZaxis()->GetBinLowEdge(i);
  }
  ptBinLimits[nBinsZ] = trd3D->GetZaxis()->GetBinUpEdge(nBinsZ);
  TH1F *efficiency = new TH1F("eff", "TRD-TPC matching efficiency", nBinsZ, ptBinLimits);
  // loop over Z bins
  for(Int_t i=1; i<=nBinsZ; i++) {
    Float_t tpcEntries = 0.0; Float_t trdEntries = 0.0;
    Proj3D(tpc3D, trdAcc, i, i, tpcEntries);
    Proj3D(trd3D, trdAcc, i, i, trdEntries);
    Float_t ratio = 0;
    if(tpcEntries>0) ratio = trdEntries/tpcEntries;
    Float_t error = 0;
    if(tpcEntries>0 && trdEntries>0) 
      error = ratio*TMath::Sqrt((1.0/tpcEntries)+(1.0/trdEntries));
    if(ratio>0.0) {
      efficiency->SetBinContent(i,ratio);
      efficiency->SetBinError(i,error);
    }
  }     // end loop over Z bins
  
  efficiency->SetLineColor((positives>0 ? 2 : 4));
  return efficiency;
}
