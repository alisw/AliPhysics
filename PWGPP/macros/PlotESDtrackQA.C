#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH3F.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include "TTree.h"
#include "TProfile.h"
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#endif

TH1D* ComputeMatchEff(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle);
TH1D* ComputeRatio(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle);
TH1D* ComputeFraction(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle);
void DrawDistribTrHyp(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TString pname, Bool_t showStat);
void DrawDistrib(TH1D* h1, TH1D* h2, TH1D* h3, Bool_t showStat);
void FillMeanAndRms(TH2F* h2d, TGraphErrors* gMean, TGraphErrors* gRms);
void InitFuncAndFit(TH1D* hm, TF1* fmass, Bool_t isK0s);


void PlotESDtrackQA(TString filename="AnalysisResults.root", TString suffix="QA", Int_t runNumber=-1){
  
  TFile* f=new TFile(filename.Data());
  TDirectoryFile* df=(TDirectoryFile*)f->Get("CheckESDTracks");
  TList* l=(TList*)df->Get(Form("clistCheckESDTracks%s",suffix.Data()));

  TH3F* hEtaPhiPtTPCsel=(TH3F*)l->FindObject("hEtaPhiPtTPCsel");
  TH3F* hEtaPhiPtTPCselITSref=(TH3F*)l->FindObject("hEtaPhiPtTPCselITSref");
  TH3F* hEtaPhiPtTPCselSPDany=(TH3F*)l->FindObject("hEtaPhiPtTPCselSPDany");
  TH3F* hEtaPhiPtPosChargeTPCsel=(TH3F*)l->FindObject("hEtaPhiPtPosChargeTPCsel");
  TH3F* hEtaPhiPtPosChargeTPCselITSref=(TH3F*)l->FindObject("hEtaPhiPtPosChargeTPCselITSref");
  TH3F* hEtaPhiPtPosChargeTPCselSPDany=(TH3F*)l->FindObject("hEtaPhiPtPosChargeTPCselSPDany");
  TH3F* hEtaPhiPtNegChargeTPCsel=(TH3F*)l->FindObject("hEtaPhiPtNegChargeTPCsel");
  TH3F* hEtaPhiPtNegChargeTPCselITSref=(TH3F*)l->FindObject("hEtaPhiPtNegChargeTPCselITSref");
  TH3F* hEtaPhiPtNegChargeTPCselSPDany=(TH3F*)l->FindObject("hEtaPhiPtNegChargeTPCselSPDany");


  TH1F* hev=(TH1F*)l->FindObject("hNEvents");
  Int_t etamin=hEtaPhiPtTPCsel->GetXaxis()->FindBin(-0.799);
  Int_t etamax=hEtaPhiPtTPCsel->GetXaxis()->FindBin(0.799);
  Int_t eta0p=hEtaPhiPtTPCsel->GetXaxis()->FindBin(0.01);
  Int_t eta0m=hEtaPhiPtTPCsel->GetXaxis()->FindBin(-0.01);
  printf("Bin %d Range %f %f\n",etamin,hEtaPhiPtTPCsel->GetXaxis()->GetBinLowEdge(etamin),hEtaPhiPtTPCsel->GetXaxis()->GetBinUpEdge(etamin));  printf("Bin %d Range %f %f\n",eta0m,hEtaPhiPtTPCsel->GetXaxis()->GetBinLowEdge(eta0m),hEtaPhiPtTPCsel->GetXaxis()->GetBinUpEdge(eta0m));
  printf("Bin %d Range %f %f\n",eta0p,hEtaPhiPtTPCsel->GetXaxis()->GetBinLowEdge(eta0p),hEtaPhiPtTPCsel->GetXaxis()->GetBinUpEdge(eta0p));
  printf("Bin %d Range %f %f\n",etamax,hEtaPhiPtTPCsel->GetXaxis()->GetBinLowEdge(etamax),hEtaPhiPtTPCsel->GetXaxis()->GetBinUpEdge(etamax));

  Int_t ptzero4=hEtaPhiPtTPCsel->GetZaxis()->FindBin(0.4001);
  Int_t ptzero7=hEtaPhiPtTPCsel->GetZaxis()->FindBin(0.6999);
  Int_t ptone=hEtaPhiPtTPCsel->GetZaxis()->FindBin(1.0001);
  Int_t ptten=hEtaPhiPtTPCsel->GetZaxis()->FindBin(9.999);

  TH1D* hPhiEtaNegTPCsel=hEtaPhiPtTPCsel->ProjectionY("hPhiEtaNegTPCsel",etamin,eta0m);
  TH1D* hPhiEtaPosTPCsel=hEtaPhiPtTPCsel->ProjectionY("hPhiEtaPosTPCsel",eta0p,etamax);
  TH1D* hPtEtaNegTPCsel=hEtaPhiPtTPCsel->ProjectionZ("hPtEtaNegTPCsel",etamin,eta0m);
  TH1D* hPtEtaPosTPCsel=hEtaPhiPtTPCsel->ProjectionZ("hPtEtaPosTPCsel",eta0p,etamax);

  TH1D* hPhiEtaNegTPCselLowPt=hEtaPhiPtTPCsel->ProjectionY("hPhiEtaNegTPCselLowPt",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselLowPt=hEtaPhiPtTPCsel->ProjectionY("hPhiEtaPosTPCselLowPt",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselHighPt=hEtaPhiPtTPCsel->ProjectionY("hPhiEtaNegTPCselHighPt",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselHighPt=hEtaPhiPtTPCsel->ProjectionY("hPhiEtaPosTPCselHighPt",eta0p,etamax,ptone,ptten);

  TH1D* hPhiEtaNegTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaNegTPCselITSref",etamin,eta0m);
  TH1D* hPhiEtaPosTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaPosTPCselITSref",eta0p,etamax);
  TH1D* hPtEtaNegTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionZ("hPtEtaNegTPCselITSref",etamin,eta0m);
  TH1D* hPtEtaPosTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionZ("hPtEtaPosTPCselITSref",eta0p,etamax);

  TH1D* hPhiEtaNegTPCselITSrefLowPt=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaNegTPCselITSrefLowPt",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselITSrefLowPt=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaPosTPCselITSrefLowPt",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselITSrefHighPt=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaNegTPCselITSrefHighPt",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselITSrefHighPt=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaPosTPCselITSrefHighPt",eta0p,etamax,ptone,ptten);

  TH1D* hPhiEtaNegTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaNegTPCselSPDany",etamin,eta0m);
  TH1D* hPhiEtaPosTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaPosTPCselSPDany",eta0p,etamax);
  TH1D* hPtEtaNegTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionZ("hPtEtaNegTPCselSPDany",etamin,eta0m);
  TH1D* hPtEtaPosTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionZ("hPtEtaPosTPCselSPDany",eta0p,etamax);

  TH1D* hPhiEtaNegTPCselSPDanyLowPt=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaNegTPCselSPDanyLowPt",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselSPDanyLowPt=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaPosTPCselSPDanyLowPt",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselSPDanyHighPt=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaNegTPCselSPDanyHighPt",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselSPDanyHighPt=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaPosTPCselSPDanyHighPt",eta0p,etamax,ptone,ptten);



  hPhiEtaNegTPCsel->SetMinimum(0);
  hPhiEtaPosTPCsel->SetMinimum(0);
  hPhiEtaNegTPCsel->SetTitle("#varphi tracks - #eta<0 - all p_{T}");
  hPhiEtaPosTPCsel->SetTitle("#varphi tracks - #eta>0 - all p_{T}");
  hPhiEtaNegTPCsel->GetXaxis()->SetTitle("#varphi (rad)");
  hPhiEtaPosTPCsel->GetXaxis()->SetTitle("#varphi (rad)");
  hPtEtaNegTPCsel->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtEtaPosTPCsel->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtEtaNegTPCsel->SetTitle("p_{T} tracks - #eta<0");
  hPtEtaPosTPCsel->SetTitle("p_{T} tracks - #eta>0");
  hPtEtaNegTPCselSPDany->SetTitle("p_{T} tracks - SPDany - #eta<0");
  hPtEtaPosTPCselSPDany->SetTitle("p_{T} tracks - SPDany - #eta>0");
  hPhiEtaNegTPCselLowPt->SetMinimum(0);
  hPhiEtaPosTPCselLowPt->SetMinimum(0);
  hPhiEtaNegTPCselLowPt->SetTitle("#varphi tracks - #eta<0 - 0.4<p_{T}<0.7 GeV/c");
  hPhiEtaPosTPCselLowPt->SetTitle("#varphi tracks - #eta>0 - 0.4<p_{T}<0.7 GeV/c");
  hPhiEtaNegTPCselLowPt->GetXaxis()->SetTitle("#varphi (rad)");
  hPhiEtaPosTPCselLowPt->GetXaxis()->SetTitle("#varphi (rad)");
  hPhiEtaNegTPCselHighPt->SetMinimum(0);
  hPhiEtaPosTPCselHighPt->SetMinimum(0);
  hPhiEtaNegTPCselHighPt->SetTitle("#varphi tracks - #eta<0 - 1<p_{T}<10 GeV/c");
  hPhiEtaPosTPCselHighPt->SetTitle("#varphi tracks - #eta>0 - 1<p_{T}<10 GeV/c");
  hPhiEtaNegTPCselHighPt->GetXaxis()->SetTitle("#varphi (rad)");
  hPhiEtaPosTPCselHighPt->GetXaxis()->SetTitle("#varphi (rad)");

  hPtEtaNegTPCsel->Sumw2();
  hPtEtaNegTPCselITSref->Sumw2();
  hPtEtaNegTPCselSPDany->Sumw2();
  hPtEtaPosTPCsel->Sumw2();
  hPtEtaPosTPCselITSref->Sumw2();
  hPtEtaPosTPCselSPDany->Sumw2();
  hPhiEtaNegTPCsel->Sumw2();
  hPhiEtaNegTPCselITSref->Sumw2();
  hPhiEtaNegTPCselSPDany->Sumw2();
  hPhiEtaPosTPCsel->Sumw2();
  hPhiEtaPosTPCselITSref->Sumw2();
  hPhiEtaPosTPCselSPDany->Sumw2();
  hPhiEtaNegTPCselLowPt->Sumw2();
  hPhiEtaNegTPCselITSrefLowPt->Sumw2();
  hPhiEtaNegTPCselSPDanyLowPt->Sumw2();
  hPhiEtaPosTPCselLowPt->Sumw2();
  hPhiEtaPosTPCselITSrefLowPt->Sumw2();
  hPhiEtaPosTPCselSPDanyLowPt->Sumw2();
  hPhiEtaNegTPCselHighPt->Sumw2();
  hPhiEtaNegTPCselITSrefHighPt->Sumw2();
  hPhiEtaNegTPCselSPDanyHighPt->Sumw2();
  hPhiEtaPosTPCselHighPt->Sumw2();
  hPhiEtaPosTPCselITSrefHighPt->Sumw2();
  hPhiEtaPosTPCselSPDanyHighPt->Sumw2();

  if(hEtaPhiPtPosChargeTPCsel){
    TH1D* hPtEtaNegPosChargeTPCsel=hEtaPhiPtPosChargeTPCsel->ProjectionZ("hPtEtaNegPosChargeTPCsel",etamin,eta0m);
    TH1D* hPtEtaPosPosChargeTPCsel=hEtaPhiPtPosChargeTPCsel->ProjectionZ("hPtEtaPosPosChargeTPCsel",eta0p,etamax);
    TH1D* hPtEtaNegNegChargeTPCsel=hEtaPhiPtNegChargeTPCsel->ProjectionZ("hPtEtaNegNegChargeTPCsel",etamin,eta0m);
    TH1D* hPtEtaPosNegChargeTPCsel=hEtaPhiPtNegChargeTPCsel->ProjectionZ("hPtEtaPosNegChargeTPCsel",eta0p,etamax);
    TH1D* hPtEtaNegPosChargeTPCselITSref=hEtaPhiPtPosChargeTPCselITSref->ProjectionZ("hPtEtaNegPosChargeTPCselITSref",etamin,eta0m);
    TH1D* hPtEtaPosPosChargeTPCselITSref=hEtaPhiPtPosChargeTPCselITSref->ProjectionZ("hPtEtaPosPosChargeTPCselITSref",eta0p,etamax);
    TH1D* hPtEtaNegNegChargeTPCselITSref=hEtaPhiPtNegChargeTPCselITSref->ProjectionZ("hPtEtaNegNegChargeTPCselITSref",etamin,eta0m);
    TH1D* hPtEtaPosNegChargeTPCselITSref=hEtaPhiPtNegChargeTPCselITSref->ProjectionZ("hPtEtaPosNegChargeTPCselITSref",eta0p,etamax);
    TH1D* hPtEtaNegPosChargeTPCselSPDany=hEtaPhiPtPosChargeTPCselSPDany->ProjectionZ("hPtEtaNegPosChargeTPCselSPDany",etamin,eta0m);
    TH1D* hPtEtaPosPosChargeTPCselSPDany=hEtaPhiPtPosChargeTPCselSPDany->ProjectionZ("hPtEtaPosPosChargeTPCselSPDany",eta0p,etamax);
    TH1D* hPtEtaNegNegChargeTPCselSPDany=hEtaPhiPtNegChargeTPCselSPDany->ProjectionZ("hPtEtaNegNegChargeTPCselSPDany",etamin,eta0m);
    TH1D* hPtEtaPosNegChargeTPCselSPDany=hEtaPhiPtNegChargeTPCselSPDany->ProjectionZ("hPtEtaPosNegChargeTPCselSPDany",eta0p,etamax);
    TH1D* hRatioPosNegEtaNegTPCsel=ComputeRatio(hPtEtaNegPosChargeTPCsel,hPtEtaNegNegChargeTPCsel,"hRatioPosNegEtaNegTPCsel",kGray+1,21,"p_{T} (GeV/c)");
    hRatioPosNegEtaNegTPCsel->GetYaxis()->SetTitle("Positive Charge / Negative Charge");
    TH1D* hRatioPosNegEtaPosTPCsel=ComputeRatio(hPtEtaPosPosChargeTPCsel,hPtEtaPosNegChargeTPCsel,"hRatioPosNegEtaPosTPCsel",kGray+1,21,"p_{T} (GeV/c)");
    hRatioPosNegEtaPosTPCsel->GetYaxis()->SetTitle("Positive Charge / Negative Charge");

    TH1D* hRatioPosNegEtaNegTPCselSPDany=ComputeRatio(hPtEtaNegPosChargeTPCselSPDany,hPtEtaNegNegChargeTPCselSPDany,"hRatioPosNegEtaNegTPCsel",kGray+1,21,"p_{T} (GeV/c)");
    hRatioPosNegEtaNegTPCselSPDany->GetYaxis()->SetTitle("Positive Charge / Negative Charge");
    TH1D* hRatioPosNegEtaPosTPCselSPDany=ComputeRatio(hPtEtaPosPosChargeTPCselSPDany,hPtEtaPosNegChargeTPCselSPDany,"hRatioPosNegEtaPosTPCsel",kGray+1,21,"p_{T} (GeV/c)");
    hRatioPosNegEtaPosTPCselSPDany->GetYaxis()->SetTitle("Positive Charge / Negative Charge");

    TCanvas* cdist=new TCanvas("cdist","Pt Distrib TPD sel",900,900);
    cdist->Divide(2,2);
    cdist->cd(1);
    gPad->SetLogy();
    DrawDistrib(hPtEtaNegTPCsel,hPtEtaNegPosChargeTPCsel,hPtEtaNegNegChargeTPCsel,kTRUE);
    TLegend* legd=new TLegend(0.46,0.66,0.79,0.89);
    legd->AddEntry(hPtEtaNegTPCsel,"All charges","L")->SetTextColor(hPtEtaNegTPCsel->GetLineColor());
    legd->AddEntry(hPtEtaNegPosChargeTPCsel,"Positive","L")->SetTextColor(hPtEtaNegTPCselITSref->GetLineColor());
    legd->AddEntry(hPtEtaNegNegChargeTPCsel,"Negative","L")->SetTextColor(hPtEtaNegTPCselSPDany->GetLineColor());
    legd->Draw();
    cdist->cd(2);
    gPad->SetLogy();
    DrawDistrib(hPtEtaPosTPCsel,hPtEtaPosPosChargeTPCsel,hPtEtaPosNegChargeTPCsel,kTRUE);
    cdist->cd(3);
    hRatioPosNegEtaNegTPCsel->Draw();
    cdist->cd(4);
    hRatioPosNegEtaPosTPCsel->Draw();
    cdist->SaveAs("TracksPtDistrib-TPCsel.png");

    TCanvas* cdists=new TCanvas("cdists","Pt Distrib SPDany",900,900);
    cdists->Divide(2,2);
    cdists->cd(1);
    gPad->SetLogy();
    DrawDistrib(hPtEtaNegTPCselSPDany,hPtEtaNegPosChargeTPCselSPDany,hPtEtaNegNegChargeTPCselSPDany,kTRUE);
    legd->Draw();
    cdists->cd(2);
    gPad->SetLogy();
    DrawDistrib(hPtEtaPosTPCselSPDany,hPtEtaPosPosChargeTPCselSPDany,hPtEtaPosNegChargeTPCselSPDany,kTRUE);
    cdists->cd(3);
    hRatioPosNegEtaNegTPCsel->Draw();
    cdists->cd(4);
    hRatioPosNegEtaPosTPCsel->Draw();
    cdists->SaveAs("TracksPtDistrib-TPCselSPDany.png");
  }else{
    TCanvas* cdist=new TCanvas("cdist","Pt+Phi Distrib",900,900);
    cdist->Divide(2,2);
    cdist->cd(1);
    gPad->SetLogy();
    DrawDistrib(hPtEtaNegTPCsel,hPtEtaNegTPCselITSref,hPtEtaNegTPCselSPDany,kTRUE);
    TLegend* legd=new TLegend(0.16,0.16,0.5,0.4);
    legd->AddEntry(hPtEtaNegTPCsel,"TPC only cuts","L")->SetTextColor(hPtEtaNegTPCsel->GetLineColor());
    legd->AddEntry(hPtEtaNegTPCselITSref,"ITSrefit","L")->SetTextColor(hPtEtaNegTPCselITSref->GetLineColor());
    legd->AddEntry(hPtEtaNegTPCselSPDany,"SPD any","L")->SetTextColor(hPtEtaNegTPCselSPDany->GetLineColor());
    legd->Draw();
    cdist->cd(2);
    gPad->SetLogy();
    DrawDistrib(hPtEtaPosTPCsel,hPtEtaPosTPCselITSref,hPtEtaPosTPCselSPDany,kTRUE);
    cdist->cd(3);
    DrawDistrib(hPhiEtaNegTPCsel,hPhiEtaNegTPCselITSref,hPhiEtaNegTPCselSPDany,kFALSE);
    legd->Draw();
    cdist->cd(4);
    DrawDistrib(hPhiEtaPosTPCsel,hPhiEtaPosTPCselITSref,hPhiEtaPosTPCselSPDany,kFALSE);
    cdist->SaveAs("TracksPtPhiDistrib.png");
  }

  TCanvas* cdist2=new TCanvas("cdist2","Phi Distrib",900,900);
  cdist2->Divide(2,2);
  cdist2->cd(1);
  DrawDistrib(hPhiEtaNegTPCselLowPt,hPhiEtaNegTPCselITSrefLowPt,hPhiEtaNegTPCselSPDanyLowPt,kFALSE);
  TLegend* legdf=new TLegend(0.16,0.16,0.5,0.4);
  legdf->AddEntry(hPhiEtaNegTPCselLowPt,"TPC only cuts","L")->SetTextColor(hPhiEtaNegTPCselLowPt->GetLineColor());
  legdf->AddEntry(hPhiEtaNegTPCselITSrefLowPt,"ITSrefit","L")->SetTextColor(hPhiEtaNegTPCselITSrefLowPt->GetLineColor());
  legdf->AddEntry(hPhiEtaNegTPCselSPDanyLowPt,"SPD any","L")->SetTextColor(hPhiEtaNegTPCselSPDanyLowPt->GetLineColor());
  legdf->Draw();
  cdist2->cd(2);
  DrawDistrib(hPhiEtaPosTPCselLowPt,hPhiEtaPosTPCselITSrefLowPt,hPhiEtaPosTPCselSPDanyLowPt,kFALSE);
  cdist2->cd(3);
  DrawDistrib(hPhiEtaNegTPCselHighPt,hPhiEtaNegTPCselITSrefHighPt,hPhiEtaNegTPCselSPDanyHighPt,kFALSE);
  cdist2->cd(4);
  DrawDistrib(hPhiEtaPosTPCselHighPt,hPhiEtaPosTPCselITSrefHighPt,hPhiEtaPosTPCselSPDanyHighPt,kFALSE);
  cdist2->SaveAs("TracksPhiDistrib-PtBins.png");


  
  const Int_t checkSpecies=2;
  TString pNames[checkSpecies]={"Pion","Proton"};
  TString trType[3]={"TPCsel","TPCselITSref","TPCselSPDany"};
  TH3F* hEtaPhiPtGoodHyp[3][checkSpecies];
  TH3F* hEtaPhiPtBadHyp[3][checkSpecies];
  TH1D* hPtGoodHyp[3][checkSpecies];
  TH1D* hPtEtaNegGoodHyp[3][checkSpecies];
  TH1D* hPtEtaPosGoodHyp[3][checkSpecies];
  TH1D* hPhiGoodHyp[3][checkSpecies];
  TH1D* hPhiEtaNegGoodHyp[3][checkSpecies];
  TH1D* hPhiEtaPosGoodHyp[3][checkSpecies];
  TH1D* hPtBadHyp[3][checkSpecies];
  TH1D* hPtEtaNegBadHyp[3][checkSpecies];
  TH1D* hPtEtaPosBadHyp[3][checkSpecies];
  TH1D* hPhiBadHyp[3][checkSpecies];
  TH1D* hPhiEtaNegBadHyp[3][checkSpecies];
  TH1D* hPhiEtaPosBadHyp[3][checkSpecies];

  TH1D* hRatioBadGoodVsPtEtaNeg[3][checkSpecies];
  TH1D* hRatioBadGoodVsPtEtaPos[3][checkSpecies];


  TH1D* hMatchEffGoodVsPtEtaNeg[2][checkSpecies];
  TH1D* hMatchEffGoodVsPtEtaPos[2][checkSpecies];
  TH1D* hMatchEffBadVsPtEtaNeg[2][checkSpecies];
  TH1D* hMatchEffBadVsPtEtaPos[2][checkSpecies];
  TH1D* hMatchEffGoodVsPhiEtaNeg[2][checkSpecies];
  TH1D* hMatchEffGoodVsPhiEtaPos[2][checkSpecies];
  TH1D* hMatchEffBadVsPhiEtaNeg[2][checkSpecies];
  TH1D* hMatchEffBadVsPhiEtaPos[2][checkSpecies];

  Double_t mCol[3]={1,kRed+1,kBlue+1};
  Double_t mStyle[3]={22,23,26};

  Double_t mColG[2]={1,kBlue-7};
  Double_t mStyleG[2]={20,33};
  Double_t mStyleB[2]={24,27};

  for(Int_t iSp=0; iSp<checkSpecies; iSp++){
    for(Int_t jTy=0; jTy<3; jTy++){    
      hEtaPhiPtGoodHyp[jTy][iSp]=(TH3F*)l->FindObject(Form("hEtaPhiPtGoodHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()));
      hPtGoodHyp[jTy][iSp]=hEtaPhiPtGoodHyp[jTy][iSp]->ProjectionZ(Form("hPtGoodHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),etamin,etamax);
      hPtEtaNegGoodHyp[jTy][iSp]=hEtaPhiPtGoodHyp[jTy][iSp]->ProjectionZ(Form("hPtEtaNegGoodHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),etamin,eta0m);
      hPtEtaPosGoodHyp[jTy][iSp]=hEtaPhiPtGoodHyp[jTy][iSp]->ProjectionZ(Form("hPtEtaPosGoodHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),eta0p,etamax);
      hPhiGoodHyp[jTy][iSp]=hEtaPhiPtGoodHyp[jTy][iSp]->ProjectionY(Form("hPhiGoodHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),etamin,etamax);
      hPhiEtaNegGoodHyp[jTy][iSp]=hEtaPhiPtGoodHyp[jTy][iSp]->ProjectionY(Form("hPhiEtaNegGoodHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),etamin,eta0m);
      hPhiEtaPosGoodHyp[jTy][iSp]=hEtaPhiPtGoodHyp[jTy][iSp]->ProjectionY(Form("hPhiEtaPosGoodHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),eta0p,etamax);

      hEtaPhiPtBadHyp[jTy][iSp]=(TH3F*)l->FindObject(Form("hEtaPhiPtBadHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()));
      hPtBadHyp[jTy][iSp]=hEtaPhiPtBadHyp[jTy][iSp]->ProjectionZ(Form("hPtBadHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),etamin,etamax);
      hPtEtaNegBadHyp[jTy][iSp]=hEtaPhiPtBadHyp[jTy][iSp]->ProjectionZ(Form("hPtEtaNegBadHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),etamin,eta0m);
      hPtEtaPosBadHyp[jTy][iSp]=hEtaPhiPtBadHyp[jTy][iSp]->ProjectionZ(Form("hPtEtaPosBadHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),eta0p,etamax);
      hPhiBadHyp[jTy][iSp]=hEtaPhiPtBadHyp[jTy][iSp]->ProjectionY(Form("hPhiBadHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),etamin,etamax);
      hPhiEtaNegBadHyp[jTy][iSp]=hEtaPhiPtBadHyp[jTy][iSp]->ProjectionY(Form("hPhiEtaNegBadHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),etamin,eta0m);
      hPhiEtaPosBadHyp[jTy][iSp]=hEtaPhiPtBadHyp[jTy][iSp]->ProjectionY(Form("hPhiEtaPosBadHyp%s%s",pNames[iSp].Data(),trType[jTy].Data()),eta0p,etamax);

      hPtEtaNegGoodHyp[jTy][iSp]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hPtEtaNegGoodHyp[jTy][iSp]->SetTitle("p_{T} tracks - #eta<0");
      hPtEtaPosGoodHyp[jTy][iSp]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hPtEtaPosGoodHyp[jTy][iSp]->SetTitle("p_{T} tracks - #eta>0");
      hPhiEtaNegGoodHyp[jTy][iSp]->GetXaxis()->SetTitle("#varphi (rad)");
      hPhiEtaNegGoodHyp[jTy][iSp]->SetTitle("#varphi tracks - #eta<0");
      hPhiEtaPosGoodHyp[jTy][iSp]->GetXaxis()->SetTitle("#varphi (rad)");
      hPhiEtaPosGoodHyp[jTy][iSp]->SetTitle("#varphi tracks - #eta>0");

      hRatioBadGoodVsPtEtaNeg[jTy][iSp]=ComputeFraction(hPtEtaNegBadHyp[jTy][iSp],hPtEtaNegGoodHyp[jTy][iSp],Form("hRatio%sBadGoodVsPtEtaNeg",pNames[iSp].Data()),mCol[jTy],mStyle[jTy],"p_{T} (GeV/c)");
      hRatioBadGoodVsPtEtaPos[jTy][iSp]=ComputeFraction(hPtEtaPosBadHyp[jTy][iSp],hPtEtaPosGoodHyp[jTy][iSp],Form("hRatio%sBadGoodVsPtEtaPos",pNames[iSp].Data()),mCol[jTy],mStyle[jTy],"p_{T} (GeV/c)");
    }

    for(Int_t jTy=1; jTy<3; jTy++){
      hMatchEffGoodVsPtEtaNeg[jTy-1][iSp]=ComputeMatchEff(hPtEtaNegGoodHyp[jTy][iSp],hPtEtaNegGoodHyp[0][iSp],Form("hMatchEff%sGood%sVsPtEtaNeg",pNames[iSp].Data(),trType[jTy].Data()),mColG[jTy-1],mStyleG[jTy-1],"p_{T} (GeV/c)");
      hMatchEffGoodVsPtEtaPos[jTy-1][iSp]=ComputeMatchEff(hPtEtaPosGoodHyp[jTy][iSp],hPtEtaPosGoodHyp[0][iSp],Form("hMatchEff%sGood%sVsPtEtaPos",pNames[iSp].Data(),trType[jTy].Data()),mColG[jTy-1],mStyleG[jTy-1],"p_{T} (GeV/c)");
      hMatchEffGoodVsPhiEtaNeg[jTy-1][iSp]=ComputeMatchEff(hPhiEtaNegGoodHyp[jTy][iSp],hPhiEtaNegGoodHyp[0][iSp],Form("hMatchEff%sGood%sVsPhiEtaNeg",pNames[iSp].Data(),trType[jTy].Data()),mColG[jTy-1],mStyleG[jTy-1],"#varphi (rad)");
      hMatchEffGoodVsPhiEtaPos[jTy-1][iSp]=ComputeMatchEff(hPhiEtaPosGoodHyp[jTy][iSp],hPhiEtaPosGoodHyp[0][iSp],Form("hMatchEff%sGood%sVsPhiEtaPos",pNames[iSp].Data(),trType[jTy].Data()),mColG[jTy-1],mStyleG[jTy-1],"#varphi (rad)");
      hMatchEffBadVsPtEtaNeg[jTy-1][iSp]=ComputeMatchEff(hPtEtaNegBadHyp[jTy][iSp],hPtEtaNegBadHyp[0][iSp],Form("hMatchEff%sBad%sVsPtEtaNeg",pNames[iSp].Data(),trType[jTy].Data()),mColG[jTy-1],mStyleB[jTy-1],"p_{T} (GeV/c)");
      hMatchEffBadVsPtEtaPos[jTy-1][iSp]=ComputeMatchEff(hPtEtaPosBadHyp[jTy][iSp],hPtEtaPosBadHyp[0][iSp],Form("hMatchEff%sBad%sVsPtEtaPos",pNames[iSp].Data(),trType[jTy].Data()),mColG[jTy-1],mStyleB[jTy-1],"p_{T} (GeV/c)");
      hMatchEffBadVsPhiEtaNeg[jTy-1][iSp]=ComputeMatchEff(hPhiEtaNegBadHyp[jTy][iSp],hPhiEtaNegBadHyp[0][iSp],Form("hMatchEff%sBad%sVsPhiEtaNeg",pNames[iSp].Data(),trType[jTy].Data()),mColG[jTy-1],mStyleB[jTy-1],"#varphi (rad)");
      hMatchEffBadVsPhiEtaPos[jTy-1][iSp]=ComputeMatchEff(hPhiEtaPosBadHyp[jTy][iSp],hPhiEtaPosBadHyp[0][iSp],Form("hMatchEff%sBad%sVsPhiEtaPos",pNames[iSp].Data(),trType[jTy].Data()),mColG[jTy-1],mStyleB[jTy-1],"#varphi (rad)");

      hMatchEffGoodVsPtEtaNeg[jTy-1][iSp]->SetTitle("#eta<0");
      hMatchEffGoodVsPtEtaPos[jTy-1][iSp]->SetTitle("#eta<0");
      hMatchEffGoodVsPhiEtaNeg[jTy-1][iSp]->SetTitle("#eta<0");
      hMatchEffGoodVsPhiEtaPos[jTy-1][iSp]->SetTitle("#eta<0");
      hMatchEffBadVsPtEtaNeg[jTy-1][iSp]->SetTitle("#eta<0");
      hMatchEffBadVsPtEtaPos[jTy-1][iSp]->SetTitle("#eta<0");
      hMatchEffBadVsPhiEtaNeg[jTy-1][iSp]->SetTitle("#eta<0");
      hMatchEffBadVsPhiEtaPos[jTy-1][iSp]->SetTitle("#eta<0");
    }

    TCanvas* cdistp=new TCanvas(Form("cdistp%d",iSp),Form("%s Distrib",pNames[iSp].Data()),900,900);
    cdistp->Divide(2,2);
    cdistp->cd(1);
    gPad->SetLogy();
    DrawDistribTrHyp(hPtEtaNegGoodHyp[0][iSp],hPtEtaNegBadHyp[0][iSp],hPtEtaNegGoodHyp[1][iSp],hPtEtaNegBadHyp[1][iSp],pNames[iSp].Data(),kTRUE);
    cdistp->cd(2);
    gPad->SetLogy();
    DrawDistribTrHyp(hPtEtaPosGoodHyp[0][iSp],hPtEtaPosBadHyp[0][iSp],hPtEtaPosGoodHyp[1][iSp],hPtEtaPosBadHyp[1][iSp],pNames[iSp].Data(),kTRUE);
    cdistp->cd(3);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    hRatioBadGoodVsPtEtaNeg[0][iSp]->SetTitle("#eta<0");
    hRatioBadGoodVsPtEtaNeg[0][iSp]->GetXaxis()->SetRangeUser(0.,10.);
    hRatioBadGoodVsPtEtaNeg[0][iSp]->Draw();
    hRatioBadGoodVsPtEtaNeg[1][iSp]->Draw("same");
    hRatioBadGoodVsPtEtaNeg[2][iSp]->Draw("same");
    TLegend* legr=new TLegend(0.4,0.2,0.83,0.4);
    legr->AddEntry(hRatioBadGoodVsPtEtaNeg[0][iSp],"TPC cuts","P");
    legr->AddEntry(hRatioBadGoodVsPtEtaNeg[1][iSp],"TPC cuts+ITS refit","P");
    legr->AddEntry(hRatioBadGoodVsPtEtaNeg[2][iSp],"TPC cuts+SPD any","P");
    legr->Draw();
    cdistp->cd(4);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    hRatioBadGoodVsPtEtaPos[0][iSp]->SetTitle("#eta>0");
    hRatioBadGoodVsPtEtaPos[0][iSp]->GetXaxis()->SetRangeUser(0.,10.);
    hRatioBadGoodVsPtEtaPos[0][iSp]->Draw();
    hRatioBadGoodVsPtEtaPos[1][iSp]->Draw("same");
    hRatioBadGoodVsPtEtaPos[2][iSp]->Draw("same");
 
    

    TCanvas* cmep=new TCanvas(Form("cmep%d",iSp),Form("MatchEff %s",pNames[iSp].Data()),900,900);
    cmep->Divide(2,2);
    cmep->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    hMatchEffGoodVsPtEtaNeg[0][iSp]->GetXaxis()->SetRangeUser(0.,10.);
    hMatchEffGoodVsPtEtaNeg[0][iSp]->Draw("PE");
    hMatchEffGoodVsPtEtaNeg[1][iSp]->Draw("samepe");
    hMatchEffBadVsPtEtaNeg[0][iSp]->Draw("samepe");
    hMatchEffBadVsPtEtaNeg[1][iSp]->Draw("samepe");
    TLegend* legm=new TLegend(0.4,0.16,0.84,0.4);
    legm->AddEntry(hMatchEffGoodVsPtEtaNeg[0][iSp],Form("%s - Good Hyp - ITSrefit",pNames[iSp].Data()),"P");
    legm->AddEntry(hMatchEffBadVsPtEtaNeg[0][iSp],Form("%s - Bad Hyp - ITSrefit",pNames[iSp].Data()),"P");
    legm->AddEntry(hMatchEffGoodVsPtEtaNeg[1][iSp],Form("%s - Good Hyp - SPDany",pNames[iSp].Data()),"P");
    legm->AddEntry(hMatchEffBadVsPtEtaNeg[1][iSp],Form("%s - Bad Hyp - SPDany",pNames[iSp].Data()),"P");
    legm->Draw();
    cmep->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    hMatchEffGoodVsPtEtaPos[0][iSp]->GetXaxis()->SetRangeUser(0.,10.);
    hMatchEffGoodVsPtEtaPos[0][iSp]->Draw("PE");
    hMatchEffGoodVsPtEtaPos[1][iSp]->Draw("samepe");
    hMatchEffBadVsPtEtaPos[0][iSp]->Draw("samepe");
    hMatchEffBadVsPtEtaPos[1][iSp]->Draw("samepe");
    cmep->cd(3);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    hMatchEffGoodVsPhiEtaNeg[0][iSp]->Draw("PE");
    hMatchEffGoodVsPhiEtaNeg[1][iSp]->Draw("samepe");
    hMatchEffBadVsPhiEtaNeg[0][iSp]->Draw("samepe");
    hMatchEffBadVsPhiEtaNeg[1][iSp]->Draw("samepe");
    cmep->cd(4);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    hMatchEffGoodVsPhiEtaPos[0][iSp]->Draw("PE");
    hMatchEffGoodVsPhiEtaPos[1][iSp]->Draw("samepe");
    hMatchEffBadVsPhiEtaPos[0][iSp]->Draw("samepe");
    hMatchEffBadVsPhiEtaPos[1][iSp]->Draw("samepe");
    cmep->SaveAs(Form("MatchingEfficiency-%s.png",pNames[iSp].Data()));
  }

  TH3F* hptresTPC3d=(TH3F*)l->FindObject("hTPCsig1ptPerClusPhiPtTPCsel");
  TH3F* hptresITS3d=(TH3F*)l->FindObject("hTPCsig1ptPerClusPhiPtTPCselITSref");
  TH3F* hptresSPD3d=(TH3F*)l->FindObject("hTPCsig1ptPerClusPhiPtTPCselSPDany");
  TH2D* hptresTPC=(TH2D*)hptresTPC3d->Project3D("xy");
  TH2D* hptresITS=(TH2D*)hptresITS3d->Project3D("xy");
  TH2D* hptresSPD=(TH2D*)hptresSPD3d->Project3D("xy");
  TProfile* pptresTPC=hptresTPC->ProfileX("pptresTPC");
  TProfile* pptresITS=hptresITS->ProfileX("pptresITS");
  TProfile* pptresSPD=hptresSPD->ProfileX("pptresSPD");
  pptresTPC->GetYaxis()->SetTitle(hptresTPC3d->GetXaxis()->GetTitle());
  pptresTPC->GetYaxis()->SetTitleOffset(1.5);
  pptresTPC->SetStats(0);
  pptresTPC->SetTitle("p_{T} resolution from cov. matrix");
  pptresTPC->SetLineColor(1);
  pptresTPC->SetMarkerStyle(20);
  pptresTPC->SetMarkerColor(1);
  pptresITS->SetLineColor(kRed+1);
  pptresITS->SetMarkerStyle(25);
  pptresITS->SetMarkerColor(kRed+1);
  pptresSPD->SetLineColor(kBlue-7);
  pptresSPD->SetMarkerStyle(33);
  pptresSPD->SetMarkerColor(kBlue-7);


  TCanvas* cptcm=new TCanvas("cptcm","Pt resol",800,800);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);

  pptresTPC->SetMinimum(0);
  pptresTPC->SetMaximum(0.08);
  pptresTPC->Draw();
  pptresITS->Draw("same");
  pptresSPD->Draw("same");
  TLegend* leg=new TLegend(0.17,0.7,0.4,0.87);
  leg->AddEntry(pptresTPC,"TPC only","P");
  leg->AddEntry(pptresITS,"ITSrefit","P");
  leg->AddEntry(pptresSPD,"SPD any","P");
  leg->Draw();
  cptcm->SaveAs("PtResolCovMat.png");

  TH2F* hPtResidVsPtTPCselITSrefPiong=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefGoodHyppi");
  TH2F* hPtResidVsPtTPCselITSrefKaong=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefGoodHypK");
  TH2F* hPtResidVsPtTPCselITSrefProtong=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefGoodHypp");
  TH2F* hPtResidVsPtTPCselITSrefPionb=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefBadHyppi");
  TH2F* hPtResidVsPtTPCselITSrefKaonb=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefBadHypK");
  TH2F* hPtResidVsPtTPCselITSrefProtonb=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefBadHypp");
  TH2F* hPtResidVsPtTPCselITSrefPion=(TH2F*)hPtResidVsPtTPCselITSrefPiong->Clone("hPtResidVsPtTPCselITSrefPion");
  hPtResidVsPtTPCselITSrefPion->Add(hPtResidVsPtTPCselITSrefPionb);
  TH2F* hPtResidVsPtTPCselITSrefKaon=(TH2F*)hPtResidVsPtTPCselITSrefKaong->Clone("hPtResidVsPtTPCselITSrefKaon");
  hPtResidVsPtTPCselITSrefKaon->Add(hPtResidVsPtTPCselITSrefKaonb);
  TH2F* hPtResidVsPtTPCselITSrefProton=(TH2F*)hPtResidVsPtTPCselITSrefProtong->Clone("hPtResidVsPtTPCselITSrefProton");
  hPtResidVsPtTPCselITSrefProton->Add(hPtResidVsPtTPCselITSrefProtonb);

  if(hPtResidVsPtTPCselITSrefPion){
    hPtResidVsPtTPCselITSrefPion->SetStats(0);
    hPtResidVsPtTPCselITSrefKaon->SetStats(0);
    hPtResidVsPtTPCselITSrefProton->SetStats(0);
    hPtResidVsPtTPCselITSrefPion->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefKaon->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefProton->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefPion->GetXaxis()->SetTitleOffset(1.1);
    hPtResidVsPtTPCselITSrefKaon->GetXaxis()->SetTitleOffset(1.1);
    hPtResidVsPtTPCselITSrefProton->GetXaxis()->SetTitleOffset(1.1);
  }

  TH2F* hOneOverPtResidVsPtTPCselITSrefPiong=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefGoodHyppi");
  TH2F* hOneOverPtResidVsPtTPCselITSrefKaong=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefGoodHypK");
  TH2F* hOneOverPtResidVsPtTPCselITSrefProtong=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefGoodHypp");
  TH2F* hOneOverPtResidVsPtTPCselITSrefPionb=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefBadHyppi");
  TH2F* hOneOverPtResidVsPtTPCselITSrefKaonb=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefBadHypK");
  TH2F* hOneOverPtResidVsPtTPCselITSrefProtonb=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefBadHypp");
  TH2F* hOneOverPtResidVsPtTPCselITSrefPion=(TH2F*)hOneOverPtResidVsPtTPCselITSrefPiong->Clone("hOneOverPtResidVsPtTPCselITSrefPion");
  hOneOverPtResidVsPtTPCselITSrefPion->Add(hOneOverPtResidVsPtTPCselITSrefPionb);
  TH2F* hOneOverPtResidVsPtTPCselITSrefKaon=(TH2F*)hOneOverPtResidVsPtTPCselITSrefKaong->Clone("hOneOverPtResidVsPtTPCselITSrefKaon");
  hOneOverPtResidVsPtTPCselITSrefKaon->Add(hOneOverPtResidVsPtTPCselITSrefKaonb);
  TH2F* hOneOverPtResidVsPtTPCselITSrefProton=(TH2F*)hOneOverPtResidVsPtTPCselITSrefProtong->Clone("hOneOverPtResidVsPtTPCselITSrefProton");
  hOneOverPtResidVsPtTPCselITSrefProton->Add(hOneOverPtResidVsPtTPCselITSrefProtonb);

  if(hOneOverPtResidVsPtTPCselITSrefPion){
    hOneOverPtResidVsPtTPCselITSrefPion->SetStats(0);
    hOneOverPtResidVsPtTPCselITSrefKaon->SetStats(0);
    hOneOverPtResidVsPtTPCselITSrefProton->SetStats(0);
    hOneOverPtResidVsPtTPCselITSrefPion->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselITSrefKaon->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselITSrefProton->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselITSrefPion->GetXaxis()->SetTitleOffset(1.1);
    hOneOverPtResidVsPtTPCselITSrefKaon->GetXaxis()->SetTitleOffset(1.1);
    hOneOverPtResidVsPtTPCselITSrefProton->GetXaxis()->SetTitleOffset(1.1);
  }

  TGraphErrors* gMeanPi=new TGraphErrors(0);
  TGraphErrors* gMeanK=new TGraphErrors(0);
  TGraphErrors* gMeanProt=new TGraphErrors(0);
  TGraphErrors* gRmsPi=new TGraphErrors(0);
  TGraphErrors* gRmsK=new TGraphErrors(0);
  TGraphErrors* gRmsProt=new TGraphErrors(0);
  TGraphErrors* gDum=new TGraphErrors(0);
  TGraphErrors* gRelPi=new TGraphErrors(0);
  TGraphErrors* gRelK=new TGraphErrors(0);
  TGraphErrors* gRelProt=new TGraphErrors(0);
  gMeanPi->SetName("gMeanProt");
  gMeanK->SetName("gMeanProt");
  gMeanProt->SetName("gMeanProt");
  gMeanPi->SetTitle("");
  gMeanK->SetTitle("");
  gMeanProt->SetTitle("");
  gRmsPi->SetName("gRmsProt");
  gRmsK->SetName("gRmsProt");
  gRmsProt->SetName("gRmsProt");
  gRmsPi->SetTitle("");
  gRmsK->SetTitle("");
  gRmsProt->SetTitle("");
  gRelPi->SetName("gRelProt");
  gRelK->SetName("gRelProt");
  gRelProt->SetName("gRelProt");
  gRelPi->SetTitle("");
  gRelK->SetTitle("");
  gRelProt->SetTitle("");

  Bool_t okRes=kFALSE;
  Bool_t okOneOverRes=kFALSE;
  if(hPtResidVsPtTPCselITSrefPion && hPtResidVsPtTPCselITSrefPion->Integral()>0){
    FillMeanAndRms(hPtResidVsPtTPCselITSrefPion,gMeanPi,gRmsPi);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefKaon,gMeanK,gRmsK);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefProton,gMeanProt,gRmsProt);
    okRes=kTRUE;
  }
  if(hOneOverPtResidVsPtTPCselITSrefPion && hOneOverPtResidVsPtTPCselITSrefPion->Integral()>0){
    FillMeanAndRms(hOneOverPtResidVsPtTPCselITSrefPion,gDum,gRelPi);
    FillMeanAndRms(hOneOverPtResidVsPtTPCselITSrefKaon,gDum,gRelK);
    FillMeanAndRms(hOneOverPtResidVsPtTPCselITSrefProton,gDum,gRelProt);
    okOneOverRes=kTRUE;
  }
  if(okRes){
    TCanvas* c2=new TCanvas("c2","Pt resol",1500,900);
    c2->Divide(3,2);
    c2->cd(1);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefPion->Draw("colz");
    // TLatex* tpi=new TLatex(0.45,0.93,"Pions");
    // tpi->SetNDC();
    // tpi->Draw();
    c2->cd(2);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefKaon->Draw("colz");
    // TLatex* tk=new TLatex(0.45,0.93,"Kaons");
    // tk->SetNDC();
    // tk->Draw();
    c2->cd(3);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefProton->Draw("colz");
    // TLatex* tpr=new TLatex(0.43,0.93,"Protons");
    // tpr->SetNDC();
    // tpr->Draw();
    c2->cd(4);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetTickx();
    gPad->SetTicky();
    gMeanPi->SetMarkerStyle(20);
    gMeanPi->SetMarkerColor(1);
    gMeanPi->SetLineColor(1);
    gMeanK->SetMarkerStyle(25);
    gMeanK->SetMarkerColor(2);
    gMeanK->SetLineColor(2);
    gMeanProt->SetMarkerStyle(27);
    gMeanProt->SetMarkerColor(4);
    gMeanProt->SetLineColor(4);
    gMeanPi->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gMeanK->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gMeanProt->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gMeanPi->GetYaxis()->SetTitle("<p_{T,reco}-p_{T,gen}> (GeV/c)");
    gMeanK->GetYaxis()->SetTitle("<p_{T,reco}-p_{T,gen}> (GeV/c)");
    gMeanProt->GetYaxis()->SetTitle("<p_{T,reco}-p_{T,gen}> (GeV/c)");
    gMeanProt->GetYaxis()->SetTitleOffset(1.6);
    gMeanProt->GetXaxis()->SetTitleOffset(1.2);
    gMeanProt->Draw("AP");
    gMeanPi->Draw("psame");
    gMeanK->Draw("psame");
    c2->cd(5);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetTickx();
    gPad->SetTicky();
    gRmsPi->SetMarkerStyle(20);
    gRmsPi->SetMarkerColor(1);
    gRmsPi->SetLineColor(1);
    gRmsK->SetMarkerStyle(25);
    gRmsK->SetMarkerColor(2);
    gRmsK->SetLineColor(2);
    gRmsProt->SetMarkerStyle(27);
    gRmsProt->SetMarkerColor(4);
    gRmsProt->SetLineColor(4);
    gRmsPi->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRmsK->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRmsProt->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRmsPi->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen}) (GeV/c)");
    gRmsK->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen}) (GeV/c)");
    gRmsProt->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen}) (GeV/c)");
    gRmsProt->GetYaxis()->SetTitleOffset(1.6);
    gRmsProt->GetXaxis()->SetTitleOffset(1.2);
    gRmsProt->Draw("AP");
    gRmsProt->Draw("AP");
    gRmsPi->Draw("psame");
    gRmsK->Draw("psame");
    TLegend* legsp=new TLegend(0.2,0.7,0.4,0.89);
    legsp->AddEntry(gRmsPi,"#pi","P");
    legsp->AddEntry(gRmsK,"K","P");
    legsp->AddEntry(gRmsProt,"p","P");
    legsp->Draw();
    c2->cd(6);
    if(okOneOverRes){
      gPad->SetLeftMargin(0.13);
      gPad->SetRightMargin(0.07);
      gPad->SetTickx();
      gPad->SetTicky();
      gRelPi->SetMarkerStyle(20);
      gRelPi->SetMarkerColor(1);
      gRelPi->SetLineColor(1);
      gRelK->SetMarkerStyle(25);
      gRelK->SetMarkerColor(2);
      gRelK->SetLineColor(2);
      gRelProt->SetMarkerStyle(27);
      gRelProt->SetMarkerColor(4);
      gRelProt->SetLineColor(4);
      gRelPi->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
      gRelK->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
      gRelProt->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
      gRelPi->GetYaxis()->SetTitle("#sigma(p_{T} (1/p_{T,reco}-1/p_{T,gen}))");
      gRelK->GetYaxis()->SetTitle("#sigma(p_{T} (1/p_{T,reco}-1/p_{T,gen}))");
      gRelPi->GetYaxis()->SetTitle("#sigma(p_{T} (1/p_{T,reco}-1/p_{T,gen}))");
      gRelPi->GetYaxis()->SetTitleOffset(1.6);
      gRelPi->GetXaxis()->SetTitleOffset(1.2);
      gRelPi->SetMinimum(0.);
      gRelPi->SetMaximum(0.05);
      gRelPi->Draw("AP");
      gRelProt->Draw("psame");
      gRelK->Draw("psame");
    }
    c2->SaveAs("PtResolRecoMinusGen.png");
  }
  TH3F*	hInvMassK0s3d=(TH3F*)l->FindObject("hInvMassK0s");
  TH3F*	hInvMassLambda3d=(TH3F*)l->FindObject("hInvMassLambda");
  TH3F*	hInvMassAntiLambda3d=(TH3F*)l->FindObject("hInvMassAntiLambda");

  // integrated histos
  TH1D* hInvMassK0s=hInvMassK0s3d->ProjectionX("hInvMassK0s1d");
  TH1D* hInvMassLambda=hInvMassLambda3d->ProjectionX("hInvMassLambda1d");
  TH1D* hInvMassAntiLambda=hInvMassAntiLambda3d->ProjectionX("hInvMassAntiLambda1d");
  hInvMassLambda->SetTitle("#Lambda - all p_{T}");
  hInvMassAntiLambda->SetTitle("#bar{#Lambda} - all p_{T}");
  hInvMassK0s->SetTitle("K^{0}_{S} - all p_{T}");
  hInvMassLambda->SetStats(0);
  hInvMassAntiLambda->SetStats(0);
  hInvMassK0s->SetStats(0);

  // lambda histos vs. radius
  Int_t z1=hInvMassLambda3d->GetZaxis()->FindBin(2.99);
  TH1D* hInvMassLambdaR1=hInvMassLambda3d->ProjectionX("hInvMassLambda1dR1",0,-1,1,z1);
  TH1D* hInvMassAntiLambdaR1=hInvMassAntiLambda3d->ProjectionX("hInvMassAntiLambda1dR1",0,-1,1,z1);
  Int_t z2=hInvMassLambda3d->GetZaxis()->FindBin(5.99);
  TH1D* hInvMassLambdaR2=hInvMassLambda3d->ProjectionX("hInvMassLambda1dR2",0,-1,z1+1,z2);
  TH1D* hInvMassAntiLambdaR2=hInvMassAntiLambda3d->ProjectionX("hInvMassAntiLambda1dR2",0,-1,z1,z2);
  Int_t z3=hInvMassLambda3d->GetZaxis()->FindBin(8.01);
  Int_t z4=hInvMassLambda3d->GetZaxis()->FindBin(22.99);
  TH1D* hInvMassLambdaR3=hInvMassLambda3d->ProjectionX("hInvMassLambda1dR3",0,-1,z3,z4);
  TH1D* hInvMassAntiLambdaR3=hInvMassAntiLambda3d->ProjectionX("hInvMassAntiLambda1dR3",0,-1,z3,z4);
  Int_t z5=hInvMassLambda3d->GetZaxis()->FindBin(28.01);
  Int_t z6=hInvMassLambda3d->GetZaxis()->FindBin(42.99);
  TH1D* hInvMassLambdaR4=hInvMassLambda3d->ProjectionX("hInvMassLambda1dR4",0,-1,z5,z6);
  TH1D* hInvMassAntiLambdaR4=hInvMassAntiLambda3d->ProjectionX("hInvMassAntiLambda1dR4",0,-1,z5,z6);


  // K0s histos vs. pt
  Int_t p1=hInvMassK0s3d->GetYaxis()->FindBin(0.499);
  TH1D* hInvMassK0sP1=hInvMassK0s3d->ProjectionX("hInvMassK0sP1",1,p1,0,-1);
  Int_t p2=hInvMassK0s3d->GetYaxis()->FindBin(0.999);
  TH1D* hInvMassK0sP2=hInvMassK0s3d->ProjectionX("hInvMassK0sP2",p1+1,p2,0,-1);
  Int_t p3=hInvMassK0s3d->GetYaxis()->FindBin(2.999);
  TH1D* hInvMassK0sP3=hInvMassK0s3d->ProjectionX("hInvMassK0sP3",p2+1,p3,0,-1);
  Int_t p4=hInvMassK0s3d->GetYaxis()->FindBin(5.01);
  Int_t p5=hInvMassK0s3d->GetYaxis()->FindBin(9.99);
  TH1D* hInvMassK0sP4=hInvMassK0s3d->ProjectionX("hInvMassK0sP4",p4,p5,0,-1);
  


  TF1* fmassk0=new TF1("fmassk0","[0]+[1]*x+[2]/sqrt(2.*TMath::Pi())/[4]*TMath::Exp(-0.5*(x-[3])*(x-[3])/[4]/[4])",0.46,0.52);
  fmassk0->SetLineWidth(2);
  fmassk0->SetLineColor(kMagenta+1);

  TF1* fmassL=new TF1("fmassL","[0]+[1]*x+[2]/sqrt(2.*TMath::Pi())/[4]*TMath::Exp(-0.5*(x-[3])*(x-[3])/[4]/[4])",1.10,1.13);
  fmassL->SetLineWidth(2);
  fmassL->SetLineColor(kRed+1);

  TCanvas* cv0=new TCanvas("cv0","V0s - all pt",1500,500);
  cv0->Divide(3,1);
  cv0->cd(1);
  hInvMassK0s->Draw();
  InitFuncAndFit(hInvMassK0s,fmassk0,kTRUE);
  cv0->cd(2);
  hInvMassLambda->Draw();
  InitFuncAndFit(hInvMassLambda,fmassL,kFALSE);
  cv0->cd(3);
  hInvMassAntiLambda->Draw();
  InitFuncAndFit(hInvMassAntiLambda,fmassL,kFALSE);
  cv0->SaveAs("MassSpectraV0-integrated.png");

  TCanvas* clam=new TCanvas("clam","Lambda vs R",1400,900);
  clam->Divide(2,2);
  clam->cd(1);
  hInvMassLambdaR1->Draw();
  InitFuncAndFit(hInvMassLambdaR1,fmassL,kFALSE);
  TLatex* tr1=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(1)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z1+1))));
  tr1->SetNDC();
  tr1->SetTextFont(43);
  tr1->SetTextSize(26);
  tr1->Draw();
  clam->cd(2);
  hInvMassLambdaR2->Draw();
  InitFuncAndFit(hInvMassLambdaR2,fmassL,kFALSE);
  TLatex* tr2=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z1+1)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z2+1))));
  tr2->SetNDC();
  tr2->SetTextFont(43);
  tr2->SetTextSize(26);
  tr2->Draw();
  clam->cd(3);
  hInvMassLambdaR3->Draw();
  InitFuncAndFit(hInvMassLambdaR3,fmassL,kFALSE);
  TLatex* tr3=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z3)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z4+1))));
  tr3->SetNDC();
  tr3->SetTextFont(43);
  tr3->SetTextSize(26);
  tr3->Draw();
  clam->cd(4);
  hInvMassLambdaR4->Draw();
  InitFuncAndFit(hInvMassLambdaR4,fmassL,kFALSE);
  TLatex* tr4=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z5)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z6+1))));
  tr4->SetNDC();
  tr4->SetTextFont(43);
  tr4->SetTextSize(26);
  tr4->Draw();
  clam->SaveAs("Lambda-MassSpectra-VsR.png");

  TCanvas* ck0=new TCanvas("ck0","K0s vs. pt",1400,900);
  ck0->Divide(2,2);
  ck0->cd(1);
  hInvMassK0sP1->Draw();
  InitFuncAndFit(hInvMassK0sP1,fmassk0,kTRUE);
  TLatex* tp1=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(1),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p1+1)));
  tp1->SetNDC();
  tp1->SetTextFont(43);
  tp1->SetTextSize(26);
  tp1->Draw();
  ck0->cd(2);
  hInvMassK0sP2->Draw();
  InitFuncAndFit(hInvMassK0sP2,fmassk0,kTRUE);
  TLatex* tp2=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p1+1),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p2+1)));
  tp2->SetNDC();
  tp2->SetTextFont(43);
  tp2->SetTextSize(26);
  tp2->Draw();
  ck0->cd(3);
  hInvMassK0sP3->Draw();
  InitFuncAndFit(hInvMassK0sP3,fmassk0,kTRUE);
  TLatex* tp3=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p2+1),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p3+1)));
  tp3->SetNDC();
  tp3->SetTextFont(43);
  tp3->SetTextSize(26);
  tp3->Draw();
  ck0->cd(4);
  hInvMassK0sP4->Draw();
  InitFuncAndFit(hInvMassK0sP4,fmassk0,kTRUE);
  TLatex* tp4=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p4),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p5+1)));
  tp4->SetNDC();
  tp4->SetTextFont(43);
  tp4->SetTextSize(26);
  tp4->Draw();
  ck0->SaveAs("K0s-MassSpectra-VsPt.png");
}

void FillMeanAndRms(TH2F* h2d, TGraphErrors* gMean, TGraphErrors* gRms){
  Int_t jpt=0;
  Int_t jptr=0;
  for(Int_t j=1; j<=h2d->GetNbinsX(); j++){
    Double_t pt=h2d->GetXaxis()->GetBinCenter(j);
    Double_t ept=0;
    Int_t ji=j;
    Int_t jf=j;
    if(pt>10){
      ji=j;
      jf=j+5;
      j=jf;
      pt=0.5*(h2d->GetXaxis()->GetBinLowEdge(ji)+h2d->GetXaxis()->GetBinUpEdge(jf));
      ept=pt-h2d->GetXaxis()->GetBinLowEdge(ji);
    }
    TH1D* htmp=h2d->ProjectionY("htmp",ji,jf);
    if(htmp->Integral()>40){
      htmp->Fit("gaus","Q0");
      TF1* fg=(TF1*)htmp->GetListOfFunctions()->FindObject("gaus");      
      Double_t m=htmp->GetMean();
      Double_t em=htmp->GetMeanError();
      gMean->SetPoint(jpt,pt,m);
      gMean->SetPointError(jpt,ept,em);
      ++jpt;
      Double_t r=fg->GetParameter(2);//htmp->GetRMS();
      Double_t er=fg->GetParError(2);//=htmp->GetRMSError();
      if(er/r<0.35){
	gRms->SetPoint(jptr,pt,r);
	gRms->SetPointError(jptr,ept,er);
	++jptr;
      }
      delete htmp;
    }
  }
}
 
void InitFuncAndFit(TH1D* hm, TF1* fmass, Bool_t isK0s){
  fmass->SetParameter(0,hm->GetBinContent(hm->FindBin(1.10)));
  if(isK0s)  fmass->SetParameter(0,hm->GetBinContent(hm->FindBin(0.45)));
  fmass->SetParameter(1,0.);
  //  fmass->SetParLimits(1,-99999999999,0.);
  fmass->SetParameter(2,100.);
  if(isK0s){
    fmass->SetParameter(3,0.5);
    fmass->SetParLimits(3,0.49,0.51);
    fmass->SetParameter(4,0.002);
    fmass->SetParLimits(4,0.0006,0.02);
  }else{
    fmass->SetParameter(3,1.116);
    fmass->SetParLimits(3,1.11,1.12);
    fmass->SetParameter(4,0.002);
    fmass->SetParLimits(4,0.0006,0.003);
  }

  hm->Fit(fmass,"R");
  TLatex* t1=new TLatex(0.14,0.8,Form("Mean = %.3f+-%.3f GeV/c^{2}",fmass->GetParameter(3),fmass->GetParError(3)));
  t1->SetTextSize(0.04);
  t1->SetNDC();
  t1->Draw();
  TLatex* t2=new TLatex(0.14,0.75,Form("Sigma = %.2f+-%.2f MeV/c^{2}",fmass->GetParameter(4)*1000.,fmass->GetParError(4)*1000.));
  t2->SetNDC();
  t2->SetTextSize(0.04);
  t2->Draw();
  TLatex* t3=new TLatex(0.14,0.7,Form("Yield = %.0f+-%.0f",fmass->GetParameter(2)/hm->GetBinWidth(1),fmass->GetParError(2)/hm->GetBinWidth(1)));
  t3->SetNDC();
  t3->SetTextSize(0.04);
  t3->Draw();

}


TH1D* ComputeMatchEff(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle){
  hnumer->Sumw2();
  hdenom->Sumw2();
  TH1D* hratio=(TH1D*)hnumer->Clone(name.Data());
  hratio->Divide(hnumer,hdenom,1.,1.,"B");
  hratio->SetLineColor(iCol);
  hratio->SetMarkerColor(iCol);
  hratio->SetMarkerStyle(iMarker);
  hratio->GetXaxis()->SetTitle(xtitle.Data());
  hratio->GetYaxis()->SetTitle("TPC-ITS matching efficiency");
  hratio->GetYaxis()->SetTitleOffset(1.25);
  hratio->SetMinimum(0);
  hratio->SetMaximum(1.05);
  hratio->SetStats(0);
  return hratio;
}

TH1D* ComputeRatio(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle){
  hnumer->Sumw2();
  hdenom->Sumw2();
  TH1D* hratio=(TH1D*)hnumer->Clone(name.Data());
  hratio->Divide(hnumer,hdenom);
  hratio->SetLineColor(iCol);
  hratio->SetMarkerColor(iCol);
  hratio->SetMarkerStyle(iMarker);
  hratio->GetXaxis()->SetTitle(xtitle.Data());
  hratio->GetYaxis()->SetTitle("Ratio");
  hratio->SetTitle(" ");
  hratio->GetYaxis()->SetTitleOffset(1.25);
  hratio->SetMinimum(0.);
  hratio->SetMaximum(2.);
  hratio->SetStats(0);
  return hratio;
}

TH1D* ComputeFraction(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle){
  hnumer->Sumw2();
  hdenom->Sumw2();
  TH1D* hratio=(TH1D*)hnumer->Clone(name.Data());
  TH1D* hsum=(TH1D*)hnumer->Clone("tmp");
  hsum->Add(hnumer,hdenom);
  hratio->Divide(hnumer,hsum,1.,1.,"B");
  hratio->SetLineColor(iCol);
  hratio->SetMarkerColor(iCol);
  hratio->SetMarkerStyle(iMarker);
  hratio->GetXaxis()->SetTitle(xtitle.Data());
  TString pname="pions";
  if(name.Contains("Prot")) pname="protons";
  hratio->GetYaxis()->SetTitle(Form("Fraction of %s tracked with wrong mass",pname.Data()));
  hratio->GetYaxis()->SetTitleOffset(1.25);
  hratio->SetMinimum(0.01);
  hratio->SetMaximum(1.05);
  hratio->SetStats(0);
  delete hsum;
  return hratio;
}

void DrawDistribTrHyp(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TString pname, Bool_t showStat){
  if(!showStat){
    h1->SetStats(0);
  }
  h1->SetLineColor(kGreen+1);
  h1->SetMarkerStyle(20);
  h1->SetMarkerColor(kGreen+1);
  h1->SetMarkerSize(0.8);
  h1->GetXaxis()->SetRangeUser(0.,10.);
  h1->Draw();
  if(showStat){
    gPad->Update();
    TPaveStats* tp1=(TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    tp1->SetY1NDC(0.73);
    tp1->SetY2NDC(0.92);
    tp1->SetTextColor(kGreen+1);
    gPad->Modified();
  }
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->SetMarkerSize(0.8);
  h2->SetMarkerStyle(21);
  if(h2->GetMaximum()>h1->GetMaximum()) h1->SetMaximum(1.2*h2->GetMaximum());
  if(h2->GetMinimum()<h1->GetMinimum() && h2->GetMinimum()>0) h1->SetMinimum(0.8*h2->GetMinimum());
  if(showStat){
    h2->Draw("sames"); 
    gPad->Update();
    TPaveStats* tp2=(TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
    tp2->SetY1NDC(0.53);
    tp2->SetY2NDC(0.72);
    tp2->SetTextColor(2);
    gPad->Modified();
  }else{
    h2->Draw("same"); 
  }
  h3->SetLineColor(kGreen+2);
  h3->SetMarkerColor(kGreen+2);
  h3->SetMarkerSize(0.8);
  h3->SetMarkerStyle(24);
  if(h3->GetMaximum()>h1->GetMaximum()) h1->SetMaximum(1.2*h3->GetMaximum());
  if(h3->GetMinimum()<h1->GetMinimum() && h3->GetMinimum()>0) h1->SetMinimum(0.8*h3->GetMinimum());
  if(showStat){
    h3->Draw("sames");
    gPad->Update();
    TPaveStats* tp3=(TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
    tp3->SetY1NDC(0.33);
    tp3->SetY2NDC(0.52);
    tp3->SetTextColor(kGreen+1);
    gPad->Modified();
  }else{
    h3->Draw("same");
  }
  h4->SetLineColor(kRed+1);
  h4->SetMarkerColor(kRed+1);
  h4->SetMarkerSize(0.8);
  h4->SetMarkerStyle(25);
  if(h4->GetMaximum()>h1->GetMaximum()) h1->SetMaximum(1.2*h4->GetMaximum());
  if(h4->GetMinimum()<h1->GetMinimum() && h4->GetMinimum()>0) h1->SetMinimum(0.8*h4->GetMinimum());
  if(showStat){
    h4->Draw("sames");
    gPad->Update();
    TPaveStats* tp4=(TPaveStats*)h4->GetListOfFunctions()->FindObject("stats");
    tp4->SetY1NDC(0.13);
    tp4->SetY2NDC(0.32);
    tp4->SetTextColor(kRed+1);
    gPad->Modified();
  }else{
    h4->Draw("same");
  }
  if(showStat){
    TLegend* leg=new TLegend(0.18,0.13,0.75,0.33);
    leg->AddEntry(h1,Form("%s, good hyp, TPC cuts",pname.Data()),"P");
    leg->AddEntry(h3,Form("%s, good hyp, TPC cuts+ITSrefit",pname.Data()),"P");
    leg->AddEntry(h2,Form("%s, bad hyp, TPC cuts",pname.Data()),"P");
    leg->AddEntry(h4,Form("%s, bad hyp, TPC cuts+ITSrefit",pname.Data()),"P");
    leg->Draw();
  }
}

void DrawDistrib(TH1D* h1, TH1D* h2, TH1D* h3, Bool_t showStat){
  if(!showStat){
    h1->SetStats(0);
  }
  h1->SetLineColor(1);
  h1->SetLineWidth(2);
  h1->Draw();
  h1->Draw("histosame");
  if(showStat){
    gPad->Update();
    TPaveStats* tp1=(TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    tp1->SetY1NDC(0.73);
    tp1->SetY2NDC(0.92);
    gPad->Modified();
  }
  h2->SetLineColor(2);
  h2->SetLineWidth(2);
  if(showStat){
    h2->Draw("sames"); 
    gPad->Update();
    TPaveStats* tp2=(TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
    tp2->SetY1NDC(0.53);
    tp2->SetY2NDC(0.72);
    tp2->SetTextColor(2);
    gPad->Modified();
  }else{
    h2->Draw("same"); 
  }
  h2->Draw("histosame");
  h3->SetLineColor(4);
  h3->SetLineWidth(2);
  if(showStat){
    h3->Draw("sames");
    gPad->Update();
    TPaveStats* tp3=(TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
    tp3->SetY1NDC(0.33);
    tp3->SetY2NDC(0.52);
    tp3->SetTextColor(4);
    gPad->Modified();
  }else{
    h3->Draw("same");
  }
    h3->Draw("histosame");
}
