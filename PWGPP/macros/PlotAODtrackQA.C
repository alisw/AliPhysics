#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH3F.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TMath.h>
#include "TTree.h"
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#endif

void InitFuncAndFit(TH1D* hm, TF1* fmass, Bool_t isK0s, Bool_t isMC=kFALSE);
void DrawDistrib(TH1D* h1, TH1D* h2, TH1D* h3, Bool_t showStat);
TH1D* ComputeMatchEff(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle);
TH1D* ComputeRatio(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle);
void FillMeanAndRms(TH2F* h2d, TGraphErrors* gMean, TGraphErrors* gRms);
Double_t fp2bkgk0(Double_t *x, Double_t *par);

const Int_t totTrending=36;
Float_t vecForTrend[totTrending];

void PlotAODtrackQA(TString filename="AnalysisResults.root", TString suffix="QA", Int_t runNumber=-1, TString outputForm="png"){

  TString pdfFileNames="";
  TString plotFileName="";
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);

  TString varForTrending[totTrending];
  for(Int_t jbit=0; jbit<12; jbit++){
    varForTrending[2*jbit]=Form("averPtFB%d",jbit);
    varForTrending[2*jbit+1]=Form("RatioFB%dFB0",jbit);
  }
  varForTrending[24]="MatchEffPt350EtaPos";
  varForTrending[25]="MatchEffPt1000EtaPos";
  varForTrending[26]="MatchEffPt4000EtaPos";
  varForTrending[27]="MatchEffPt350EtaNeg";
  varForTrending[28]="MatchEffPt1000EtaNeg";
  varForTrending[29]="MatchEffPt4000EtaNeg";
  varForTrending[30]="MatchEffSPDPt350EtaPos";
  varForTrending[31]="MatchEffSPDPt1000EtaPos";
  varForTrending[32]="MatchEffSPDPt4000EtaPos";
  varForTrending[33]="MatchEffSPDPt350EtaNeg";
  varForTrending[34]="MatchEffSPDPt1000EtaNeg";
  varForTrending[35]="MatchEffSPDPt4000EtaNeg";

  TTree* trtree=new TTree("trendingTrack","tree of trending variables");
  trtree->Branch("nrun",&runNumber,"nrun/I");
  for(Int_t j=0; j<totTrending; j++){
    trtree->Branch(varForTrending[j].Data(),&vecForTrend[j],Form("%s/F",varForTrending[j].Data()));
    vecForTrend[j]=-99.;
  }
  Bool_t isMC=kFALSE; // set automatically based on histos filled

  TFile* f=new TFile(filename.Data());
  TDirectoryFile* df=(TDirectoryFile*)f->Get("CheckAODTracks");
  if(!df){
    printf("Directory CheckAODTracks not found in file %s\n",filename.Data());
    return;
  }
  TList* l=(TList*)df->Get(Form("clistCheckAODTracks%s",suffix.Data()));
  if(!l){
    printf("TList clistCheckAODTracks%s not found in file %s\n",suffix.Data(),filename.Data());
    return;    
  }

  
  TH1F* hNEvents=(TH1F*)l->FindObject("hNEvents");
  Int_t nSelectedEvents=hNEvents->GetBinContent(6);

  TH3F* hEtaPhiPtTPCsel=(TH3F*)l->FindObject("hEtaPhiPtTPCsel");
  TH3F* hEtaPhiPtTPCselITSref=(TH3F*)l->FindObject("hEtaPhiPtTPCselITSref");
  TH3F* hEtaPhiPtTPCselSPDany=(TH3F*)l->FindObject("hEtaPhiPtTPCselSPDany");
  TH3F* hEtaPhiPtTPCselTOFbc=(TH3F*)l->FindObject("hEtaPhiPtTPCselTOFbc");
  TH3F* hEtaPhiPtTPCselITSrefTOFbc=(TH3F*)l->FindObject("hEtaPhiPtTPCselITSrefTOFbc");
  TH3F* hEtaPhiPtTPCselSPDanyTOFbc=(TH3F*)l->FindObject("hEtaPhiPtTPCselSPDanyTOFbc");
  TH3F* hEtaPhiPtPosChargeTPCsel=(TH3F*)l->FindObject("hEtaPhiPtPosChargeTPCsel");
  TH3F* hEtaPhiPtPosChargeTPCselITSref=(TH3F*)l->FindObject("hEtaPhiPtPosChargeTPCselITSref");
  TH3F* hEtaPhiPtPosChargeTPCselSPDany=(TH3F*)l->FindObject("hEtaPhiPtPosChargeTPCselSPDany");
  TH3F* hEtaPhiPtNegChargeTPCsel=(TH3F*)l->FindObject("hEtaPhiPtNegChargeTPCsel");
  TH3F* hEtaPhiPtNegChargeTPCselITSref=(TH3F*)l->FindObject("hEtaPhiPtNegChargeTPCselITSref");
  TH3F* hEtaPhiPtNegChargeTPCselSPDany=(TH3F*)l->FindObject("hEtaPhiPtNegChargeTPCselSPDany");

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
  TH1D* hPhiEtaNegTPCselLowPtTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionY("hPhiEtaNegTPCselLowPtTOFbc",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselLowPtTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionY("hPhiEtaPosTPCselLowPtTOFbc",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselHighPtTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionY("hPhiEtaNegTPCselHighPtTOFbc",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselHighPtTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionY("hPhiEtaPosTPCselHighPtTOFbc",eta0p,etamax,ptone,ptten);


  TH1D* hPhiEtaNegTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaNegTPCselITSref",etamin,eta0m);
  TH1D* hPhiEtaPosTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaPosTPCselITSref",eta0p,etamax);
  TH1D* hPtEtaNegTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionZ("hPtEtaNegTPCselITSref",etamin,eta0m);
  TH1D* hPtEtaPosTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionZ("hPtEtaPosTPCselITSref",eta0p,etamax);
  TH1D* hPhiEtaNegTPCselITSrefLowPt=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaNegTPCselITSrefLowPt",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselITSrefLowPt=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaPosTPCselITSrefLowPt",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselITSrefHighPt=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaNegTPCselITSrefHighPt",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselITSrefHighPt=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaPosTPCselITSrefHighPt",eta0p,etamax,ptone,ptten);
  TH1D* hPhiEtaNegTPCselITSrefLowPtTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionY("hPhiEtaNegTPCselITSrefLowPtTOFbc",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselITSrefLowPtTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionY("hPhiEtaPosTPCselITSrefLowPtTOFbc",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselITSrefHighPtTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionY("hPhiEtaNegTPCselITSrefHighPtTOFbc",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselITSrefHighPtTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionY("hPhiEtaPosTPCselITSrefHighPtTOFbc",eta0p,etamax,ptone,ptten);

  TH1D* hPhiEtaNegTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaNegTPCselSPDany",etamin,eta0m);
  TH1D* hPhiEtaPosTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaPosTPCselSPDany",eta0p,etamax);
  TH1D* hPtEtaNegTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionZ("hPtEtaNegTPCselSPDany",etamin,eta0m);
  TH1D* hPtEtaPosTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionZ("hPtEtaPosTPCselSPDany",eta0p,etamax);
  TH1D* hPhiEtaNegTPCselSPDanyLowPt=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaNegTPCselSPDanyLowPt",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselSPDanyLowPt=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaPosTPCselSPDanyLowPt",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselSPDanyHighPt=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaNegTPCselSPDanyHighPt",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselSPDanyHighPt=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaPosTPCselSPDanyHighPt",eta0p,etamax,ptone,ptten);
  TH1D* hPhiEtaNegTPCselSPDanyLowPtTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionY("hPhiEtaNegTPCselSPDanyLowPtTOFbc",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselSPDanyLowPtTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionY("hPhiEtaPosTPCselSPDanyLowPtTOFbc",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselSPDanyHighPtTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionY("hPhiEtaNegTPCselSPDanyHighPtTOFbc",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselSPDanyHighPtTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionY("hPhiEtaPosTPCselSPDanyHighPtTOFbc",eta0p,etamax,ptone,ptten);

  TH1D* hPtEtaNegTPCselTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionZ("hPtEtaNegTPCselTOFbc",etamin,eta0m);
  TH1D* hPtEtaPosTPCselTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionZ("hPtEtaPosTPCselTOFbc",eta0p,etamax);
  TH1D* hPtEtaNegTPCselITSrefTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionZ("hPtEtaNegTPCselITSrefTOFbc",etamin,eta0m);
  TH1D* hPtEtaPosTPCselITSrefTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionZ("hPtEtaPosTPCselITSrefTOFbc",eta0p,etamax);
  TH1D* hPtEtaNegTPCselSPDanyTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionZ("hPtEtaNegTPCselSPDanyTOFbc",etamin,eta0m);
  TH1D* hPtEtaPosTPCselSPDanyTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionZ("hPtEtaPosTPCselSPDanyTOFbc",eta0p,etamax);

  hPhiEtaNegTPCsel->SetMinimum(0);
  hPhiEtaPosTPCsel->SetMinimum(0);
  hPhiEtaNegTPCsel->SetTitle("#varphi tracks - #eta<0 - all p_{T}");
  hPhiEtaPosTPCsel->SetTitle("#varphi tracks - #eta>0 - all p_{T}");
  hPhiEtaNegTPCsel->GetXaxis()->SetTitle("#varphi (rad)");
  hPhiEtaPosTPCsel->GetXaxis()->SetTitle("#varphi (rad)");
  hPtEtaNegTPCsel->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtEtaPosTPCsel->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtEtaNegTPCsel->SetTitle("p_{T} tracks - TPC cuts - #eta<0");
  hPtEtaPosTPCsel->SetTitle("p_{T} tracks - TPC cuts - #eta>0");
  hPtEtaNegTPCselSPDany->SetTitle("p_{T} tracks - TPC cuts, SPD any - #eta<0");
  hPtEtaPosTPCselSPDany->SetTitle("p_{T} tracks - TPC cuts, SPD any - #eta>0");
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

    TCanvas* cdist=new TCanvas("cdist","Pt Distrib TPC sel",900,900);
    cdist->Divide(2,2);
    cdist->cd(1);
    gPad->SetLogy();
    DrawDistrib(hPtEtaNegTPCsel,hPtEtaNegPosChargeTPCsel,hPtEtaNegNegChargeTPCsel,kTRUE);
    TLegend* legd=new TLegend(0.46,0.66,0.79,0.89);
    legd->AddEntry(hPtEtaNegTPCsel,"All charges","L")->SetTextColor(hPtEtaNegTPCsel->GetLineColor());
    legd->AddEntry(hPtEtaNegPosChargeTPCsel,"Positive","L")->SetTextColor(hPtEtaNegPosChargeTPCsel->GetLineColor());
    legd->AddEntry(hPtEtaNegNegChargeTPCsel,"Negative","L")->SetTextColor(hPtEtaNegNegChargeTPCsel->GetLineColor());
    legd->Draw();
    cdist->cd(2);
    gPad->SetLogy();
    DrawDistrib(hPtEtaPosTPCsel,hPtEtaPosPosChargeTPCsel,hPtEtaPosNegChargeTPCsel,kTRUE);
    cdist->cd(3);
    hRatioPosNegEtaNegTPCsel->Draw();
    cdist->cd(4);
    hRatioPosNegEtaPosTPCsel->Draw();
    plotFileName=Form("TracksPtDistrib-TPCsel.%s",outputForm.Data());
    cdist->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

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
    hRatioPosNegEtaNegTPCselSPDany->Draw();
    cdists->cd(4);
    hRatioPosNegEtaPosTPCselSPDany->Draw();
    plotFileName=Form("TracksPtDistrib-TPCselSPDany.%s",outputForm.Data());
    cdists->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());
    TFile* outtrsp=new TFile("TrackPtSpectra.root","recreate");
    hPtEtaNegTPCsel->Write();
    hPtEtaNegPosChargeTPCsel->Write();
    hPtEtaNegNegChargeTPCsel->Write();
    hPtEtaPosTPCsel->Write();
    hPtEtaPosPosChargeTPCsel->Write();
    hPtEtaPosNegChargeTPCsel->Write();
    hPtEtaNegTPCselSPDany->Write();
    hPtEtaNegPosChargeTPCselSPDany->Write();
    hPtEtaNegNegChargeTPCselSPDany->Write();
    hPtEtaPosTPCselSPDany->Write();
    hPtEtaPosPosChargeTPCselSPDany->Write();
    hPtEtaPosNegChargeTPCselSPDany->Write();
    hNEvents->Write();
    outtrsp->Close();
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
    plotFileName=Form("TracksPtPhiDistrib.%s",outputForm.Data());
    cdist->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());
  }

  TCanvas* cdist2=new TCanvas("cdist2","Phi Distrib",900,900);
  cdist2->Divide(2,2);
  cdist2->cd(1);
  DrawDistrib(hPhiEtaNegTPCselLowPt,hPhiEtaNegTPCselITSrefLowPt,hPhiEtaNegTPCselSPDanyLowPt,kFALSE);
  TLegend* legf=new TLegend(0.16,0.16,0.5,0.4);
  legf->AddEntry(hPhiEtaNegTPCselLowPt,"TPC only cuts","L")->SetTextColor(hPhiEtaNegTPCselLowPt->GetLineColor());
  legf->AddEntry(hPhiEtaNegTPCselITSrefLowPt,"ITSrefit","L")->SetTextColor(hPhiEtaNegTPCselITSrefLowPt->GetLineColor());
  legf->AddEntry(hPhiEtaNegTPCselSPDanyLowPt,"SPD any","L")->SetTextColor(hPhiEtaNegTPCselSPDanyLowPt->GetLineColor());
  legf->Draw();  
  cdist2->cd(2);
  DrawDistrib(hPhiEtaPosTPCselLowPt,hPhiEtaPosTPCselITSrefLowPt,hPhiEtaPosTPCselSPDanyLowPt,kFALSE);
  cdist2->cd(3);
  DrawDistrib(hPhiEtaNegTPCselHighPt,hPhiEtaNegTPCselITSrefHighPt,hPhiEtaNegTPCselSPDanyHighPt,kFALSE);
  cdist2->cd(4);
  DrawDistrib(hPhiEtaPosTPCselHighPt,hPhiEtaPosTPCselITSrefHighPt,hPhiEtaPosTPCselSPDanyHighPt,kFALSE);
  plotFileName=Form("TracksPhiDistrib-PtBins.%s",outputForm.Data());
  cdist2->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  TH1D* hMatchEffVsPtNegEta=ComputeMatchEff(hPtEtaNegTPCselITSref,hPtEtaNegTPCsel,"hMatchEffVsPtNegEta",1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEta=ComputeMatchEff(hPtEtaPosTPCselITSref,hPtEtaPosTPCsel,"hMatchEffVsPtPosEta",1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtNegEtaSPDany=ComputeMatchEff(hPtEtaNegTPCselSPDany,hPtEtaNegTPCsel,"hMatchEffVsPtNegEtaSPDAny",kBlue-7,33,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaSPDany=ComputeMatchEff(hPtEtaPosTPCselSPDany,hPtEtaPosTPCsel,"hMatchEffVsPtPosEtaSPDAny",kBlue-7,33,"p_{T} (GeV/c)");


  TH1D* hMatchEffVsPtNegEtaTOFbc=ComputeMatchEff(hPtEtaNegTPCselITSrefTOFbc,hPtEtaNegTPCselTOFbc,"hMatchEffVsPtNegEtaTOFbc",kRed+1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaTOFbc=ComputeMatchEff(hPtEtaPosTPCselITSrefTOFbc,hPtEtaPosTPCselTOFbc,"hMatchEffVsPtPosEtaTOFbc",kRed+1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtNegEtaSPDanyTOFbc=ComputeMatchEff(hPtEtaNegTPCselSPDanyTOFbc,hPtEtaNegTPCselTOFbc,"hMatchEffVsPtNegEtaSPDAnyTOFbc",kGreen+2,33,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaSPDanyTOFbc=ComputeMatchEff(hPtEtaPosTPCselSPDanyTOFbc,hPtEtaPosTPCselTOFbc,"hMatchEffVsPtPosEtaSPDAnyTOFbc",kGreen+2,33,"p_{T} (GeV/c)");


  TH1D* hMatchEffVsPhiNegEta=ComputeMatchEff(hPhiEtaNegTPCselITSref,hPhiEtaNegTPCsel,"hMatchEffVsPhiNegEta",1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEta=ComputeMatchEff(hPhiEtaPosTPCselITSref,hPhiEtaPosTPCsel,"hMatchEffVsPhiPosEta",1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiNegEtaSPDany=ComputeMatchEff(hPhiEtaNegTPCselSPDany,hPhiEtaNegTPCsel,"hMatchEffVsPhiNegEtaSPDAny",kBlue-7,33,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaSPDany=ComputeMatchEff(hPhiEtaPosTPCselSPDany,hPhiEtaPosTPCsel,"hMatchEffVsPhiPosEtaSPDAny",kBlue-7,33,"#varphi (rad)");

  TH1D* hMatchEffVsPhiNegEtaLowPt=ComputeMatchEff(hPhiEtaNegTPCselITSrefLowPt,hPhiEtaNegTPCselLowPt,"hMatchEffVsPhiNegEtaLowPt",1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaLowPt=ComputeMatchEff(hPhiEtaPosTPCselITSrefLowPt,hPhiEtaPosTPCselLowPt,"hMatchEffVsPhiPosEtaLowPt",1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiNegEtaSPDanyLowPt=ComputeMatchEff(hPhiEtaNegTPCselSPDanyLowPt,hPhiEtaNegTPCselLowPt,"hMatchEffVsPhiNegEtaSPDAnyLowPt",kBlue-7,33,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaSPDanyLowPt=ComputeMatchEff(hPhiEtaPosTPCselSPDanyLowPt,hPhiEtaPosTPCselLowPt,"hMatchEffVsPhiPosEtaSPDAnyLowPt",kBlue-7,33,"#varphi (rad)");
  TH1D* hMatchEffVsPhiNegEtaHighPt=ComputeMatchEff(hPhiEtaNegTPCselITSrefHighPt,hPhiEtaNegTPCselHighPt,"hMatchEffVsPhiNegEtaHighPt",1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaHighPt=ComputeMatchEff(hPhiEtaPosTPCselITSrefHighPt,hPhiEtaPosTPCselHighPt,"hMatchEffVsPhiPosEtaHighPt",1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiNegEtaSPDanyHighPt=ComputeMatchEff(hPhiEtaNegTPCselSPDanyHighPt,hPhiEtaNegTPCselHighPt,"hMatchEffVsPhiNegEtaSPDAnyHighPt",kBlue-7,33,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaSPDanyHighPt=ComputeMatchEff(hPhiEtaPosTPCselSPDanyHighPt,hPhiEtaPosTPCselHighPt,"hMatchEffVsPhiPosEtaSPDAnyHighPt",kBlue-7,33,"#varphi (rad)");

  TH1D* hMatchEffVsPhiNegEtaLowPtTOFbc=ComputeMatchEff(hPhiEtaNegTPCselITSrefLowPtTOFbc,hPhiEtaNegTPCselLowPtTOFbc,"hMatchEffVsPhiNegEtaLowPtTOFbc",kRed+1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaLowPtTOFbc=ComputeMatchEff(hPhiEtaPosTPCselITSrefLowPtTOFbc,hPhiEtaPosTPCselLowPtTOFbc,"hMatchEffVsPhiPosEtaLowPtTOFbc",kRed+1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiNegEtaSPDanyLowPtTOFbc=ComputeMatchEff(hPhiEtaNegTPCselSPDanyLowPtTOFbc,hPhiEtaNegTPCselLowPtTOFbc,"hMatchEffVsPhiNegEtaSPDAnyLowPtTOFbc",kGreen+2,33,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaSPDanyLowPtTOFbc=ComputeMatchEff(hPhiEtaPosTPCselSPDanyLowPtTOFbc,hPhiEtaPosTPCselLowPtTOFbc,"hMatchEffVsPhiPosEtaSPDAnyLowPtTOFbc",kGreen+2,33,"#varphi (rad)");
  TH1D* hMatchEffVsPhiNegEtaHighPtTOFbc=ComputeMatchEff(hPhiEtaNegTPCselITSrefHighPtTOFbc,hPhiEtaNegTPCselHighPtTOFbc,"hMatchEffVsPhiNegEtaHighPtTOFbc",kRed+1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaHighPtTOFbc=ComputeMatchEff(hPhiEtaPosTPCselITSrefHighPtTOFbc,hPhiEtaPosTPCselHighPtTOFbc,"hMatchEffVsPhiPosEtaHighPtTOFbc",kRed+1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiNegEtaSPDanyHighPtTOFbc=ComputeMatchEff(hPhiEtaNegTPCselSPDanyHighPtTOFbc,hPhiEtaNegTPCselHighPtTOFbc,"hMatchEffVsPhiNegEtaSPDAnyHighPtTOFbc",kGreen+2,33,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaSPDanyHighPtTOFbc=ComputeMatchEff(hPhiEtaPosTPCselSPDanyHighPtTOFbc,hPhiEtaPosTPCselHighPtTOFbc,"hMatchEffVsPhiPosEtaSPDAnyHighPtTOFbc",kGreen+2,33,"#varphi (rad)");


  hMatchEffVsPtNegEta->SetTitle("#eta<0");
  hMatchEffVsPtPosEta->SetTitle("#eta>0");
  hMatchEffVsPhiNegEta->SetTitle("#eta<0 - all p_{T}");
  hMatchEffVsPhiPosEta->SetTitle("#eta>0 - all p_{T}");
  hMatchEffVsPtNegEtaSPDany->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaSPDany->SetTitle("#eta>0");
  hMatchEffVsPhiNegEtaSPDany->SetTitle("#eta<0 - all p_{T}");
  hMatchEffVsPhiPosEtaSPDany->SetTitle("#eta>0 - all p_{T}");

  hMatchEffVsPhiNegEtaLowPt->SetTitle("#eta<0 - 0.4<p_{T}<0.7 GeV/c");
  hMatchEffVsPhiPosEtaLowPt->SetTitle("#eta>0 - 0.4<p_{T}<0.7 GeV/c");
  hMatchEffVsPhiNegEtaHighPt->SetTitle("#eta<0 - 1<p_{T}<10 GeV/c");
  hMatchEffVsPhiPosEtaHighPt->SetTitle("#eta>0 - 1<p_{T}<10 GeV/c");

  hMatchEffVsPtNegEtaTOFbc->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaTOFbc->SetTitle("#eta>0");
  hMatchEffVsPtNegEtaSPDanyTOFbc->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaSPDanyTOFbc->SetTitle("#eta>0");

  Int_t theBin350=hMatchEffVsPtNegEta->GetXaxis()->FindBin(0.35);
  Int_t theBin950=hMatchEffVsPtNegEta->GetXaxis()->FindBin(0.95);
  Int_t theBin3950=hMatchEffVsPtNegEta->GetXaxis()->FindBin(3.95);
  vecForTrend[24]=hMatchEffVsPtPosEta->GetBinContent(theBin350);
  vecForTrend[25]=hMatchEffVsPtPosEta->GetBinContent(theBin950);
  vecForTrend[26]=hMatchEffVsPtPosEta->GetBinContent(theBin3950);
  vecForTrend[27]=hMatchEffVsPtNegEta->GetBinContent(theBin350);
  vecForTrend[28]=hMatchEffVsPtNegEta->GetBinContent(theBin950);
  vecForTrend[29]=hMatchEffVsPtNegEta->GetBinContent(theBin3950);
  vecForTrend[30]=hMatchEffVsPtPosEtaSPDany->GetBinContent(theBin350);
  vecForTrend[31]=hMatchEffVsPtPosEtaSPDany->GetBinContent(theBin950);
  vecForTrend[32]=hMatchEffVsPtPosEtaSPDany->GetBinContent(theBin3950);
  vecForTrend[33]=hMatchEffVsPtNegEtaSPDany->GetBinContent(theBin350);
  vecForTrend[34]=hMatchEffVsPtNegEtaSPDany->GetBinContent(theBin950);
  vecForTrend[35]=hMatchEffVsPtNegEtaSPDany->GetBinContent(theBin3950);
 
  TCanvas* cme=new TCanvas("cme","MatchEff Pt",900,900);
  cme->Divide(2,2);
  cme->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtNegEta->Draw("PE");
  hMatchEffVsPtNegEtaSPDany->Draw("samepe");
  TLegend* leg=new TLegend(0.27,0.17,0.6,0.39);
  leg->AddEntry(hMatchEffVsPtNegEta,"ITSrefit","P");
  leg->AddEntry(hMatchEffVsPtNegEtaSPDany,"SPD any","P");
  leg->Draw();
  cme->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtPosEta->Draw("PE");
  hMatchEffVsPtPosEtaSPDany->Draw("samepe");
  cme->cd(3);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtNegEtaTOFbc->Draw("PE");
  hMatchEffVsPtNegEtaSPDanyTOFbc->Draw("samepe");
  TLegend* legt=new TLegend(0.27,0.17,0.89,0.39);
  legt->AddEntry(hMatchEffVsPtNegEtaTOFbc,"ITSrefit, TOF bc=0","P");
  legt->AddEntry(hMatchEffVsPtNegEtaSPDanyTOFbc,"SPD any, TOF bc=0","P");
  legt->Draw();
  cme->cd(4);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtPosEtaTOFbc->Draw("PE");
  hMatchEffVsPtPosEtaSPDanyTOFbc->Draw("samepe");
  plotFileName=Form("MatchEff.%s",outputForm.Data());
  cme->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  TCanvas* cme2=new TCanvas("cme2","MatchEff Phi",900,900);
  cme2->Divide(2,2);
  cme2->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPhiNegEtaLowPt->Draw("PE");
  hMatchEffVsPhiNegEtaSPDanyLowPt->Draw("samepe");
  hMatchEffVsPhiNegEtaLowPtTOFbc->Draw("samepe");
  hMatchEffVsPhiNegEtaSPDanyLowPtTOFbc->Draw("samepe");
  TLegend* legt2=new TLegend(0.17,0.14,0.89,0.29);
  legt2->SetNColumns(2);
  legt2->SetMargin(0.1);
  legt2->AddEntry(hMatchEffVsPhiNegEtaLowPt,"ITSrefit","P");
  legt2->AddEntry(hMatchEffVsPhiNegEtaSPDanyLowPt,"SPD any","P");
  legt2->AddEntry(hMatchEffVsPhiNegEtaLowPtTOFbc,"ITSrefit, TOF bc=0","P");
  legt2->AddEntry(hMatchEffVsPhiNegEtaSPDanyLowPtTOFbc,"SPD any, TOF bc=0","P");
  legt2->Draw();
  cme2->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPhiPosEtaLowPt->Draw("PE");
  hMatchEffVsPhiPosEtaSPDanyLowPt->Draw("samepe");
  hMatchEffVsPhiPosEtaLowPtTOFbc->Draw("samepe");
  hMatchEffVsPhiPosEtaSPDanyLowPtTOFbc->Draw("samepe");
  cme2->cd(3);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPhiNegEtaHighPt->Draw("PE");
  hMatchEffVsPhiNegEtaSPDanyHighPt->Draw("samepe");
  hMatchEffVsPhiNegEtaHighPtTOFbc->Draw("samepe");
  hMatchEffVsPhiNegEtaSPDanyHighPtTOFbc->Draw("samepe");
  cme2->cd(4);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPhiPosEtaHighPt->Draw("PE");
  hMatchEffVsPhiPosEtaSPDanyHighPt->Draw("samepe");
  hMatchEffVsPhiPosEtaHighPtTOFbc->Draw("samepe");
  hMatchEffVsPhiPosEtaSPDanyHighPtTOFbc->Draw("samepe");
  plotFileName=Form("MatchEffVsPhi.%s",outputForm.Data());
  cme2->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  TFile* outME=new TFile("MatchingEff.root","recreate");
  hMatchEffVsPtNegEta->Write();
  hMatchEffVsPtPosEta->Write();
  hMatchEffVsPtNegEtaSPDany->Write();
  hMatchEffVsPtPosEtaSPDany->Write();
  hMatchEffVsPtNegEtaTOFbc->Write();
  hMatchEffVsPtPosEtaTOFbc->Write();
  hMatchEffVsPtNegEtaSPDanyTOFbc->Write();
  hMatchEffVsPtPosEtaSPDanyTOFbc->Write();
  hMatchEffVsPhiNegEtaLowPt->Write();
  hMatchEffVsPhiNegEtaSPDanyLowPt->Write();
  hMatchEffVsPhiNegEtaLowPtTOFbc->Write();
  hMatchEffVsPhiNegEtaSPDanyLowPtTOFbc->Write();
  hMatchEffVsPhiPosEtaLowPt->Write();
  hMatchEffVsPhiPosEtaSPDanyLowPt->Write();
  hMatchEffVsPhiPosEtaLowPtTOFbc->Write();
  hMatchEffVsPhiPosEtaSPDanyLowPtTOFbc->Write();
  hMatchEffVsPhiNegEtaHighPt->Write();
  hMatchEffVsPhiNegEtaSPDanyHighPt->Write();
  hMatchEffVsPhiNegEtaHighPtTOFbc->Write();
  hMatchEffVsPhiNegEtaSPDanyHighPtTOFbc->Write();
  hMatchEffVsPhiPosEtaHighPt->Write();
  hMatchEffVsPhiPosEtaSPDanyHighPt->Write();
  hMatchEffVsPhiPosEtaHighPtTOFbc->Write();
  hMatchEffVsPhiPosEtaSPDanyHighPtTOFbc->Write();
  outME->Close();
  delete outME;

  TH3F* hEtaPhiPtTPCselITSrefGood=(TH3F*)l->FindObject("hEtaPhiPtTPCselITSrefGood");
  TH3F* hEtaPhiPtTPCselITSrefFake=(TH3F*)l->FindObject("hEtaPhiPtTPCselITSrefFake");
  TH1D* hPtGood=hEtaPhiPtTPCselITSrefGood->ProjectionZ("hPtGood",etamin,eta0m);
  TH1D* hPtFake=hEtaPhiPtTPCselITSrefFake->ProjectionZ("hPtFake",etamin,eta0m);
  TH1D* hPtAll=(TH1D*)hPtGood->Clone("hPtAll");
  TH1F* hratiofake=(TH1F*)hPtFake->Clone("hratiofake");
  if(hPtFake->GetEntries()>0){
    isMC=kTRUE;
    hPtAll->Add(hPtFake);
    hPtAll->SetLineColor(1);
    hPtGood->Sumw2();
    hPtAll->Sumw2();
    hPtFake->Sumw2();
    hratiofake->Divide(hPtFake,hPtAll,1,1,"B");
    hratiofake->SetStats(0);
    hPtAll->SetMinimum(hPtFake->GetMinimum()*0.9);
  }

  TH3F* hImpParXYPtMulTPCselSPDanyGood=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanyGood");
  TH3F* hImpParXYPtMulTPCselSPDanyFake=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanyFake");
  TH1D* hImpParGood=hImpParXYPtMulTPCselSPDanyGood->ProjectionY("hImpParGood");
  TH1D* hImpParFake=hImpParXYPtMulTPCselSPDanyFake->ProjectionY("hImpParFake");
  hImpParGood->SetLineColor(kGreen+1);
  hImpParFake->SetLineColor(2);
  TH1D* hImpParAllGF=(TH1D*)hImpParGood->Clone("hImpParAllGF");
  hImpParAllGF->Add(hImpParFake);
  hImpParAllGF->SetLineColor(1);
  hImpParAllGF->Sumw2();
  hImpParGood->Sumw2();
  hImpParFake->Sumw2();
  TH1F* hratiofakeip=(TH1F*)hImpParFake->Clone("hratiofakeip");
  hratiofakeip->Divide(hImpParFake,hImpParAllGF,1,1,"B");
  hratiofakeip->SetLineColor(1);
  hratiofakeip->SetStats(0);
  if(hImpParFake->Integral()>0 && hImpParGood->Integral()>0 ){
    isMC=kTRUE;
    TCanvas* c1=new TCanvas("c1","FakeGood",1200,900);
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetLogy();
    hPtAll->SetMinimum(0.5);
    hPtAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPtAll->Draw("histo");
    gPad->Update();
    TPaveStats* stp1=(TPaveStats*)hPtAll->GetListOfFunctions()->FindObject("stats");
    if(stp1){
      stp1->SetTextColor(1);
      stp1->SetY1NDC(0.73);
      stp1->SetY2NDC(0.92);
    }
    gPad->Modified();
    hPtGood->SetLineColor(kGreen+1);
    hPtGood->Draw("histosames");
    gPad->Update();
    TPaveStats* stp2=(TPaveStats*)hPtGood->GetListOfFunctions()->FindObject("stats");
    stp2->SetTextColor(kGreen+1);
    stp2->SetY1NDC(0.72);
    stp2->SetY2NDC(0.53);
    gPad->Modified();
    hPtFake->SetLineColor(2);
    hPtFake->Draw("histosames");
    gPad->Update();
    TPaveStats* stp3=(TPaveStats*)hPtFake->GetListOfFunctions()->FindObject("stats");
    stp3->SetTextColor(2);
    stp3->SetY1NDC(0.52);
    stp3->SetY2NDC(0.33);
    gPad->Modified();
    c1->cd(2);
    gPad->SetLogy();
    hImpParGood->Draw("histo");
    gPad->Update();
    TPaveStats* st1=(TPaveStats*)hImpParGood->GetListOfFunctions()->FindObject("stats");
    st1->SetTextColor(kGreen+1);
    st1->SetY1NDC(0.73);
    st1->SetY2NDC(0.92);
    gPad->Modified();
    hImpParFake->Draw("histosames"); 
    gPad->Update();
    TPaveStats* st2=(TPaveStats*)hImpParFake->GetListOfFunctions()->FindObject("stats");
    st2->SetTextColor(2);
    st2->SetY1NDC(0.72);
    st2->SetY2NDC(0.53);
    gPad->Modified();
    c1->cd(3);
    hratiofake->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hratiofake->GetYaxis()->SetTitle("Fraction of fakes");
    hratiofake->GetYaxis()->SetTitleOffset(1.2);
    hratiofake->Draw();
    c1->cd(4);
    hratiofakeip->GetYaxis()->SetTitle("Fraction of fakes");
    hratiofakeip->GetYaxis()->SetTitleOffset(1.2);
    hratiofakeip->Draw();
    plotFileName=Form("GoodFakeTracks.%s",outputForm.Data());
    c1->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());
  }
  
  TH3F* hImpParXYPtMulTPCselSPDanyPrim=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanyPrim");
  TH3F* hImpParXYPtMulTPCselSPDanySecDec=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanySecDec");
  TH3F* hImpParXYPtMulTPCselSPDanySecMat=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanySecMat");
  TH1D* hPtPrim=hImpParXYPtMulTPCselSPDanyPrim->ProjectionX("hPtPrim");
  TH1D* hPtSecDec=hImpParXYPtMulTPCselSPDanySecDec->ProjectionX("hPtSecDec");
  TH1D* hPtSecMat=hImpParXYPtMulTPCselSPDanySecMat->ProjectionX("hPtSecMat");
  TH1D* hImpParPrim=hImpParXYPtMulTPCselSPDanyPrim->ProjectionY("hImpParPrim");
  TH1D* hImpParSecDec=hImpParXYPtMulTPCselSPDanySecDec->ProjectionY("hImpParSecDec");
  TH1D* hImpParSecMat=hImpParXYPtMulTPCselSPDanySecMat->ProjectionY("hImpParSecMat");
    
  TH1D* hImpParAll=(TH1D*)hImpParSecDec->Clone("hImpParAll");
  hImpParAll->Add(hImpParSecMat);
  hImpParAll->Add(hImpParPrim);
  TH1D* hPtAllPS=(TH1D*)hPtSecDec->Clone("hPtAllPS");
  hPtAllPS->Add(hPtSecMat);
  hPtAllPS->Add(hPtPrim);
  hImpParAll->SetLineColor(1);
  hImpParPrim->SetLineColor(4);
  hImpParSecDec->SetLineColor(kOrange+1);
  hImpParSecMat->SetLineColor(6);
  hImpParPrim->SetLineWidth(2);
  hImpParSecDec->SetLineWidth(2);
  hImpParSecMat->SetLineWidth(2);
  hPtPrim->SetLineColor(4);
  hPtSecDec->SetLineColor(kOrange+1);
  hPtSecMat->SetLineColor(6);
  hPtPrim->SetLineWidth(2);
  hPtSecDec->SetLineWidth(2);
  hPtSecMat->SetLineWidth(2);
  hPtAllPS->Sumw2();
  hPtSecDec->Sumw2();
  hPtSecMat->Sumw2();
  hPtPrim->Sumw2();
  hImpParAll->Sumw2();
  hImpParSecDec->Sumw2();
  hImpParSecMat->Sumw2();
  hImpParPrim->Sumw2();
  TH1F* hratiosecdec=(TH1F*)hPtSecDec->Clone("hratiosecdec");
  hratiosecdec->Divide(hPtSecDec,hPtAllPS,1,1,"B");
  hratiosecdec->SetStats(0);
  TH1F* hratiosecdecip=(TH1F*)hImpParSecDec->Clone("hratiosecdecpi");
  hratiosecdecip->Divide(hImpParSecDec,hImpParAll,1,1,"B");
  hratiosecdecip->SetStats(0);
  TH1F* hratiosecmat=(TH1F*)hPtSecMat->Clone("hratiosecmat");
  hratiosecmat->Divide(hPtSecMat,hPtAllPS,1,1,"B");
  hratiosecmat->SetStats(0);
  TH1F* hratiosecmatip=(TH1F*)hImpParSecMat->Clone("hratiosecmatpi");
  hratiosecmatip->Divide(hImpParSecMat,hImpParAll,1,1,"B");
  hratiosecmatip->SetStats(0);
  hPtAllPS->Scale(1.,"width");
  hPtPrim->Scale(1.,"width");
  hPtSecDec->Scale(1.,"width");
  hPtSecMat->Scale(1.,"width");

  if(hImpParSecDec->Integral()>0 && hImpParPrim->Integral()>0 ){
    isMC=kTRUE;
    TCanvas* cps1=new TCanvas("cps1","SecPrim",1200,900);
    cps1->Divide(2,2);
    cps1->cd(1);
    gPad->SetLogy();
    hPtAllPS->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPtAllPS->Draw("histo");
    gPad->Update();
    TPaveStats* stp1=(TPaveStats*)hPtAllPS->GetListOfFunctions()->FindObject("stats");
    if(stp1){
      stp1->SetTextColor(1);
      stp1->SetY1NDC(0.73);
      stp1->SetY2NDC(0.92);
    }
    gPad->Modified();
    hPtPrim->Draw("histosames");
    gPad->Update();
    TPaveStats* stp2=(TPaveStats*)hPtPrim->GetListOfFunctions()->FindObject("stats");
    stp2->SetTextColor(hPtPrim->GetLineColor());
    stp2->SetY1NDC(0.72);
    stp2->SetY2NDC(0.53);
    gPad->Modified();
    hPtSecDec->Draw("histosames");
    gPad->Update();
    TPaveStats* stp3=(TPaveStats*)hPtSecDec->GetListOfFunctions()->FindObject("stats");
    stp3->SetTextColor(hPtSecDec->GetLineColor());
    stp3->SetY1NDC(0.52);
    stp3->SetY2NDC(0.33);
    gPad->Modified();
    hPtSecMat->Draw("histosames");
    gPad->Update();
    TPaveStats* stp4=(TPaveStats*)hPtSecMat->GetListOfFunctions()->FindObject("stats");
    stp4->SetTextColor(hPtSecMat->GetLineColor());
    stp4->SetY1NDC(0.32);
    stp4->SetY2NDC(0.13);
    gPad->Modified();
    cps1->cd(2);
    gPad->SetLogy();
    hImpParPrim->Draw("histo");
    gPad->Update();
    TPaveStats* sti1=(TPaveStats*)hImpParPrim->GetListOfFunctions()->FindObject("stats");
    sti1->SetTextColor(hImpParPrim->GetLineColor());
    sti1->SetY1NDC(0.73);
    sti1->SetY2NDC(0.92);
    gPad->Modified();
    hImpParSecDec->Draw("histosames"); 
    gPad->Update();
    TPaveStats* sti2=(TPaveStats*)hImpParSecDec->GetListOfFunctions()->FindObject("stats");
    sti2->SetTextColor(hImpParSecDec->GetLineColor());
    sti2->SetY1NDC(0.72);
    sti2->SetY2NDC(0.53);
    gPad->Modified();
    hImpParSecMat->Draw("histosames"); 
    gPad->Update();
    TPaveStats* sti3=(TPaveStats*)hImpParSecMat->GetListOfFunctions()->FindObject("stats");
    sti3->SetTextColor(hImpParSecMat->GetLineColor());
    sti3->SetY1NDC(0.52);
    sti3->SetY2NDC(0.33);
    gPad->Modified();
    cps1->cd(3);
    hratiosecdec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hratiosecdec->GetYaxis()->SetTitle("Fraction of secondaries");
    hratiosecdec->GetYaxis()->SetTitleOffset(1.2);
    hratiosecdec->Draw();
    hratiosecmat->Draw("same");
    cps1->cd(4);
    //    hratiosecip->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hratiosecdecip->GetYaxis()->SetTitle("Fraction of secondaries");
    hratiosecdecip->GetYaxis()->SetTitleOffset(1.2);
    hratiosecdecip->Draw();
    hratiosecmatip->Draw("same");
    plotFileName=Form("PrimSecTracks.%s",outputForm.Data());
    cps1->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());
  }
  
  TH2F* hPtResidVsPtTPCselITSrefPrim=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefPrim");
  TH2F* hPtResidVsPtTPCselITSrefSecDec=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefSecDec");
  TH2F* hPtResidVsPtTPCselITSrefSecMat=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefSecMat");  
  TH2F* hOneOverPtResidVsPtTPCselITSrefPrim=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefPrim");
  TH2F* hOneOverPtResidVsPtTPCselITSrefSecDec=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefSecDec");
  TH2F* hOneOverPtResidVsPtTPCselITSrefSecMat=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefSecMat");  
  if(hPtResidVsPtTPCselITSrefPrim){
    hPtResidVsPtTPCselITSrefPrim->SetStats(0);
    hPtResidVsPtTPCselITSrefSecDec->SetStats(0);
    hPtResidVsPtTPCselITSrefSecMat->SetStats(0);
    hPtResidVsPtTPCselITSrefPrim->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefSecDec->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefSecMat->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefPrim->GetXaxis()->SetTitleOffset(1.1);
    hPtResidVsPtTPCselITSrefSecDec->GetXaxis()->SetTitleOffset(1.1);
    hPtResidVsPtTPCselITSrefSecMat->GetXaxis()->SetTitleOffset(1.1);
  }
  if(hOneOverPtResidVsPtTPCselITSrefPrim){
    hOneOverPtResidVsPtTPCselITSrefPrim->SetStats(0);
    hOneOverPtResidVsPtTPCselITSrefSecDec->SetStats(0);
    hOneOverPtResidVsPtTPCselITSrefSecMat->SetStats(0);
    hOneOverPtResidVsPtTPCselITSrefPrim->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselITSrefSecDec->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselITSrefSecMat->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselITSrefPrim->GetXaxis()->SetTitleOffset(1.1);
    hOneOverPtResidVsPtTPCselITSrefSecDec->GetXaxis()->SetTitleOffset(1.1);
    hOneOverPtResidVsPtTPCselITSrefSecMat->GetXaxis()->SetTitleOffset(1.1);
  }
  TGraphErrors* gMeanPrim=new TGraphErrors(0);
  TGraphErrors* gMeanSecDec=new TGraphErrors(0);
  TGraphErrors* gMeanSecMat=new TGraphErrors(0);
  TGraphErrors* gRmsPrim=new TGraphErrors(0);
  TGraphErrors* gRmsSecDec=new TGraphErrors(0);
  TGraphErrors* gRmsSecMat=new TGraphErrors(0);
  TGraphErrors* gDum=new TGraphErrors(0);
  TGraphErrors* gRelPrim=new TGraphErrors(0);
  TGraphErrors* gRelSecDec=new TGraphErrors(0);
  TGraphErrors* gRelSecMat=new TGraphErrors(0);
  gMeanPrim->SetName("gMeanPrim");
  gMeanSecDec->SetName("gMeanSecDec");
  gMeanSecMat->SetName("gMeanSecMat");
  gMeanPrim->SetTitle("");
  gMeanSecDec->SetTitle("");
  gMeanSecMat->SetTitle("");
  gRmsPrim->SetName("gRmsPrim");
  gRmsSecDec->SetName("gRmsSecDec");
  gRmsSecMat->SetName("gRmsSecMat");
  gRmsPrim->SetTitle("");
  gRmsSecDec->SetTitle("");
  gRmsSecMat->SetTitle("");
  gRelPrim->SetName("gRelPrim");
  gRelSecDec->SetName("gRelSecDec");
  gRelSecMat->SetName("gRelSecMat");
  gRelPrim->SetTitle("");
  gRelSecDec->SetTitle("");
  gRelSecMat->SetTitle("");

  Bool_t okRes=kFALSE;
  Bool_t okOneOverRes=kFALSE;
  if(hPtResidVsPtTPCselITSrefPrim && hPtResidVsPtTPCselITSrefPrim->Integral()>0){
    isMC=kTRUE;
    FillMeanAndRms(hPtResidVsPtTPCselITSrefPrim,gMeanPrim,gRmsPrim);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefSecDec,gMeanSecDec,gRmsSecDec);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefSecMat,gMeanSecMat,gRmsSecMat);
    okRes=kTRUE;
  }
  if(hOneOverPtResidVsPtTPCselITSrefPrim && hOneOverPtResidVsPtTPCselITSrefPrim->Integral()>0){
    isMC=kTRUE;
    FillMeanAndRms(hOneOverPtResidVsPtTPCselITSrefPrim,gDum,gRelPrim);
    FillMeanAndRms(hOneOverPtResidVsPtTPCselITSrefSecDec,gDum,gRelSecDec);
    FillMeanAndRms(hOneOverPtResidVsPtTPCselITSrefSecMat,gDum,gRelSecMat);
    okOneOverRes=kTRUE;
  }
  if(okRes){
    TCanvas* cps2=new TCanvas("cps2","Pt resol - prim/sec",1500,900);
    cps2->Divide(3,2);
    cps2->cd(1);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefPrim->GetXaxis()->SetRangeUser(0.,20.);
    hPtResidVsPtTPCselITSrefPrim->Draw("colz");
    cps2->cd(2);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefSecDec->GetXaxis()->SetRangeUser(0.,20.);
    hPtResidVsPtTPCselITSrefSecDec->Draw("colz");
    cps2->cd(3);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefSecMat->GetXaxis()->SetRangeUser(0.,20.);
    hPtResidVsPtTPCselITSrefSecMat->Draw("colz");
    cps2->cd(4);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetTickx();
    gPad->SetTicky();
    gMeanPrim->SetMarkerStyle(21);
    gMeanPrim->SetMarkerColor(1);
    gMeanPrim->SetLineColor(1);
    gMeanSecDec->SetMarkerStyle(24);
    gMeanSecDec->SetMarkerColor(kRed+1);
    gMeanSecDec->SetLineColor(kRed+1);
    gMeanSecMat->SetMarkerStyle(28);
    gMeanSecMat->SetMarkerColor(kBlue+1);
    gMeanSecMat->SetLineColor(kBlue+1);
    gMeanPrim->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gMeanSecDec->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gMeanSecMat->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gMeanPrim->GetYaxis()->SetTitle("<p_{T,reco}-p_{T,gen}> (GeV/c)");
    gMeanSecDec->GetYaxis()->SetTitle("<p_{T,reco}-p_{T,gen}> (GeV/c)");
    gMeanSecMat->GetYaxis()->SetTitle("<p_{T,reco}-p_{T,gen}> (GeV/c)");
    gMeanPrim->GetYaxis()->SetTitleOffset(1.6);
    gMeanPrim->GetXaxis()->SetTitleOffset(1.2);
    gMeanPrim->SetMinimum(-0.04);
    gMeanPrim->SetMaximum(0.04);
    gMeanPrim->GetXaxis()->SetLimits(0.,20.);
    gMeanPrim->Draw("AP");
    gMeanSecMat->Draw("psame");
    gMeanSecDec->Draw("psame");
    cps2->cd(5);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetTickx();
    gPad->SetTicky();
    gRmsPrim->SetMarkerStyle(21);
    gRmsPrim->SetMarkerColor(1);
    gRmsPrim->SetLineColor(1);
    gRmsSecDec->SetMarkerStyle(24);
    gRmsSecDec->SetMarkerColor(kRed+1);
    gRmsSecDec->SetLineColor(kRed+1);
    gRmsSecMat->SetMarkerStyle(28);
    gRmsSecMat->SetMarkerColor(kBlue+1);
    gRmsSecMat->SetLineColor(kBlue+1);
    gRmsPrim->SetMinimum(0.);
    gRmsPrim->SetMaximum(0.4);
    gRmsPrim->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRmsSecDec->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRmsSecMat->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRmsPrim->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen}) (GeV/c)");
    gRmsSecDec->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen}) (GeV/c)");
    gRmsSecMat->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen}) (GeV/c)");
    gRmsPrim->GetYaxis()->SetTitleOffset(1.6);
    gRmsPrim->GetXaxis()->SetTitleOffset(1.2);
    gRmsPrim->GetXaxis()->SetLimits(0.,20.);
    gRmsPrim->Draw("AP");
    gRmsSecMat->Draw("psame");
    gRmsSecDec->Draw("psame");
    TLegend* legps=new TLegend(0.2,0.7,0.4,0.89);
    legps->AddEntry(gRmsPrim,"Primary","P");
    legps->AddEntry(gRmsSecDec,"Sec from decay","P");
    legps->AddEntry(gRmsSecMat,"Sec from mat","P");
    legps->Draw();
    cps2->cd(6);
    if(okOneOverRes){
      gPad->SetLeftMargin(0.13);
      gPad->SetRightMargin(0.07);
      gPad->SetTickx();
      gPad->SetTicky();
      gRelPrim->SetMarkerStyle(21);
      gRelPrim->SetMarkerColor(1);
      gRelPrim->SetLineColor(1);
      gRelSecDec->SetMarkerStyle(24);
      gRelSecDec->SetMarkerColor(kRed+1);
      gRelSecDec->SetLineColor(kRed+1);
      gRelSecMat->SetMarkerStyle(28);
      gRelSecMat->SetMarkerColor(kBlue+1);
      gRelSecMat->SetLineColor(kBlue+1);
      gRelPrim->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
      gRelSecDec->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
      gRelSecMat->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
      gRelPrim->GetYaxis()->SetTitle("#sigma(p_{T} (1/p_{T,reco}-1/p_{T,gen}))");
      gRelSecDec->GetYaxis()->SetTitle("#sigma(p_{T} (1/p_{T,reco}-1/p_{T,gen}))");
      gRelPrim->GetYaxis()->SetTitle("#sigma(p_{T} (1/p_{T,reco}-1/p_{T,gen}))");
      gRelPrim->GetYaxis()->SetTitleOffset(1.6);
      gRelPrim->GetXaxis()->SetTitleOffset(1.2);
      gRelPrim->GetXaxis()->SetLimits(0.,20.);
      gRelPrim->SetMinimum(0.);
      gRelPrim->SetMaximum(0.04);
      gRelPrim->Draw("AP");
      gRelSecMat->Draw("psame");
      gRelSecDec->Draw("psame");
    }
    plotFileName=Form("PtResidualsPrimSec.%s",outputForm.Data());
    cps2->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());
  }

  TH2F* hPtResidVsPtTPCselITSrefPion=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefpi");
  TH2F* hPtResidVsPtTPCselITSrefKaon=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefK");
  TH2F* hPtResidVsPtTPCselITSrefProton=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefp");  
  TH2F* hOneOverPtResidVsPtTPCselITSrefPion=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefpi");
  TH2F* hOneOverPtResidVsPtTPCselITSrefKaon=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefK");
  TH2F* hOneOverPtResidVsPtTPCselITSrefProton=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselITSrefp");  
  TH2F* hPtResidVsPtTPCselPion=(TH2F*)l->FindObject("hPtResidVsPtTPCselpi");
  TH2F* hOneOverPtResidVsPtTPCselPion=(TH2F*)l->FindObject("hOneOverPtResidVsPtTPCselpi");
  if(hPtResidVsPtTPCselITSrefPion){
    hPtResidVsPtTPCselITSrefPion->SetStats(0);
    hPtResidVsPtTPCselITSrefKaon->SetStats(0);
    hPtResidVsPtTPCselITSrefProton->SetStats(0);
    hPtResidVsPtTPCselPion->SetStats(0);
    hPtResidVsPtTPCselITSrefPion->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefKaon->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefProton->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselPion->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefPion->GetXaxis()->SetTitleOffset(1.1);
    hPtResidVsPtTPCselITSrefKaon->GetXaxis()->SetTitleOffset(1.1);
    hPtResidVsPtTPCselITSrefProton->GetXaxis()->SetTitleOffset(1.1);
    hPtResidVsPtTPCselPion->GetXaxis()->SetTitleOffset(1.1);
  }
  if(hOneOverPtResidVsPtTPCselITSrefPion){
    hOneOverPtResidVsPtTPCselITSrefPion->SetStats(0);
    hOneOverPtResidVsPtTPCselITSrefKaon->SetStats(0);
    hOneOverPtResidVsPtTPCselITSrefProton->SetStats(0);
    hOneOverPtResidVsPtTPCselPion->SetStats(0);
    hOneOverPtResidVsPtTPCselITSrefPion->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselITSrefKaon->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselITSrefProton->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselPion->GetYaxis()->SetTitleOffset(1.5);
    hOneOverPtResidVsPtTPCselITSrefPion->GetXaxis()->SetTitleOffset(1.1);
    hOneOverPtResidVsPtTPCselITSrefKaon->GetXaxis()->SetTitleOffset(1.1);
    hOneOverPtResidVsPtTPCselITSrefProton->GetXaxis()->SetTitleOffset(1.1);
    hOneOverPtResidVsPtTPCselPion->GetXaxis()->SetTitleOffset(1.1);
  }

  TGraphErrors* gMeanPi=new TGraphErrors(0);
  TGraphErrors* gMeanK=new TGraphErrors(0);
  TGraphErrors* gMeanProt=new TGraphErrors(0);
  TGraphErrors* gMeanPiTPC=new TGraphErrors(0);
  TGraphErrors* gRmsPi=new TGraphErrors(0);
  TGraphErrors* gRmsK=new TGraphErrors(0);
  TGraphErrors* gRmsProt=new TGraphErrors(0);
  TGraphErrors* gRmsPiTPC=new TGraphErrors(0);
  TGraphErrors* gRelPi=new TGraphErrors(0);
  TGraphErrors* gRelK=new TGraphErrors(0);
  TGraphErrors* gRelProt=new TGraphErrors(0);
  TGraphErrors* gRelPiTPC=new TGraphErrors(0);
  gMeanPi->SetName("gMeanPi");
  gMeanK->SetName("gMeanK");
  gMeanProt->SetName("gMeanProt");
  gMeanPiTPC->SetName("gMeanPiTPC");
  gMeanPi->SetTitle("");
  gMeanK->SetTitle("");
  gMeanProt->SetTitle("");
  gMeanPiTPC->SetTitle("");
  gRmsPi->SetName("gRmsPi");
  gRmsK->SetName("gRmsK");
  gRmsProt->SetName("gRmsProt");
  gRmsPiTPC->SetName("gRmsPiTPC");
  gRmsPi->SetTitle("");
  gRmsK->SetTitle("");
  gRmsProt->SetTitle("");
  gRmsPiTPC->SetTitle("");
  gRelPi->SetName("gRelPi");
  gRelK->SetName("gRelK");
  gRelProt->SetName("gRelProt");
  gRelPiTPC->SetName("gRelPiTPC");
  gRelPi->SetTitle("");
  gRelK->SetTitle("");
  gRelProt->SetTitle("");
  gRelPiTPC->SetTitle("");

  okRes=kFALSE;
  okOneOverRes=kFALSE;
  if(hPtResidVsPtTPCselITSrefPion && hPtResidVsPtTPCselITSrefPion->Integral()>0){
    FillMeanAndRms(hPtResidVsPtTPCselITSrefPion,gMeanPi,gRmsPi);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefKaon,gMeanK,gRmsK);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefProton,gMeanProt,gRmsProt);
    FillMeanAndRms(hPtResidVsPtTPCselPion,gMeanPiTPC,gRmsPiTPC);
    okRes=kTRUE;
  }
  if(hOneOverPtResidVsPtTPCselITSrefPion && hOneOverPtResidVsPtTPCselITSrefPion->Integral()>0){
    FillMeanAndRms(hOneOverPtResidVsPtTPCselITSrefPion,gDum,gRelPi);
    FillMeanAndRms(hOneOverPtResidVsPtTPCselITSrefKaon,gDum,gRelK);
    FillMeanAndRms(hOneOverPtResidVsPtTPCselITSrefProton,gDum,gRelProt);
    FillMeanAndRms(hOneOverPtResidVsPtTPCselPion,gDum,gRelPiTPC);
    okOneOverRes=kTRUE;
  }
  if(okRes){

    TCanvas* c2=new TCanvas("c2","Pt resol - species",1500,900);
    c2->Divide(3,2);
    c2->cd(1);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefPion->GetXaxis()->SetRangeUser(0.,20.);
    hPtResidVsPtTPCselITSrefPion->Draw("colz");
    // TLatex* tpi=new TLatex(0.45,0.93,"Pions");
    // tpi->SetNDC();
    // tpi->Draw();
    c2->cd(2);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefKaon->GetXaxis()->SetRangeUser(0.,20.);
    hPtResidVsPtTPCselITSrefKaon->Draw("colz");
    // TLatex* tk=new TLatex(0.45,0.93,"Kaons");
    // tk->SetNDC();
    // tk->Draw();
    c2->cd(3);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefProton->GetXaxis()->SetRangeUser(0.,20.);
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
    gMeanProt->SetMinimum(-0.04);
    gMeanProt->SetMaximum(0.04);
    gMeanProt->GetXaxis()->SetLimits(0.,20.);
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
    gRmsProt->SetMinimum(0.);
    gRmsProt->SetMaximum(0.4);
    gRmsProt->GetXaxis()->SetLimits(0.,20.);
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
      gRelPi->SetMaximum(0.04);
      gRelPi->GetXaxis()->SetLimits(0.,20.);
      gRelPi->Draw("AP");
      gRelProt->Draw("psame");
      gRelK->Draw("psame");
    }
    plotFileName=Form("PtResiduals.%s",outputForm.Data());
    c2->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

    TCanvas* c2p=new TCanvas("c2p","Pt resol - pion",1500,450);
    c2p->Divide(3,1);
    c2p->cd(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetTickx();
    gPad->SetTicky();
    gMeanPiTPC->SetMarkerStyle(26);
    gMeanPiTPC->SetMarkerColor(kMagenta+1);
    gMeanPiTPC->SetLineColor(kMagenta+1);
    gMeanPi->GetYaxis()->SetTitleOffset(1.6);
    gMeanPi->GetXaxis()->SetTitleOffset(1.2);
    gMeanPi->SetMinimum(-0.04);
    gMeanPi->SetMaximum(0.04);
    gMeanPi->Draw("AP");
    gMeanPiTPC->Draw("psame");
    c2p->cd(2);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetTickx();
    gPad->SetTicky();
    gRmsPiTPC->SetMarkerStyle(26);
    gRmsPiTPC->SetMarkerColor(kMagenta+1);
    gRmsPiTPC->SetLineColor(kMagenta+1);
    gRmsPi->SetMinimum(-0.04);
    gRmsPi->SetMaximum(0.4);
    gRmsPi->GetYaxis()->SetTitleOffset(1.6);
    gRmsPi->GetXaxis()->SetTitleOffset(1.2);
    gRmsPi->GetXaxis()->SetLimits(0.,20.);
    gRmsPi->Draw("AP");
    gRmsPiTPC->Draw("psame");
    TLegend* legpit=new TLegend(0.2,0.7,0.5,0.89);
    legpit->AddEntry(gRmsPi,"#pi, TPC+ITS","P");
    legpit->AddEntry(gRmsPiTPC,"#pi, TPC","P");
    legpit->Draw();
    c2p->cd(3);
    if(okOneOverRes){
      gPad->SetLeftMargin(0.13);
      gPad->SetRightMargin(0.07);
      gPad->SetTickx();
      gPad->SetTicky();
      gRelPiTPC->SetMarkerStyle(26);
      gRelPiTPC->SetMarkerColor(kMagenta+1);
      gRelPiTPC->SetLineColor(kMagenta+1);
      gRelPi->GetYaxis()->SetTitleOffset(1.6);
      gRelPi->GetXaxis()->SetTitleOffset(1.2);
      gRelPi->SetMinimum(0.);
      gRelPi->SetMaximum(0.04);
      gRelPi->GetXaxis()->SetLimits(0.,20.);
      gRelPi->Draw("AP");
      gRelPiTPC->Draw("psame");
    }
    plotFileName=Form("PtResidualsPionsTPCITS.%s",outputForm.Data());
    c2p->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

    TFile* outres=new TFile("PtResol.root","recreate");
    gMeanPrim->Write();
    gMeanSecDec->Write();
    gMeanSecMat->Write();
    gMeanPi->Write();
    gMeanK->Write();
    gMeanProt->Write();
    gMeanPiTPC->Write();
    gRmsPrim->Write();
    gRmsSecDec->Write();
    gRmsSecMat->Write();
    gRmsPi->Write();
    gRmsK->Write();
    gRmsProt->Write();
    gRmsPiTPC->Write();
    gRelPrim->Write();
    gRelSecDec->Write();
    gRelSecMat->Write();
    gRelPi->Write();
    gRelK->Write();
    gRelProt->Write();
    gRelPiTPC->Write();
  }


  TH2F* hFilterBits=(TH2F*)l->FindObject("hFilterBits");
  hFilterBits->SetStats(0);
  Double_t cnt0=hFilterBits->GetBinContent(1,1)+hFilterBits->GetBinContent(1,2);
  if(cnt0>0){
    for(Int_t jbit=0; jbit<12; jbit++){
      Double_t cntj=hFilterBits->GetBinContent(jbit+1,1)+hFilterBits->GetBinContent(jbit+1,2);
      vecForTrend[2*jbit+1]=cntj/cnt0;
    }
  }

  TDirectory* curDir=gDirectory;
  TFile* outfb=new TFile("FiltBitsHistos.root","recreate");
  curDir->cd();
  Int_t colors[12]={kRed+1,kRed-7,kOrange+1,kYellow+1,kGreen+1,kGreen,kCyan,kBlue+1,kMagenta,kMagenta+1,kGray+1,1};
  Int_t lstyl[12]={1,9,1,3,1,8,2,5,7,1,1,9};
  Int_t lwid[12]={2,2,2,3,2,2,3,3,3,2,2,2};

  TH1F* hNtracksFB0=(TH1F*)l->FindObject("hNtracksFb0");
  if(hNtracksFB0){
    TH1F* hMean=new TH1F("hMeanMultVsBF"," ; Filter bit ; <N_{Tracks}> (kINT7 events)",9,-0.5,8.5);
    TH1F* hRMS=new TH1F("hRMS"," ; Filter bit ; r.m.s.(N_{Tracks})",9,-0.5,8.5);
    TLegend* legmfb=new TLegend(0.6,0.45,0.89,0.89);
    TH1F* hNtracksFB[9];
    TLatex* tmeanmult[9];
    Double_t maxYaxis=1;
    for(Int_t kb=0; kb<9; kb++){
      hNtracksFB[kb]=(TH1F*)l->FindObject(Form("hNtracksFb%d",kb));
      tmeanmult[kb]=0x0;
      if(hNtracksFB[kb]){
	hNtracksFB[kb]->SetLineColor(colors[kb]);
	hNtracksFB[kb]->SetLineStyle(lstyl[kb]);
	hNtracksFB[kb]->SetLineWidth(lwid[kb]);
	if(hNtracksFB[kb]->GetMaximum()>maxYaxis) maxYaxis=1.5*hNtracksFB[kb]->GetMaximum();
	legmfb->AddEntry(hNtracksFB[kb],Form("FiltBit %d",kb),"L")->SetTextColor(colors[kb]);
	hMean->SetBinContent(hMean->FindBin(kb),hNtracksFB[kb]->GetMean());
	hMean->SetBinError(hMean->FindBin(kb),hNtracksFB[kb]->GetMeanError());	
	hNtracksFB[kb]->SetStats(0);
	hNtracksFB[kb]->GetXaxis()->SetTitle("N_{tracks}");
	hNtracksFB[kb]->GetYaxis()->SetTitle("kINT7 Events");
	hNtracksFB[kb]->SetTitle("");
	tmeanmult[kb]=new TLatex(kb-0.4,hNtracksFB[kb]->GetMean()+0.02*hNtracksFB0->GetMean(),Form("%.2f",hNtracksFB[kb]->GetMean()));
	tmeanmult[kb]->SetTextFont(43);
	tmeanmult[kb]->SetTextSize(20);
	tmeanmult[kb]->SetTextColor(colors[kb]);
      }
    }
    TCanvas* cmult=new TCanvas("cmult","Track Multipl",1200,500);
    cmult->Divide(2,1);
    cmult->cd(1);
    gPad->SetLogy();
    hNtracksFB[0]->SetMaximum(maxYaxis);
    hNtracksFB[0]->GetYaxis()->SetTitleOffset(1.2);
    hNtracksFB[0]->Draw();
    for(Int_t kb=1; kb<9; kb++){
      if(hNtracksFB[kb]) hNtracksFB[kb]->Draw("same");
    }
    legmfb->Draw();
    cmult->cd(2);
    gStyle->SetPaintTextFormat(".2f");
    hMean->SetStats(0);
    hMean->SetMarkerStyle(20);
    hMean->GetYaxis()->SetTitleOffset(1.2);
    hMean->Draw();
    for(Int_t kb=0; kb<9; kb++){
      if(tmeanmult[kb]){tmeanmult[kb]->Draw();}
    }
    cmult->SaveAs("TrackMultDistPerFilterBit.png");
  }
  TCanvas* cip=new TCanvas("cip","FiltBits",1100,900);
  cip->Divide(2,2);
  cip->cd(1);
  gPad->SetRightMargin(0.13);
  hFilterBits->Draw("colz");
  cip->cd(2);
  gPad->SetLogy();
  TLegend* leg2=new TLegend(0.4,0.5,0.89,0.89);
  leg2->SetMargin(0.3);
  leg2->SetNColumns(2);
  for(Int_t jbit=0; jbit<12; jbit++){
    TString hname=Form("hImpParXYPtMulPionFiltBit%d",jbit);
    TH3F* h=(TH3F*)l->FindObject(hname.Data());
    TH1D* h1=h->ProjectionY(Form("hImpParXYFiltBit%d",jbit));
    h1->SetLineColor(colors[jbit]);
    h1->SetLineStyle(lstyl[jbit]);
    h1->SetLineWidth(lwid[jbit]);
    h1->SetStats(0);
    if(jbit==0){ 
      h1->DrawNormalized();
    }
    else h1->DrawNormalized("same");
    leg2->AddEntry(h1,Form("Filt Bit %d",jbit),"L")->SetTextColor(colors[jbit]);
  }

  for(Int_t jbit=0; jbit<12; jbit++){
    TString hname=Form("hEtaPhiPtFiltBit%d",jbit);
    TH3F* h=(TH3F*)l->FindObject(hname.Data());
    TH1D* hphi=h->ProjectionY(Form("hPhiFiltBit%d",jbit));
    hphi->SetLineColor(colors[jbit]);
    hphi->SetLineStyle(lstyl[jbit]);
    hphi->SetLineWidth(lwid[jbit]);
    hphi->SetStats(0);
    TH1D* hpt1=h->ProjectionZ(Form("hPtFiltBit%d",jbit));
    vecForTrend[2*jbit]=hpt1->GetMean();
    hpt1->SetLineColor(colors[jbit]);
    hpt1->SetLineStyle(lstyl[jbit]);
    hpt1->SetLineWidth(lwid[jbit]);
    hpt1->SetStats(0);
    cip->cd(3);
    if(jbit==0){ 
      hphi->SetMinimum(0.);
      hphi->SetMaximum(1.2*hphi->GetMaximum());
      hphi->Draw();
    }else{ 
      hphi->Draw("same");
    }
    cip->cd(4);
    gPad->SetLogy();
    if(jbit==0){ 
      hpt1->SetMinimum(0.9);
      hpt1->GetXaxis()->SetRangeUser(0.,15.);
      hpt1->SetMaximum(4.*hpt1->GetMaximum());
      hpt1->Draw();
    }else{ 
      hpt1->Draw("same");
    }
    leg2->Draw();
    outfb->cd();
    hphi->Write();
    hpt1->Write();
    curDir->cd();
  }
  plotFileName=Form("FilterBits.%s",outputForm.Data());
  cip->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());


  for(Int_t jb=0; jb<12; jb++){
    TCanvas* ccc=new TCanvas(Form("cfb%d",jb),Form("cfb%d",jb),1500,800);
    ccc->Divide(3,2);
    ccc->cd(1);
    gPad->SetLogz();
    TH2F* hdist1=(TH2F*)l->FindObject(Form("hITScluPtFiltBit%d",jb));
    hdist1->SetStats(0);
    hdist1->SetTitle(Form("Filter bit %d",jb ));
    hdist1->Draw("colz");
    TProfile* hpr1=hdist1->ProfileX(Form("%s-Profile",hdist1->GetName()));
    hpr1->SetMarkerStyle(21);
    hpr1->SetMarkerSize(0.8);
    hpr1->Draw("same");
    ccc->cd(2);
    gPad->SetLogz();
    TH2F* hdist2=(TH2F*)l->FindObject(Form("hSPDcluPtFiltBit%d",jb));
    hdist2->SetTitle(Form("Filter bit %d",jb ));
    hdist2->SetStats(0);
    hdist2->Draw("colz");
    TProfile* hpr2=hdist2->ProfileX(Form("%s-Profile",hdist2->GetName()));
    hpr2->SetMarkerStyle(21);
    hpr2->SetMarkerSize(0.8);
    hpr2->Draw("same");
    ccc->cd(3);
    gPad->SetLogz();
    TH2F* hdist3=(TH2F*)l->FindObject(Form("hTPCcluPtFiltBit%d",jb));
    hdist3->SetTitle(Form("Filter bit %d",jb ));
    hdist3->SetStats(0);
    hdist3->Draw("colz");
    TProfile* hpr3=hdist3->ProfileX(Form("%s-Profile",hdist3->GetName()));
    hpr3->SetMarkerStyle(21);
    hpr3->SetMarkerSize(0.8);
    hpr3->Draw("same");
    ccc->cd(4);
    gPad->SetLogz();
    TH2F* hdist4=(TH2F*)l->FindObject(Form("hTPCcrrowsPtFiltBit%d",jb));
    hdist4->SetTitle(Form("Filter bit %d",jb ));
    hdist4->SetStats(0);
    hdist4->Draw("colz");
    TProfile* hpr4=hdist4->ProfileX(Form("%s-Profile",hdist4->GetName()));
    hpr4->SetMarkerStyle(21);
    hpr4->SetMarkerSize(0.8);
    hpr4->Draw("same");
    ccc->cd(5);
    gPad->SetLogz();
    TH2F* hdist5=(TH2F*)l->FindObject(Form("hTPCCrowOverFindPtFiltBit%d",jb));
    hdist5->SetTitle(Form("Filter bit %d",jb ));
    hdist5->SetStats(0);
    hdist5->Draw("colz");
    TProfile* hpr5=hdist5->ProfileX(Form("%s-Profile",hdist5->GetName()));
    hpr5->SetMarkerStyle(21);
    hpr5->SetMarkerSize(0.8);
    hpr5->Draw("same");
    ccc->cd(6);
    gPad->SetLogz();
    TH2F* hdist6=(TH2F*)l->FindObject(Form("hTPCChi2ndfPtFiltBit%d",jb));
    if(!hdist6) hdist6=(TH2F*)l->FindObject(Form("hTPCChi2clusPtFiltBit%d",jb));
    hdist6->SetTitle(Form("Filter bit %d",jb ));
    hdist6->SetStats(0);
    hdist6->Draw("colz");
    TProfile* hpr6=hdist6->ProfileX(Form("%s-Profile",hdist6->GetName()));
    hpr6->SetMarkerStyle(21);
    hpr6->SetMarkerSize(0.8);
    hpr6->Draw("same");
    outfb->cd();
    hpr1->Write();
    hpr2->Write();
    hpr3->Write();
    hpr4->Write();
    hpr5->Write();
    hpr6->Write();
    curDir->cd();

    plotFileName=Form("VarDistFiltBit%d.%s",jb,outputForm.Data());
    ccc->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());
  }
  outfb->Close();
  
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

  // K0s histos vs. radius
  z1=hInvMassK0s3d->GetZaxis()->FindBin(2.99);
  TH1D* hInvMassK0sR1=hInvMassK0s3d->ProjectionX("hInvMassK0s1dR1",0,-1,1,z1);
  z2=hInvMassK0s3d->GetZaxis()->FindBin(5.99);
  TH1D* hInvMassK0sR2=hInvMassK0s3d->ProjectionX("hInvMassK0s1dR2",0,-1,z1+1,z2);
  z3=hInvMassK0s3d->GetZaxis()->FindBin(8.01);
  z4=hInvMassK0s3d->GetZaxis()->FindBin(22.99);
  TH1D* hInvMassK0sR3=hInvMassK0s3d->ProjectionX("hInvMassK0s1dR3",0,-1,z3,z4);
  z5=hInvMassK0s3d->GetZaxis()->FindBin(28.01);
  z6=hInvMassK0s3d->GetZaxis()->FindBin(42.99);
  TH1D* hInvMassK0sR4=hInvMassK0s3d->ProjectionX("hInvMassK0s1dR4",0,-1,z5,z6);


  TF1* fmassk0=new TF1("fmassk0","[0]+[1]*x+[2]*x*x+[3]/sqrt(2.*TMath::Pi())/[5]*TMath::Exp(-0.5*(x-[4])*(x-[4])/[5]/[5])",0.44,0.56);
  fmassk0->SetLineWidth(2);
  fmassk0->SetLineColor(kMagenta+1);

  TF1* fmassL=new TF1("fmassL","[0]+[1]*x+[2]*x*x+[3]/sqrt(2.*TMath::Pi())/[5]*TMath::Exp(-0.5*(x-[4])*(x-[4])/[5]/[5])",1.10,1.13);
  fmassL->SetLineWidth(2);
  fmassL->SetLineColor(kRed+1);

  TCanvas* cv0=new TCanvas("cv0","V0s - all pt",1500,500);
  cv0->Divide(3,1);
  cv0->cd(1);
  hInvMassK0s->Draw();
  InitFuncAndFit(hInvMassK0s,fmassk0,kTRUE,isMC);
  cv0->cd(2);
  hInvMassLambda->Draw();
  InitFuncAndFit(hInvMassLambda,fmassL,kFALSE,isMC);
  cv0->cd(3);
  hInvMassAntiLambda->Draw();
  InitFuncAndFit(hInvMassAntiLambda,fmassL,kFALSE,isMC);
  plotFileName=Form("MassSpectraV0-integrated.%s",outputForm.Data());
  cv0->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  TCanvas* clam=new TCanvas("clam","Lambda vs R",1400,900);
  clam->Divide(2,2);
  clam->cd(1);
  hInvMassLambdaR1->Draw();
  InitFuncAndFit(hInvMassLambdaR1,fmassL,kFALSE,isMC);
  TLatex* tr1=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(1)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z1+1))));
  tr1->SetNDC();
  tr1->SetTextFont(43);
  tr1->SetTextSize(26);
  tr1->Draw();
  clam->cd(2);
  hInvMassLambdaR2->Draw();
  InitFuncAndFit(hInvMassLambdaR2,fmassL,kFALSE,isMC);
  TLatex* tr2=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z1+1)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z2+1))));
  tr2->SetNDC();
  tr2->SetTextFont(43);
  tr2->SetTextSize(26);
  tr2->Draw();
  clam->cd(3);
  hInvMassLambdaR3->Draw();
  InitFuncAndFit(hInvMassLambdaR3,fmassL,kFALSE,isMC);
  TLatex* tr3=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z3)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z4+1))));
  tr3->SetNDC();
  tr3->SetTextFont(43);
  tr3->SetTextSize(26);
  tr3->Draw();
  clam->cd(4);
  hInvMassLambdaR4->Draw();
  InitFuncAndFit(hInvMassLambdaR4,fmassL,kFALSE,isMC);
  TLatex* tr4=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z5)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z6+1))));
  tr4->SetNDC();
  tr4->SetTextFont(43);
  tr4->SetTextSize(26);
  tr4->Draw();
  plotFileName=Form("Lambda-MassSpectra-VsR.%s",outputForm.Data());
  clam->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());


  TCanvas* ck0=new TCanvas("ck0","K0s vs. pt",1400,900);
  ck0->Divide(2,2);
  ck0->cd(1);
  hInvMassK0sP1->Draw();
  InitFuncAndFit(hInvMassK0sP1,fmassk0,kTRUE,isMC);
  TLatex* tp1=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(1),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p1+1)));
  tp1->SetNDC();
  tp1->SetTextFont(43);
  tp1->SetTextSize(26);
  tp1->Draw();
  ck0->cd(2);
  hInvMassK0sP2->Draw();
  InitFuncAndFit(hInvMassK0sP2,fmassk0,kTRUE,isMC);
  TLatex* tp2=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p1+1),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p2+1)));
  tp2->SetNDC();
  tp2->SetTextFont(43);
  tp2->SetTextSize(26);
  tp2->Draw();
  ck0->cd(3);
  hInvMassK0sP3->Draw();
  InitFuncAndFit(hInvMassK0sP3,fmassk0,kTRUE,isMC);
  TLatex* tp3=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p2+1),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p3+1)));
  tp3->SetNDC();
  tp3->SetTextFont(43);
  tp3->SetTextSize(26);
  tp3->Draw();
  ck0->cd(4);
  hInvMassK0sP4->Draw();
  InitFuncAndFit(hInvMassK0sP4,fmassk0,kTRUE,isMC);
  TLatex* tp4=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p4),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p5+1)));
  tp4->SetNDC();
  tp4->SetTextFont(43);
  tp4->SetTextSize(26);
  tp4->Draw();
  plotFileName=Form("K0s-MassSpectra-VsPt.%s",outputForm.Data());
  ck0->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  TCanvas* ck0r=new TCanvas("ck0r","K0s vs R",1400,900);
  ck0r->Divide(2,2);
  ck0r->cd(1);
  hInvMassK0sR1->Draw();
  InitFuncAndFit(hInvMassK0sR1,fmassk0,kTRUE,isMC);
  TLatex* tr1k=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassK0s3d->GetZaxis()->GetBinLowEdge(1)),TMath::Nint(hInvMassK0s3d->GetZaxis()->GetBinLowEdge(z1+1))));
  tr1k->SetNDC();
  tr1k->SetTextFont(43);
  tr1k->SetTextSize(26);
  tr1k->Draw();
  ck0r->cd(2);
  hInvMassK0sR2->Draw();
  InitFuncAndFit(hInvMassK0sR2,fmassk0,kTRUE,isMC);
  TLatex* tr2k=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassK0s3d->GetZaxis()->GetBinLowEdge(z1+1)),TMath::Nint(hInvMassK0s3d->GetZaxis()->GetBinLowEdge(z2+1))));
  tr2k->SetNDC();
  tr2k->SetTextFont(43);
  tr2k->SetTextSize(26);
  tr2k->Draw();
  ck0r->cd(3);
  hInvMassK0sR3->Draw();
  InitFuncAndFit(hInvMassK0sR3,fmassk0,kTRUE,isMC);
  TLatex* tr3k=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassK0s3d->GetZaxis()->GetBinLowEdge(z3)),TMath::Nint(hInvMassK0s3d->GetZaxis()->GetBinLowEdge(z4+1))));
  tr3k->SetNDC();
  tr3k->SetTextFont(43);
  tr3k->SetTextSize(26);
  tr3k->Draw();
  ck0r->cd(4);
  hInvMassK0sR4->Draw();
  InitFuncAndFit(hInvMassK0sR4,fmassk0,kTRUE,isMC);
  TLatex* tr4k=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassK0s3d->GetZaxis()->GetBinLowEdge(z5)),TMath::Nint(hInvMassK0s3d->GetZaxis()->GetBinLowEdge(z6+1))));
  tr4k->SetNDC();
  tr4k->SetTextFont(43);
  tr4k->SetTextSize(26);
  tr4k->Draw();
  plotFileName=Form("K0s-MassSpectra-VsR.%s",outputForm.Data());
  ck0r->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  // K0 pt resolution vs. pt
  const Int_t nPtBinsK0=10;
  Double_t ptbinlimsK0[nPtBinsK0+1]={0.,0.4,0.8,1.2,2.0,3.,4.,5.,6.,8.,10.};
  TH1F* hSigmaK0AllR=new TH1F("hSigmaK0AllR"," ; p_{T} (GeV/c) ; #sigma_{K0} (MeV/c^{2})",nPtBinsK0,ptbinlimsK0);
  TH1F* hSigmaK0R4=new TH1F("hSigmaK0R4"," ; p_{T} (GeV/c) ; #sigma_{K0} (MeV/c^{2})",nPtBinsK0,ptbinlimsK0);
  TH1F* hMassK0AllR=new TH1F("hMassK0AllR"," ; p_{T} (GeV/c) ; #mu_{K0} (GeV/c^{2})",nPtBinsK0,ptbinlimsK0);
  TH1F* hYieldK0AllR=new TH1F("hYieldK0AllR"," ; p_{T} (GeV/c) ; N_{K0}/event",nPtBinsK0,ptbinlimsK0);

  TCanvas* ctmpk0=new TCanvas("ctmpk0","K0s vs. pt R<4",1600,900);
  ctmpk0->Divide(5,4);
  for(Int_t ipt=0; ipt<nPtBinsK0; ipt++){
    Int_t pfine1=hInvMassK0s3d->GetYaxis()->FindBin(ptbinlimsK0[ipt]+0.001);
    Int_t pfine2=hInvMassK0s3d->GetYaxis()->FindBin(ptbinlimsK0[ipt+1]-0.001);
    Int_t r1=hInvMassK0s3d->GetZaxis()->FindBin(0.001);
    Int_t r2=hInvMassK0s3d->GetZaxis()->FindBin(3.999);
    TH1D* hTmpInvMassK0sR4=hInvMassK0s3d->ProjectionX(Form("hInvMassK0sR4PtFine%d",ipt),pfine1,pfine2,r1,r2);
    TH1D* hTmpInvMassK0sAllR=hInvMassK0s3d->ProjectionX(Form("hInvMassK0sAllRPtFine%d",ipt),pfine1,pfine2,0,-1);
    ctmpk0->cd(ipt+1);
    hTmpInvMassK0sR4->Draw();
    InitFuncAndFit(hTmpInvMassK0sR4,fmassk0,kTRUE,isMC);
    TLatex* tpfine=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(pfine1),hInvMassK0s3d->GetYaxis()->GetBinUpEdge(pfine2)));
    tpfine->SetNDC();
    tpfine->SetTextFont(43);
    tpfine->SetTextSize(24);
    tpfine->Draw();
    hSigmaK0R4->SetBinContent(ipt+1,fmassk0->GetParameter(5)*1000.);
    hSigmaK0R4->SetBinError(ipt+1,fmassk0->GetParError(5)*1000.);
    ctmpk0->cd(ipt+11);
    hTmpInvMassK0sAllR->Draw();
    InitFuncAndFit(hTmpInvMassK0sAllR,fmassk0,kTRUE,isMC);
    tpfine->Draw();
    hMassK0AllR->SetBinContent(ipt+1,fmassk0->GetParameter(4));
    hMassK0AllR->SetBinError(ipt+1,fmassk0->GetParError(4));
    hSigmaK0AllR->SetBinContent(ipt+1,fmassk0->GetParameter(5)*1000.);
    hSigmaK0AllR->SetBinError(ipt+1,fmassk0->GetParError(5)*1000.);
    Double_t yield=fmassk0->GetParameter(3)/hTmpInvMassK0sAllR->GetBinWidth(1)/nSelectedEvents;
    Double_t eyield=fmassk0->GetParError(3)/hTmpInvMassK0sAllR->GetBinWidth(1)/nSelectedEvents;
    hYieldK0AllR->SetBinContent(ipt+1,yield);
    hYieldK0AllR->SetBinError(ipt+1,eyield);
  }

  TCanvas* cK0signal=new TCanvas("cK0signal","K0 width and yield vs pt",1600,500);
  cK0signal->Divide(3,1);
  cK0signal->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hMassK0AllR->SetMinimum(0.495);
  hMassK0AllR->SetMaximum(0.500);
  hMassK0AllR->SetStats(0);
  hMassK0AllR->SetMarkerStyle(20);
  hMassK0AllR->SetLineWidth(2);
  hMassK0AllR->Draw();
  cK0signal->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hSigmaK0AllR->SetMinimum(0);
  hSigmaK0AllR->SetMaximum(10);
  hSigmaK0AllR->SetStats(0);
  hSigmaK0AllR->SetMarkerStyle(20);
  hSigmaK0AllR->SetLineWidth(2);
  hSigmaK0R4->SetMarkerStyle(25);
  hSigmaK0R4->SetMarkerColor(kRed+1);
  hSigmaK0R4->SetLineColor(kRed+1);
  hSigmaK0R4->SetLineWidth(2);
  hSigmaK0AllR->Draw();
  hSigmaK0R4->Draw("same");
  TLegend* lk=new TLegend(0.18,0.18,0.5,0.3);
  lk->AddEntry(hSigmaK0AllR,"All decay radii","P")->SetTextColor(hSigmaK0AllR->GetMarkerColor());
  lk->AddEntry(hSigmaK0R4,"R < 4 cm","P")->SetTextColor(hSigmaK0R4->GetMarkerColor());
  lk->Draw();
  cK0signal->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hYieldK0AllR->SetStats(0);
  hYieldK0AllR->SetMarkerStyle(20);
  hYieldK0AllR->SetLineWidth(2);
  hYieldK0AllR->Draw();
  plotFileName=Form("K0s-SignalVsPt.%s",outputForm.Data());
  cK0signal->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  trtree->Fill();

  if(runNumber>0){
    TFile* foutfile=new TFile("trendingAODtracks.root","recreate");
    trtree->Write();
    TDirectory* outdir=foutfile->mkdir(df->GetName());
    outdir->cd();
    l->Write(l->GetName(),1);
    foutfile->Close();
    delete foutfile;
  }

  TFile* foutk0=new TFile("SigmaK0s.root","recreate");
  hSigmaK0AllR->Write();
  hSigmaK0R4->Write();
  foutk0->Close();
  delete foutk0;
  
  if(outputForm=="pdf") gSystem->Exec(Form("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=PlotsAODTrackQA.pdf %s",pdfFileNames.Data()));

  printf("SUMMARY:\n");
  printf("Number of events used in the plots = %d\n",nSelectedEvents);

}

void InitFuncAndFit(TH1D* hm, TF1* fmass, Bool_t isK0s, Bool_t isMC){

  // first estimate of background
  Double_t cntpeak,expSigma;
  if(isK0s){
    TF1* ffb = new TF1("fp2bkgk0",fp2bkgk0,0.44,0.56,3);
    hm->Fit("fp2bkgk0","R");
    for(Int_t k=0; k<3; k++) fmass->SetParameter(k,ffb->GetParameter(k));
    cntpeak=hm->GetBinContent(hm->FindBin(0.498))-ffb->Integral(0.498-0.5*hm->GetBinWidth(1),0.498+hm->GetBinWidth(1));
    delete ffb;
    expSigma=0.004;
  }else{
    fmass->SetParameter(0,hm->GetBinContent(hm->FindBin(1.10)));
    fmass->SetParameter(1,0.);
    fmass->FixParameter(2,0.);
    cntpeak=hm->GetBinContent(hm->FindBin(1.116))-hm->GetBinContent(hm->FindBin(1.14));
    expSigma=0.0015;
  }
  //  fmass->SetParLimits(1,-99999999999,0.);
  fmass->SetParameter(3,cntpeak*TMath::Sqrt(2*TMath::Pi())*expSigma);
  //  fmass->SetParLimits(3,0.,9999999999.);
  if(isK0s){
    fmass->SetParameter(4,0.5);
    fmass->SetParLimits(4,0.49,0.51);
    fmass->SetParameter(5,0.002);
    fmass->SetParLimits(5,0.0006,0.02);
  }else{
    fmass->SetParameter(4,1.116);
    fmass->SetParLimits(4,1.11,1.12);
    fmass->SetParameter(5,0.0015);
    fmass->SetParLimits(5,0.0006,0.003);
  }
  if(isMC){
    fmass->FixParameter(0,0.);
    fmass->FixParameter(1,0.);
    fmass->FixParameter(2,0.);
  }

  if(isMC)hm->Fit(fmass,"R");
  else hm->Fit(fmass,"RL");
  TLatex* t1=new TLatex(0.14,0.8,Form("Mean = %.3f+-%.3f GeV/c^{2}",fmass->GetParameter(4),fmass->GetParError(4)));
  t1->SetTextSize(0.04);
  t1->SetNDC();
  t1->Draw();
  TLatex* t2=new TLatex(0.14,0.75,Form("Sigma = %.2f+-%.2f MeV/c^{2}",fmass->GetParameter(5)*1000.,fmass->GetParError(5)*1000.));
  t2->SetNDC();
  t2->SetTextSize(0.04);
  t2->Draw();
  TLatex* t3=new TLatex(0.14,0.7,Form("Yield = %.0f+-%.0f",fmass->GetParameter(3)/hm->GetBinWidth(1),fmass->GetParError(3)/hm->GetBinWidth(1)));
  t3->SetNDC();
  t3->SetTextSize(0.04);
  t3->Draw();

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

TH1D* ComputeMatchEff(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle){
  if(hnumer->GetSumw2N()==0) hnumer->Sumw2();
  if(hdenom->GetSumw2N()==0) hdenom->Sumw2();
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

void FillMeanAndRms(TH2F* h2d, TGraphErrors* gMean, TGraphErrors* gRms){
  Int_t jpt=0;
  Int_t jptr=0;
  for(Int_t j=1; j<=h2d->GetNbinsX(); j++){
    Double_t pt=h2d->GetXaxis()->GetBinCenter(j);
    Double_t ept=h2d->GetXaxis()->GetBinWidth(j)/2.;
    Int_t ji=j;
    Int_t jf=j;
    if(pt>5 && pt<=10){
      ji=j;
      jf=j+1;
      j=jf;
      pt=0.5*(h2d->GetXaxis()->GetBinLowEdge(ji)+h2d->GetXaxis()->GetBinUpEdge(jf));
      ept=pt-h2d->GetXaxis()->GetBinLowEdge(ji);
    }else if(pt>10 && pt<=16){
      ji=j;
      jf=j+3;
      j=jf;
      pt=0.5*(h2d->GetXaxis()->GetBinLowEdge(ji)+h2d->GetXaxis()->GetBinUpEdge(jf));
      ept=pt-h2d->GetXaxis()->GetBinLowEdge(ji);
    }else if(pt>16 && pt<=20){
      ji=j;
      jf=j+7;
      j=jf;
      pt=0.5*(h2d->GetXaxis()->GetBinLowEdge(ji)+h2d->GetXaxis()->GetBinUpEdge(jf));
      ept=pt-h2d->GetXaxis()->GetBinLowEdge(ji);
     }else if(pt>20){
      ji=j;
      jf=j+19;
      j=jf;
      pt=0.5*(h2d->GetXaxis()->GetBinLowEdge(ji)+h2d->GetXaxis()->GetBinUpEdge(jf));
      ept=pt-h2d->GetXaxis()->GetBinLowEdge(ji);
    }
    TH1D* htmp=h2d->ProjectionY("htmp",ji,jf);
    if(htmp->Integral()>40){
      htmp->Fit("gaus","Q0");
      TF1* fg=(TF1*)htmp->GetListOfFunctions()->FindObject("gaus");      
      Double_t m=fg->GetParameter(1);//htmp->GetMean();
      Double_t em=fg->GetParError(1);//htmp->GetMeanError();
      Double_t r=fg->GetParameter(2);//htmp->GetRMS();
      Double_t er=fg->GetParError(2);//=htmp->GetRMSError();
      if(er/r<0.35){
	gMean->SetPoint(jpt,pt,m);
	gMean->SetPointError(jpt,ept,em);
	++jpt;
	gRms->SetPoint(jptr,pt,r);
	gRms->SetPointError(jptr,ept,er);
	++jptr;
      }
      delete htmp;
    }
  }
}
 

Double_t fp2bkgk0(Double_t *x, Double_t *par){
  if (x[0] > 0.47 && x[0] < 0.53) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}
