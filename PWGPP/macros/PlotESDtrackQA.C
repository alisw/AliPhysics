#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH3F.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TSystem.h>
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
void InitFuncAndFit(TH1D* hm, TF1* fmass, Bool_t isK0s, Bool_t isMC=kFALSE);
Double_t fp2bkgk0(Double_t *x, Double_t *par);

Double_t maxPtHypoPlots=5.;
Double_t maxPtMEPlots=20.;

void PlotESDtrackQA(TString filename="QAresults.root", TString suffix="QA", Int_t runNumber=-1, TString outputForm="png"){


  TString pdfFileNames="";
  TString plotFileName="";
  
  TTree* trtree=new TTree("trending","tree of trending variables");
  trtree->Branch("nrun",&runNumber,"nrun/I");
  Bool_t isMC=kFALSE; // set automatically based on histos filled

  // match eff variables
  Double_t ptForTrend[3]={0.35,1.,4.};
  Double_t ptForBadHyp[3]={0.35,0.6,0.9};
  Double_t ptForResol[4]={0.35,1.,4.,10.};

  Float_t vecMatchEff[24];
  Float_t vecErrMatchEff[24];
  for(Int_t itof=0; itof<2; itof++){
    TString tof="";
    if(itof==1) tof="TOFbc";
    for(Int_t ispd=0; ispd<2; ispd++){
      TString spd="";
      if(ispd==1) spd="SPDany";
      for(Int_t isid=0; isid<2; isid++){
	TString side="Pos";
	if(isid==1) side="Neg";
	for(Int_t ipt=0; ipt<3; ipt++){
	  TString bName=Form("MatchEff%sPt%dEta%s%s",spd.Data(),TMath::Nint(ptForTrend[ipt]*1000.),side.Data(),tof.Data());
	  Int_t index=itof*12+ispd*6+isid*3+ipt;
	  TString errbName=Form("err%s",bName.Data());
	  trtree->Branch(bName.Data(),&vecMatchEff[index],Form("%s/F",bName.Data()));
	  trtree->Branch(errbName.Data(),&vecErrMatchEff[index],Form("%s/F",errbName.Data()));
	}
      }
    }
  }
  Float_t vecPosNeg[12];
  Float_t vecErrPosNeg[12];
  for(Int_t ispd=0; ispd<2; ispd++){
    TString spd="";
    if(ispd==1) spd="SPDany";
    for(Int_t isid=0; isid<2; isid++){
      TString side="Pos";
      if(isid==1) side="Neg";
      for(Int_t ipt=0; ipt<3; ipt++){
	TString bName=Form("PosNegCharge%sPt%dEta%s",spd.Data(),TMath::Nint(ptForTrend[ipt]*1000.),side.Data());
	Int_t index=ispd*6+isid*3+ipt;
	trtree->Branch(bName.Data(),&vecPosNeg[index],Form("%s/F",bName.Data()));
	TString errbName=Form("err%s",bName.Data());
	trtree->Branch(errbName.Data(),&vecErrPosNeg[index],Form("%s/F",errbName.Data()));
      }
    }
  }


  TFile* f=new TFile(filename.Data());
  TDirectoryFile* df=(TDirectoryFile*)f->Get("CheckESDTracks");
  if(!df){
    printf("Directory CheckESDTracks not found in file %s\n",filename.Data());
    f->ls();
    return;
  }
  TList* l=(TList*)df->Get(Form("clistCheckESDTracks%s",suffix.Data()));
  if(!l){
    printf("TList clistCheckESDTracks%s not found in file %s\n",suffix.Data(),filename.Data());
    df->ls();
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
  TH3F* hEtaPhiPtTPCselTPCpt=(TH3F*)l->FindObject("hEtaPhiPtInnerTPCsel");
  TH3F* hEtaPhiPtTPCselITSrefTPCpt=(TH3F*)l->FindObject("hEtaPhiPtInnerTPCselITSref");
  TH3F* hEtaPhiPtTPCselSPDanyTPCpt=(TH3F*)l->FindObject("hEtaPhiPtInnerTPCselSPDany");
  TH3F* hEtaPhiPtTPCselTPCptTOFbc=(TH3F*)l->FindObject("hEtaPhiPtInnerTPCselTOFbc");
  TH3F* hEtaPhiPtTPCselITSrefTPCptTOFbc=(TH3F*)l->FindObject("hEtaPhiPtInnerTPCselITSrefTOFbc");
  TH3F* hEtaPhiPtTPCselSPDanyTPCptTOFbc=(TH3F*)l->FindObject("hEtaPhiPtInnerTPCselSPDanyTOFbc");
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

  TH1D* hPtEtaNegTPCselTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionZ("hPtEtaNegTPCselTOFbc",etamin,eta0m);
  TH1D* hPtEtaPosTPCselTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionZ("hPtEtaPosTPCselTOFbc",eta0p,etamax);
  TH1D* hPtEtaNegTPCselITSrefTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionZ("hPtEtaNegTPCselITSrefTOFbc",etamin,eta0m);
  TH1D* hPtEtaPosTPCselITSrefTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionZ("hPtEtaPosTPCselITSrefTOFbc",eta0p,etamax);
  TH1D* hPtEtaNegTPCselSPDanyTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionZ("hPtEtaNegTPCselSPDanyTOFbc",etamin,eta0m);
  TH1D* hPtEtaPosTPCselSPDanyTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionZ("hPtEtaPosTPCselSPDanyTOFbc",eta0p,etamax);
  TH1D* hPhiEtaNegTPCselLowPtTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionY("hPhiEtaNegTPCselLowPtTOFbc",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselLowPtTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionY("hPhiEtaPosTPCselLowPtTOFbc",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselHighPtTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionY("hPhiEtaNegTPCselHighPtTOFbc",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselHighPtTOFbc=hEtaPhiPtTPCselTOFbc->ProjectionY("hPhiEtaPosTPCselHighPtTOFbc",eta0p,etamax,ptone,ptten);
  TH1D* hPhiEtaNegTPCselITSrefLowPtTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionY("hPhiEtaNegTPCselITSrefLowPtTOFbc",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselITSrefLowPtTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionY("hPhiEtaPosTPCselITSrefLowPtTOFbc",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselITSrefHighPtTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionY("hPhiEtaNegTPCselITSrefHighPtTOFbc",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselITSrefHighPtTOFbc=hEtaPhiPtTPCselITSrefTOFbc->ProjectionY("hPhiEtaPosTPCselITSrefHighPtTOFbc",eta0p,etamax,ptone,ptten);
  TH1D* hPhiEtaNegTPCselSPDanyLowPtTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionY("hPhiEtaNegTPCselSPDanyLowPtTOFbc",etamin,eta0m,ptzero4,ptzero7);
  TH1D* hPhiEtaPosTPCselSPDanyLowPtTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionY("hPhiEtaPosTPCselSPDanyLowPtTOFbc",eta0p,etamax,ptzero4,ptzero7);
  TH1D* hPhiEtaNegTPCselSPDanyHighPtTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionY("hPhiEtaNegTPCselSPDanyHighPtTOFbc",etamin,eta0m,ptone,ptten);
  TH1D* hPhiEtaPosTPCselSPDanyHighPtTOFbc=hEtaPhiPtTPCselSPDanyTOFbc->ProjectionY("hPhiEtaPosTPCselSPDanyHighPtTOFbc",eta0p,etamax,ptone,ptten);

  TH1D* hPtEtaNegTPCselTPCpt=hEtaPhiPtTPCselTPCpt->ProjectionZ("hPtEtaNegTPCselTPCpt",etamin,eta0m);
  TH1D* hPtEtaPosTPCselTPCpt=hEtaPhiPtTPCselTPCpt->ProjectionZ("hPtEtaPosTPCselTPCpt",eta0p,etamax);
  TH1D* hPtEtaNegTPCselITSrefTPCpt=hEtaPhiPtTPCselITSrefTPCpt->ProjectionZ("hPtEtaNegTPCselITSrefTPCpt",etamin,eta0m);
  TH1D* hPtEtaPosTPCselITSrefTPCpt=hEtaPhiPtTPCselITSrefTPCpt->ProjectionZ("hPtEtaPosTPCselITSrefTPCpt",eta0p,etamax);
  TH1D* hPtEtaNegTPCselSPDanyTPCpt=hEtaPhiPtTPCselSPDanyTPCpt->ProjectionZ("hPtEtaNegTPCselSPDanyTPCpt",etamin,eta0m);
  TH1D* hPtEtaPosTPCselSPDanyTPCpt=hEtaPhiPtTPCselSPDanyTPCpt->ProjectionZ("hPtEtaPosTPCselSPDanyTPCpt",eta0p,etamax);

  TH1D* hPtEtaNegTPCselTPCptTOFbc=hEtaPhiPtTPCselTPCptTOFbc->ProjectionZ("hPtEtaNegTPCselTPCptTOFbc",etamin,eta0m);
  TH1D* hPtEtaPosTPCselTPCptTOFbc=hEtaPhiPtTPCselTPCptTOFbc->ProjectionZ("hPtEtaPosTPCselTPCptTOFbc",eta0p,etamax);
  TH1D* hPtEtaNegTPCselITSrefTPCptTOFbc=hEtaPhiPtTPCselITSrefTPCptTOFbc->ProjectionZ("hPtEtaNegTPCselITSrefTPCptTOFbc",etamin,eta0m);
  TH1D* hPtEtaPosTPCselITSrefTPCptTOFbc=hEtaPhiPtTPCselITSrefTPCptTOFbc->ProjectionZ("hPtEtaPosTPCselITSrefTPCptTOFbc",eta0p,etamax);
  TH1D* hPtEtaNegTPCselSPDanyTPCptTOFbc=hEtaPhiPtTPCselSPDanyTPCptTOFbc->ProjectionZ("hPtEtaNegTPCselSPDanyTPCptTOFbc",etamin,eta0m);
  TH1D* hPtEtaPosTPCselSPDanyTPCptTOFbc=hEtaPhiPtTPCselSPDanyTPCptTOFbc->ProjectionZ("hPtEtaPosTPCselSPDanyTPCptTOFbc",eta0p,etamax);

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
  hPtEtaNegTPCselSPDany->SetTitle("p_{T} tracks - TPC cuts, SPDany - #eta<0");
  hPtEtaPosTPCselSPDany->SetTitle("p_{T} tracks - TPC cuts, SPDany - #eta>0");
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
  hPtEtaNegTPCselTPCpt->Sumw2();
  hPtEtaNegTPCselITSrefTPCpt->Sumw2();
  hPtEtaNegTPCselSPDanyTPCpt->Sumw2();
  hPtEtaPosTPCselTPCpt->Sumw2();
  hPtEtaPosTPCselITSrefTPCpt->Sumw2();
  hPtEtaPosTPCselSPDanyTPCpt->Sumw2();
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

  TString partNames[9]={"Elec","Muon","Pion","Kaon","Proton","Deuteron","Triton","He3","Alpha"};
  TH2F* hdEdxVsPTPCsel[9];
  TH2F* hdEdxVsPTPCselITSref[9];
  TH2F* hdEdxVsPTPCselAll=0x0;
  TCanvas* cdedxa=new TCanvas("cdedxall","dEdx Vs. hypo",1500,700);
  cdedxa->Divide(2,1);
  TLegend * legtrhyp=new TLegend(0.6,0.5,0.89,0.89);
  legtrhyp->SetHeader("Mass Hypo in tracking");
  for(Int_t jsp=0; jsp<9; jsp++){ 
    hdEdxVsPTPCsel[jsp]=(TH2F*)l->FindObject(Form("hdEdxVsPTPCsel%s",partNames[jsp].Data()));
    hdEdxVsPTPCsel[jsp]->GetXaxis()->SetTitle("p_{TPC} (GeV/c)");
    hdEdxVsPTPCsel[jsp]->GetYaxis()->SetTitle("TPC dE/dx");
    hdEdxVsPTPCsel[jsp]->SetTitle(Form("Tracked with %s mass hypothesis - TPC cuts",partNames[jsp].Data()));
    hdEdxVsPTPCsel[jsp]->SetStats(0);
    hdEdxVsPTPCselITSref[jsp]=(TH2F*)l->FindObject(Form("hdEdxVsPTPCselITSref%s",partNames[jsp].Data()));
    hdEdxVsPTPCselITSref[jsp]->GetXaxis()->SetTitle("p_{TPC} (GeV/c)");
    hdEdxVsPTPCselITSref[jsp]->GetYaxis()->SetTitle("TPC dE/dx");
    hdEdxVsPTPCselITSref[jsp]->SetTitle(Form("Tracked with %s mass hypothesis - TPC cuts+ITS refit",partNames[jsp].Data()));
    hdEdxVsPTPCselITSref[jsp]->SetStats(0);
    if(jsp==0){
      hdEdxVsPTPCselAll=(TH2F*)hdEdxVsPTPCsel[0]->Clone("hdEdxVsPTPCselAll");
      hdEdxVsPTPCselAll->SetTitle("All mass hypotheses");
    }else{
      if(hdEdxVsPTPCselAll) hdEdxVsPTPCselAll->Add(hdEdxVsPTPCsel[jsp]);
    }

    TCanvas* cdedx=new TCanvas(Form("cdedx%s",partNames[jsp].Data()),Form("dEdx Hypo %s",partNames[jsp].Data()),1500,700);
    cdedx->Divide(2,1);
    cdedx->cd(1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    hdEdxVsPTPCsel[jsp]->GetYaxis()->SetTitleOffset(1.3);
    hdEdxVsPTPCsel[jsp]->GetXaxis()->SetTitleOffset(1.1);
    hdEdxVsPTPCsel[jsp]->Draw("colz");
    cdedx->cd(2);
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    hdEdxVsPTPCselITSref[jsp]->GetYaxis()->SetTitleOffset(1.3);
    hdEdxVsPTPCselITSref[jsp]->GetXaxis()->SetTitleOffset(1.1);
    hdEdxVsPTPCselITSref[jsp]->Draw("colz");
    plotFileName=Form("dEdx-TrackRecoWith%sHypo.%s",partNames[jsp].Data(),outputForm.Data());
    cdedx->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());
    
    TH2F* htmp1=(TH2F*)hdEdxVsPTPCsel[jsp]->Clone(Form("%s_1",hdEdxVsPTPCselITSref[jsp]->GetName()));
    htmp1->SetMarkerStyle(1);
    htmp1->SetMarkerColor(jsp+1);
    TH2F* htmp2=(TH2F*)hdEdxVsPTPCselITSref[jsp]->Clone(Form("%s_2",hdEdxVsPTPCselITSref[jsp]->GetName()));
    htmp2->SetMarkerStyle(1);
    htmp2->SetMarkerColor(jsp+1);
    cdedxa->cd(1);
    if(jsp==0) htmp1->Draw("P");
    else htmp1->Draw("PSAME");
    cdedxa->cd(2);
    if(jsp==0) htmp2->Draw("P");
    else htmp2->Draw("PSAME");
    legtrhyp->AddEntry(htmp1,partNames[jsp].Data(),"")->SetTextColor(htmp1->GetMarkerColor());
  }
  cdedxa->cd(1);
  legtrhyp->Draw();
  cdedxa->cd(2);
  legtrhyp->Draw();
  plotFileName=Form("dEdxVsTrackRecoHypo.%s",outputForm.Data());
  cdedxa->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());


  TH1D* hMatchEffVsPtNegEta=ComputeMatchEff(hPtEtaNegTPCselITSref,hPtEtaNegTPCsel,"hMatchEffVsPtNegEta",1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEta=ComputeMatchEff(hPtEtaPosTPCselITSref,hPtEtaPosTPCsel,"hMatchEffVsPtPosEta",1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtNegEtaSPDany=ComputeMatchEff(hPtEtaNegTPCselSPDany,hPtEtaNegTPCsel,"hMatchEffVsPtNegEtaSPDAny",kBlue-7,33,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaSPDany=ComputeMatchEff(hPtEtaPosTPCselSPDany,hPtEtaPosTPCsel,"hMatchEffVsPtPosEtaSPDAny",kBlue-7,33,"p_{T} (GeV/c)");

  TH1D* hMatchEffVsPtNegEtaTOFbc=ComputeMatchEff(hPtEtaNegTPCselITSrefTOFbc,hPtEtaNegTPCselTOFbc,"hMatchEffVsPtNegEtaTOFbc",kRed+1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaTOFbc=ComputeMatchEff(hPtEtaPosTPCselITSrefTOFbc,hPtEtaPosTPCselTOFbc,"hMatchEffVsPtPosEtaTOFbc",kRed+1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtNegEtaSPDanyTOFbc=ComputeMatchEff(hPtEtaNegTPCselSPDanyTOFbc,hPtEtaNegTPCselTOFbc,"hMatchEffVsPtNegEtaSPDAnyTOFbc",kGreen+2,33,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaSPDanyTOFbc=ComputeMatchEff(hPtEtaPosTPCselSPDanyTOFbc,hPtEtaPosTPCselTOFbc,"hMatchEffVsPtPosEtaSPDAnyTOFbc",kGreen+2,33,"p_{T} (GeV/c)");

  TH1D* hMatchEffVsPtNegEtaTPCpt=ComputeMatchEff(hPtEtaNegTPCselITSrefTPCpt,hPtEtaNegTPCselTPCpt,"hMatchEffVsPtNegEtaTPCpt",kGray+1,20,"p_{T}^{TPC} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaTPCpt=ComputeMatchEff(hPtEtaPosTPCselITSrefTPCpt,hPtEtaPosTPCselTPCpt,"hMatchEffVsPtPosEtaTPCpt",kGray+1,20,"p_{T}^{TPC} (GeV/c)");
  TH1D* hMatchEffVsPtNegEtaSPDanyTPCpt=ComputeMatchEff(hPtEtaNegTPCselSPDanyTPCpt,hPtEtaNegTPCselTPCpt,"hMatchEffVsPtNegEtaSPDAnyTPCpt",kBlue,33,"p_{T}^{TPC} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaSPDanyTPCpt=ComputeMatchEff(hPtEtaPosTPCselSPDanyTPCpt,hPtEtaPosTPCselTPCpt,"hMatchEffVsPtPosEtaSPDAnyTPCpt",kBlue,33,"p_{T}^{TPC} (GeV/c)");

  TH1D* hMatchEffVsPtNegEtaTPCptTOFbc=ComputeMatchEff(hPtEtaNegTPCselITSrefTPCptTOFbc,hPtEtaNegTPCselTPCptTOFbc,"hMatchEffVsPtNegEtaTPCptTOFbc",kRed-7,20,"p_{T}^{TPC} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaTPCptTOFbc=ComputeMatchEff(hPtEtaPosTPCselITSrefTPCptTOFbc,hPtEtaPosTPCselTPCptTOFbc,"hMatchEffVsPtPosEtaTPCptTOFbc",kRed-7,20,"p_{T}^{TPC} (GeV/c)");
  TH1D* hMatchEffVsPtNegEtaSPDanyTPCptTOFbc=ComputeMatchEff(hPtEtaNegTPCselSPDanyTPCptTOFbc,hPtEtaNegTPCselTPCptTOFbc,"hMatchEffVsPtNegEtaSPDAnyTPCptTOFbc",kSpring-7,33,"p_{T}^{TPC} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaSPDanyTPCptTOFbc=ComputeMatchEff(hPtEtaPosTPCselSPDanyTPCptTOFbc,hPtEtaPosTPCselTPCptTOFbc,"hMatchEffVsPtPosEtaSPDAnyTPCptTOFbc",kSpring-7,33,"p_{T}^{TPC} (GeV/c)");

  hMatchEffVsPtNegEta->SetTitle("#eta<0");
  hMatchEffVsPtPosEta->SetTitle("#eta>0");
  hMatchEffVsPtNegEtaSPDany->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaSPDany->SetTitle("#eta>0");
  hMatchEffVsPtNegEtaTOFbc->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaTOFbc->SetTitle("#eta>0");
  hMatchEffVsPtNegEtaSPDanyTOFbc->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaSPDanyTOFbc->SetTitle("#eta>0");

  hMatchEffVsPtNegEtaTPCpt->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaTPCpt->SetTitle("#eta>0");
  hMatchEffVsPtNegEtaSPDanyTPCpt->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaSPDanyTPCpt->SetTitle("#eta>0");
  hMatchEffVsPtNegEtaTPCptTOFbc->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaTPCptTOFbc->SetTitle("#eta>0");
  hMatchEffVsPtNegEtaSPDanyTPCptTOFbc->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaSPDanyTPCptTOFbc->SetTitle("#eta>0");


  if(maxPtMEPlots<100){
    hMatchEffVsPtNegEta->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtPosEta->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtNegEtaSPDany->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtPosEtaSPDany->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtNegEtaTOFbc->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtPosEtaTOFbc->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtNegEtaSPDanyTOFbc->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtPosEtaSPDanyTOFbc->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtNegEtaTPCpt->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtPosEtaTPCpt->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtNegEtaSPDanyTPCpt->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtPosEtaSPDanyTPCpt->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtNegEtaTPCptTOFbc->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtPosEtaTPCptTOFbc->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtNegEtaSPDanyTPCptTOFbc->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
    hMatchEffVsPtPosEtaSPDanyTPCptTOFbc->GetXaxis()->SetRangeUser(0,maxPtMEPlots);
  }
  for(Int_t ipt=0; ipt<3; ipt++){
    Int_t thePtBin=hMatchEffVsPtNegEta->GetXaxis()->FindBin(ptForTrend[ipt]*0.9999);
    vecMatchEff[ipt]=hMatchEffVsPtPosEta->GetBinContent(thePtBin);
    vecMatchEff[ipt+3]=hMatchEffVsPtNegEta->GetBinContent(thePtBin);
    vecMatchEff[ipt+6]=hMatchEffVsPtPosEtaSPDany->GetBinContent(thePtBin);
    vecMatchEff[ipt+9]=hMatchEffVsPtNegEtaSPDany->GetBinContent(thePtBin);
    vecMatchEff[ipt+12]=hMatchEffVsPtPosEtaTOFbc->GetBinContent(thePtBin);
    vecMatchEff[ipt+15]=hMatchEffVsPtNegEtaTOFbc->GetBinContent(thePtBin);
    vecMatchEff[ipt+18]=hMatchEffVsPtPosEtaSPDanyTOFbc->GetBinContent(thePtBin);
    vecMatchEff[ipt+21]=hMatchEffVsPtNegEtaSPDanyTOFbc->GetBinContent(thePtBin);
    vecErrMatchEff[ipt]=hMatchEffVsPtPosEta->GetBinError(thePtBin);
    vecErrMatchEff[ipt+3]=hMatchEffVsPtNegEta->GetBinError(thePtBin);
    vecErrMatchEff[ipt+6]=hMatchEffVsPtPosEtaSPDany->GetBinError(thePtBin);
    vecErrMatchEff[ipt+9]=hMatchEffVsPtNegEtaSPDany->GetBinError(thePtBin);
    vecErrMatchEff[ipt+12]=hMatchEffVsPtPosEtaTOFbc->GetBinError(thePtBin);
    vecErrMatchEff[ipt+15]=hMatchEffVsPtNegEtaTOFbc->GetBinError(thePtBin);
    vecErrMatchEff[ipt+18]=hMatchEffVsPtPosEtaSPDanyTOFbc->GetBinError(thePtBin);
    vecErrMatchEff[ipt+21]=hMatchEffVsPtNegEtaSPDanyTOFbc->GetBinError(thePtBin);
  }

  TCanvas* cme=new TCanvas("cme","MatchEff All",900,900);
  cme->Divide(2,2);
  cme->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtNegEta->Draw("PE");
  hMatchEffVsPtNegEtaSPDany->Draw("samepe");
  TLegend* legn=new TLegend(0.27,0.17,0.6,0.39);
  legn->AddEntry(hMatchEffVsPtNegEta,"ITSrefit","P");
  legn->AddEntry(hMatchEffVsPtNegEtaSPDany,"SPD any","P");
  legn->Draw();
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
  plotFileName=Form("MatchingEfficiency-AllCharged.%s",outputForm.Data());
  cme->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  TCanvas* cmet=new TCanvas("cmet","MatchEff TPCpt",900,900);
  cmet->Divide(2,2);
  cmet->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtNegEtaTPCpt->Draw("PE");
  hMatchEffVsPtNegEtaSPDanyTPCpt->Draw("samepe");
  TLegend* legnt=new TLegend(0.27,0.17,0.6,0.39);
  legnt->AddEntry(hMatchEffVsPtNegEtaTPCpt,"ITSrefit","P");
  legnt->AddEntry(hMatchEffVsPtNegEtaSPDanyTPCpt,"SPD any","P");
  legnt->Draw();
  cmet->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtPosEtaTPCpt->Draw("PE");
  hMatchEffVsPtPosEtaSPDanyTPCpt->Draw("samepe");
  cmet->cd(3);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtNegEtaTPCptTOFbc->Draw("PE");
  hMatchEffVsPtNegEtaSPDanyTPCptTOFbc->Draw("samepe");
  TLegend* legtt=new TLegend(0.27,0.17,0.89,0.39);
  legtt->AddEntry(hMatchEffVsPtNegEtaTPCptTOFbc,"ITSrefit, TOF bc=0","P");
  legtt->AddEntry(hMatchEffVsPtNegEtaSPDanyTPCptTOFbc,"SPD any, TOF bc=0","P");
  legtt->Draw();
  cmet->cd(4);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtPosEtaTPCptTOFbc->Draw("PE");
  hMatchEffVsPtPosEtaSPDanyTPCptTOFbc->Draw("samepe");
  plotFileName=Form("MatchingEfficiency-TPCpt-AllCharged.%s",outputForm.Data());
  cmet->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

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

  hMatchEffVsPhiNegEtaLowPt->SetTitle("#eta<0 - 0.4<p_{T}<0.7 GeV/c");
  hMatchEffVsPhiPosEtaLowPt->SetTitle("#eta>0 - 0.4<p_{T}<0.7 GeV/c");
  hMatchEffVsPhiNegEtaHighPt->SetTitle("#eta<0 - 1<p_{T}<10 GeV/c");
  hMatchEffVsPhiPosEtaHighPt->SetTitle("#eta>0 - 1<p_{T}<10 GeV/c");

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
  TLegend* legt2=new TLegend(0.17,0.14,0.89,0.29);
  legt2->SetNColumns(2);
  legt2->SetMargin(0.1);
  legt2->AddEntry(hMatchEffVsPhiNegEtaLowPt,"ITSrefit","P");
  legt2->AddEntry(hMatchEffVsPhiNegEtaSPDanyLowPt,"SPD any","P");
  legt2->AddEntry(hMatchEffVsPhiNegEtaLowPtTOFbc,"ITSrefit, TOF bc=0","P");
  legt2->AddEntry(hMatchEffVsPhiNegEtaSPDanyLowPtTOFbc,"SPD any, TOF bc=0","P");
  legt2->Draw();
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

    for(Int_t ipt=0; ipt<3; ipt++){
      Int_t thePtBin=hRatioPosNegEtaPosTPCsel->GetXaxis()->FindBin(ptForTrend[ipt]*0.9999);
      vecPosNeg[ipt]=hRatioPosNegEtaPosTPCsel->GetBinContent(thePtBin);
      vecPosNeg[ipt+3]=hRatioPosNegEtaNegTPCsel->GetBinContent(thePtBin);
      vecPosNeg[ipt+6]=hRatioPosNegEtaPosTPCselSPDany->GetBinContent(thePtBin);
      vecPosNeg[ipt+9]=hRatioPosNegEtaNegTPCselSPDany->GetBinContent(thePtBin);
      vecErrPosNeg[ipt]=hRatioPosNegEtaPosTPCsel->GetBinError(thePtBin);
      vecErrPosNeg[ipt+3]=hRatioPosNegEtaNegTPCsel->GetBinError(thePtBin);
      vecErrPosNeg[ipt+6]=hRatioPosNegEtaPosTPCselSPDany->GetBinError(thePtBin);
      vecErrPosNeg[ipt+9]=hRatioPosNegEtaNegTPCselSPDany->GetBinError(thePtBin);
    }
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
  plotFileName=Form("TracksPhiDistrib-PtBins.%s",outputForm.Data());
  cdist2->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());


  
  const Int_t checkSpecies=2;
  Float_t vecFracBadHyp[9*checkSpecies];
  Float_t vecErrFracBadHyp[9*checkSpecies];

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
      hPtEtaNegGoodHyp[jTy][iSp]->SetTitle(Form("%s - p_{T} - #eta<0",pNames[iSp].Data()));
      hPtEtaPosGoodHyp[jTy][iSp]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hPtEtaPosGoodHyp[jTy][iSp]->SetTitle(Form("%s - p_{T} - #eta>0",pNames[iSp].Data()));
      hPhiEtaNegGoodHyp[jTy][iSp]->GetXaxis()->SetTitle("#varphi (rad)");
      hPhiEtaNegGoodHyp[jTy][iSp]->SetTitle(Form("%s - #varphi - #eta<0",pNames[iSp].Data()));
      hPhiEtaPosGoodHyp[jTy][iSp]->GetXaxis()->SetTitle("#varphi (rad)");
      hPhiEtaPosGoodHyp[jTy][iSp]->SetTitle(Form("%s - #varphi - #eta>0",pNames[iSp].Data()));

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

      hMatchEffGoodVsPtEtaNeg[jTy-1][iSp]->SetTitle(Form("%s -  #eta<0",pNames[iSp].Data()));
      hMatchEffGoodVsPtEtaPos[jTy-1][iSp]->SetTitle(Form("%s - #eta>0",pNames[iSp].Data()));
      hMatchEffGoodVsPhiEtaNeg[jTy-1][iSp]->SetTitle(Form("%s - #eta<0",pNames[iSp].Data()));
      hMatchEffGoodVsPhiEtaPos[jTy-1][iSp]->SetTitle(Form("%s - #eta>0",pNames[iSp].Data()));
      hMatchEffBadVsPtEtaNeg[jTy-1][iSp]->SetTitle(Form("%s - #eta<0",pNames[iSp].Data()));
      hMatchEffBadVsPtEtaPos[jTy-1][iSp]->SetTitle(Form("%s - #eta>0",pNames[iSp].Data()));
      hMatchEffBadVsPhiEtaNeg[jTy-1][iSp]->SetTitle(Form("%s - #eta<0",pNames[iSp].Data()));
      hMatchEffBadVsPhiEtaPos[jTy-1][iSp]->SetTitle(Form("%s - #eta>0",pNames[iSp].Data()));
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
    hRatioBadGoodVsPtEtaNeg[0][iSp]->GetXaxis()->SetRangeUser(0.,maxPtHypoPlots);
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
    hRatioBadGoodVsPtEtaPos[0][iSp]->GetXaxis()->SetRangeUser(0.,maxPtHypoPlots);
    hRatioBadGoodVsPtEtaPos[0][iSp]->Draw();
    hRatioBadGoodVsPtEtaPos[1][iSp]->Draw("same");
    hRatioBadGoodVsPtEtaPos[2][iSp]->Draw("same");
    plotFileName=Form("MassHypoInTracking-%s.%s",pNames[iSp].Data(),outputForm.Data());
    cdistp->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());


    for(Int_t ipt=0; ipt<3; ipt++){
      TString bName1=Form("frac%sBadHypPt%d",pNames[iSp].Data(),TMath::Nint(ptForBadHyp[ipt]*1000.));
      Int_t index1=iSp*9+ipt*3;
      trtree->Branch(bName1.Data(),&vecFracBadHyp[index1],Form("%s/F",bName1.Data()));
      TString errbName1=Form("err%s",bName1.Data());
      trtree->Branch(errbName1.Data(),&vecErrFracBadHyp[index1],Form("%s/F",errbName1.Data()));
      TString bName2=Form("frac%sBadHypPt%dITSref",pNames[iSp].Data(),TMath::Nint(ptForBadHyp[ipt]*1000.));
      Int_t index2=index1+1;
      trtree->Branch(bName2.Data(),&vecFracBadHyp[index2],Form("%s/F",bName2.Data()));
      TString errbName2=Form("err%s",bName2.Data());
      trtree->Branch(errbName2.Data(),&vecErrFracBadHyp[index2],Form("%s/F",errbName2.Data()));
      TString bName3=Form("frac%sBadHypPt%dSPDany",pNames[iSp].Data(),TMath::Nint(ptForBadHyp[ipt]*1000.));
      Int_t index3=index2+1;
      trtree->Branch(bName3.Data(),&vecFracBadHyp[index3],Form("%s/F",bName3.Data()));
      TString errbName3=Form("err%s",bName3.Data());
      trtree->Branch(errbName3.Data(),&vecErrFracBadHyp[index3],Form("%s/F",errbName3.Data()));
      Int_t thePtBin=hRatioBadGoodVsPtEtaPos[0][iSp]->GetXaxis()->FindBin(ptForBadHyp[ipt]*0.9999);
      vecFracBadHyp[index1]=hRatioBadGoodVsPtEtaPos[0][iSp]->GetBinContent(thePtBin);
      vecFracBadHyp[index2]=hRatioBadGoodVsPtEtaPos[1][iSp]->GetBinContent(thePtBin);
      vecFracBadHyp[index3]=hRatioBadGoodVsPtEtaPos[2][iSp]->GetBinContent(thePtBin);
      vecErrFracBadHyp[index1]=hRatioBadGoodVsPtEtaPos[0][iSp]->GetBinError(thePtBin);
      vecErrFracBadHyp[index2]=hRatioBadGoodVsPtEtaPos[1][iSp]->GetBinError(thePtBin);
      vecErrFracBadHyp[index3]=hRatioBadGoodVsPtEtaPos[2][iSp]->GetBinError(thePtBin);
    }

    TCanvas* cmep=new TCanvas(Form("cmep%d",iSp),Form("MatchEff %s",pNames[iSp].Data()),900,900);
    cmep->Divide(2,2);
    cmep->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetTickx();
    gPad->SetTicky();
    if(maxPtMEPlots<100) hMatchEffGoodVsPtEtaNeg[0][iSp]->GetXaxis()->SetRangeUser(0.,maxPtMEPlots);
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
    if(maxPtMEPlots<100) hMatchEffGoodVsPtEtaPos[0][iSp]->GetXaxis()->SetRangeUser(0.,maxPtMEPlots);
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
    plotFileName=Form("MatchingEfficiency-%s.%s",pNames[iSp].Data(),outputForm.Data());
    cmep->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());
  }

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
  
  TH3F* hptresTPC3d=(TH3F*)l->FindObject("hSig1ptCovMatPhiPtTPCsel");
  if(!hptresTPC3d) hptresTPC3d=(TH3F*)l->FindObject("hTPCsig1ptPerClusPhiPtTPCsel"); // old name
  TH3F* hptresITS3d=(TH3F*)l->FindObject("hSig1ptCovMatPhiPtTPCselITSref");
  if(!hptresITS3d) hptresITS3d=(TH3F*)l->FindObject("hTPCsig1ptPerClusPhiPtTPCselITSref");
  TH3F* hptresSPD3d=(TH3F*)l->FindObject("hSig1ptCovMatPhiPtTPCselSPDany");
  if(!hptresSPD3d) hptresSPD3d=(TH3F*)l->FindObject("hTPCsig1ptPerClusPhiPtTPCselSPDany");
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
  plotFileName=Form("PtResolCovMat.%s",outputForm.Data());
  cptcm->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  Float_t vecPtResol[12];
  Float_t vecErrPtResol[12];

  for(Int_t ipt=0; ipt<4; ipt++){
    TString bName1=Form("PtResolTPCCovMatPt%d",TMath::Nint(ptForResol[ipt]*1000.));
    Int_t index1=ipt*3;
    trtree->Branch(bName1.Data(),&vecPtResol[index1],Form("%s/F",bName1.Data()));
    TString errbName1=Form("err%s",bName1.Data());
    trtree->Branch(errbName1.Data(),&vecErrPtResol[index1],Form("%s/F",errbName1.Data()));
    TString bName2=Form("PtResolITSrefCovMatPt%d",TMath::Nint(ptForResol[ipt]*1000.));
    Int_t index2=index1+1;
    trtree->Branch(bName2.Data(),&vecPtResol[index2],Form("%s/F",bName2.Data()));
    TString errbName2=Form("err%s",bName2.Data());
    trtree->Branch(errbName2.Data(),&vecErrPtResol[index2],Form("%s/F",errbName2.Data()));
    TString bName3=Form("PtResolSPDanyCovMatPt%d",TMath::Nint(ptForResol[ipt]*1000.));
    Int_t index3=index2+1;
    trtree->Branch(bName3.Data(),&vecPtResol[index3],Form("%s/F",bName3.Data()));
    TString errbName3=Form("err%s",bName3.Data());
    trtree->Branch(errbName3.Data(),&vecErrPtResol[index3],Form("%s/F",errbName3.Data()));
    Int_t thePtBin=pptresTPC->GetXaxis()->FindBin(ptForResol[ipt]*0.9999);
    vecPtResol[index1]=pptresTPC->GetBinContent(thePtBin);
    vecPtResol[index2]=pptresITS->GetBinContent(thePtBin);
    vecPtResol[index3]=pptresSPD->GetBinContent(thePtBin);
    vecErrPtResol[index1]=pptresTPC->GetBinError(thePtBin);
    vecErrPtResol[index2]=pptresITS->GetBinError(thePtBin);
    vecErrPtResol[index3]=pptresSPD->GetBinError(thePtBin);
  }

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
    isMC=kTRUE;
    FillMeanAndRms(hPtResidVsPtTPCselITSrefPion,gMeanPi,gRmsPi);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefKaon,gMeanK,gRmsK);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefProton,gMeanProt,gRmsProt);
    okRes=kTRUE;
  }
  if(hOneOverPtResidVsPtTPCselITSrefPion && hOneOverPtResidVsPtTPCselITSrefPion->Integral()>0){
    isMC=kTRUE;
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
    plotFileName=Form("PtResolRecoMinusGen.%s",outputForm.Data());
    c2->SaveAs(plotFileName.Data());
    if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());
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
  Float_t mK0=fmassk0->GetParameter(4);
  Float_t emK0=fmassk0->GetParError(4);
  Float_t sigK0=fmassk0->GetParameter(5);
  Float_t esigK0=fmassk0->GetParError(5);
  trtree->Branch("massK0",&mK0,"massK0/F");
  trtree->Branch("errmassK0",&emK0,"errmassK0/F");
  trtree->Branch("sigmaK0",&sigK0,"sigmaK0/F");
  trtree->Branch("errsigmaK0",&esigK0,"errsigmaK0/F");
  cv0->cd(2);
  hInvMassLambda->Draw();
  InitFuncAndFit(hInvMassLambda,fmassL,kFALSE,isMC);
  Float_t mL=fmassL->GetParameter(4);
  Float_t emL=fmassL->GetParError(4);
  Float_t sigL=fmassL->GetParameter(5);
  Float_t esigL=fmassL->GetParError(5);
  trtree->Branch("massLambda",&mL,"massLambda/F");
  trtree->Branch("errmassLambda",&emL,"errmassLambda/F");
  trtree->Branch("sigmaLambda",&sigL,"sigmaLambda/F");
  trtree->Branch("errsigmaLambda",&esigL,"errsigmaLambda/F");
  cv0->cd(3);
  hInvMassAntiLambda->Draw();
  InitFuncAndFit(hInvMassAntiLambda,fmassL,kFALSE,isMC);
  Float_t mLb=fmassL->GetParameter(4);
  Float_t emLb=fmassL->GetParError(4);
  Float_t sigLb=fmassL->GetParameter(5);
  Float_t esigLb=fmassL->GetParError(5);
  trtree->Branch("massLambdabar",&mLb,"massLambdabar/F");
  trtree->Branch("errmassLambdabar",&emLb,"errmassLambdabar/F");
  trtree->Branch("sigmaLambdabar",&sigLb,"sigmaLambdabar/F");
  trtree->Branch("errsigmaLambdabar",&esigLb,"errsigmaLambdabar/F");
  plotFileName=Form("MassSpectraV0-integrated.%s",outputForm.Data());
  cv0->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  TCanvas* clam=new TCanvas("clam","Lambda vs R",1400,900);
  clam->Divide(2,2);
  clam->cd(1);
  hInvMassLambdaR1->Draw();
  InitFuncAndFit(hInvMassLambdaR1,fmassL,kFALSE,isMC);
  Float_t mLRad1=fmassL->GetParameter(4);
  Float_t emLRad1=fmassL->GetParError(4);
  Float_t sigLRad1=fmassL->GetParameter(5);
  Float_t esigLRad1=fmassL->GetParError(5);
  TString radbininfo=Form("_%dRad%d",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(1)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z1+1)));
  trtree->Branch(Form("massLambda%s",radbininfo.Data()),&mLRad1,Form("massLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("errmassLambda%s",radbininfo.Data()),&emLRad1,Form("errmassLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("sigmaLambda%s",radbininfo.Data()),&sigLRad1,Form("sigmaLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("errsigmaLambda%s",radbininfo.Data()),&esigLRad1,Form("errsigmaLambda%s/F",radbininfo.Data()));
  TLatex* tr1=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(1)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z1+1))));
  tr1->SetNDC();
  tr1->SetTextFont(43);
  tr1->SetTextSize(26);
  tr1->Draw();
  clam->cd(2);
  hInvMassLambdaR2->Draw();
  InitFuncAndFit(hInvMassLambdaR2,fmassL,kFALSE,isMC);
  Float_t mLRad2=fmassL->GetParameter(4);
  Float_t emLRad2=fmassL->GetParError(4);
  Float_t sigLRad2=fmassL->GetParameter(5);
  Float_t esigLRad2=fmassL->GetParError(5);
  radbininfo=Form("_%dRad%d",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z1+1)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z2+1)));
  trtree->Branch(Form("massLambda%s",radbininfo.Data()),&mLRad2,Form("massLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("errmassLambda%s",radbininfo.Data()),&emLRad2,Form("errmassLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("sigmaLambda%s",radbininfo.Data()),&sigLRad2,Form("sigmaLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("errsigmaLambda%s",radbininfo.Data()),&esigLRad2,Form("errsigmaLambda%s/F",radbininfo.Data()));
  TLatex* tr2=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z1+1)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z2+1))));
  tr2->SetNDC();
  tr2->SetTextFont(43);
  tr2->SetTextSize(26);
  tr2->Draw();
  clam->cd(3);
  hInvMassLambdaR3->Draw();
  InitFuncAndFit(hInvMassLambdaR3,fmassL,kFALSE,isMC);
  Float_t mLRad3=fmassL->GetParameter(4);
  Float_t emLRad3=fmassL->GetParError(4);
  Float_t sigLRad3=fmassL->GetParameter(5);
  Float_t esigLRad3=fmassL->GetParError(5);
  radbininfo=Form("_%dRad%d",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z3)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z4+1)));
  trtree->Branch(Form("massLambda%s",radbininfo.Data()),&mLRad3,Form("massLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("errmassLambda%s",radbininfo.Data()),&emLRad3,Form("errmassLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("sigmaLambda%s",radbininfo.Data()),&sigLRad3,Form("sigmaLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("errsigmaLambda%s",radbininfo.Data()),&esigLRad3,Form("errsigmaLambda%s/F",radbininfo.Data()));
  TLatex* tr3=new TLatex(0.65,0.6,Form("%d<R<%d cm",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z3)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z4+1))));
  tr3->SetNDC();
  tr3->SetTextFont(43);
  tr3->SetTextSize(26);
  tr3->Draw();
  clam->cd(4);
  hInvMassLambdaR4->Draw();
  InitFuncAndFit(hInvMassLambdaR4,fmassL,kFALSE,isMC);
  Float_t mLRad4=fmassL->GetParameter(4);
  Float_t emLRad4=fmassL->GetParError(4);
  Float_t sigLRad4=fmassL->GetParameter(5);
  Float_t esigLRad4=fmassL->GetParError(5);
  radbininfo=Form("_%dRad%d",TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z5)),TMath::Nint(hInvMassLambda3d->GetZaxis()->GetBinLowEdge(z6+1)));
  trtree->Branch(Form("massLambda%s",radbininfo.Data()),&mLRad4,Form("massLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("errmassLambda%s",radbininfo.Data()),&emLRad4,Form("errmassLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("sigmaLambda%s",radbininfo.Data()),&sigLRad4,Form("sigmaLambda%s/F",radbininfo.Data()));
  trtree->Branch(Form("errsigmaLambda%s",radbininfo.Data()),&esigLRad4,Form("errsigmaLambda%s/F",radbininfo.Data()));
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
  Float_t mK0Pt1=fmassk0->GetParameter(4);
  Float_t emK0Pt1=fmassk0->GetParError(4);
  Float_t sigK0Pt1=fmassk0->GetParameter(5);
  Float_t esigK0Pt1=fmassk0->GetParError(5);
  TString ptbininfo=Form("_%dPt%d",TMath::Nint(hInvMassK0s3d->GetYaxis()->GetBinLowEdge(1)*1000.),TMath::Nint(hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p1+1)*1000.));
  trtree->Branch(Form("massK0%s",ptbininfo.Data()),&mK0Pt1,Form("massK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("errmassK0%s",ptbininfo.Data()),&emK0Pt1,Form("errmassK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("sigmaK0%s",ptbininfo.Data()),&sigK0Pt1,Form("sigmaK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("errsigmaK0%s",ptbininfo.Data()),&esigK0Pt1,Form("errsigmaK0%s/F",ptbininfo.Data()));
  TLatex* tp1=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(1),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p1+1)));
  tp1->SetNDC();
  tp1->SetTextFont(43);
  tp1->SetTextSize(26);
  tp1->Draw();
  ck0->cd(2);
  hInvMassK0sP2->Draw();
  InitFuncAndFit(hInvMassK0sP2,fmassk0,kTRUE,isMC);
  Float_t mK0Pt2=fmassk0->GetParameter(4);
  Float_t emK0Pt2=fmassk0->GetParError(4);
  Float_t sigK0Pt2=fmassk0->GetParameter(5);
  Float_t esigK0Pt2=fmassk0->GetParError(5);
  ptbininfo=Form("_%dPt%d",TMath::Nint(hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p1+1)*1000.),TMath::Nint(hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p2+1)*1000.));
  trtree->Branch(Form("massK0%s",ptbininfo.Data()),&mK0Pt2,Form("massK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("errmassK0%s",ptbininfo.Data()),&emK0Pt2,Form("errmassK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("sigmaK0%s",ptbininfo.Data()),&sigK0Pt2,Form("sigmaK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("errsigmaK0%s",ptbininfo.Data()),&esigK0Pt2,Form("errsigmaK0%s/F",ptbininfo.Data()));
  TLatex* tp2=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p1+1),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p2+1)));
  tp2->SetNDC();
  tp2->SetTextFont(43);
  tp2->SetTextSize(26);
  tp2->Draw();
  ck0->cd(3);
  hInvMassK0sP3->Draw();
  InitFuncAndFit(hInvMassK0sP3,fmassk0,kTRUE,isMC);
  Float_t mK0Pt3=fmassk0->GetParameter(4);
  Float_t emK0Pt3=fmassk0->GetParError(4);
  Float_t sigK0Pt3=fmassk0->GetParameter(5);
  Float_t esigK0Pt3=fmassk0->GetParError(5);
  ptbininfo=Form("_%dPt%d",TMath::Nint(hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p2+1)*1000.),TMath::Nint(hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p3+1)*1000.));
  trtree->Branch(Form("massK0%s",ptbininfo.Data()),&mK0Pt3,Form("massK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("errmassK0%s",ptbininfo.Data()),&emK0Pt3,Form("errmassK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("sigmaK0%s",ptbininfo.Data()),&sigK0Pt3,Form("sigmaK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("errsigmaK0%s",ptbininfo.Data()),&esigK0Pt3,Form("errsigmaK0%s/F",ptbininfo.Data()));
  TLatex* tp3=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p2+1),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p3+1)));
  tp3->SetNDC();
  tp3->SetTextFont(43);
  tp3->SetTextSize(26);
  tp3->Draw();
  ck0->cd(4);
  hInvMassK0sP4->Draw();
  InitFuncAndFit(hInvMassK0sP4,fmassk0,kTRUE,isMC);
  Float_t mK0Pt4=fmassk0->GetParameter(4);
  Float_t emK0Pt4=fmassk0->GetParError(4);
  Float_t sigK0Pt4=fmassk0->GetParameter(5);
  Float_t esigK0Pt4=fmassk0->GetParError(5);
  ptbininfo=Form("_%dPt%d",TMath::Nint(hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p4)*1000.),TMath::Nint(hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p5+1)*1000.));
  trtree->Branch(Form("massK0%s",ptbininfo.Data()),&mK0Pt4,Form("massK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("errmassK0%s",ptbininfo.Data()),&emK0Pt4,Form("errmassK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("sigmaK0%s",ptbininfo.Data()),&sigK0Pt4,Form("sigmaK0%s/F",ptbininfo.Data()));
  trtree->Branch(Form("errsigmaK0%s",ptbininfo.Data()),&esigK0Pt4,Form("errsigmaK0%s/F",ptbininfo.Data()));
  TLatex* tp4=new TLatex(0.6,0.6,Form("%.1f<p_{T}<%.1f GeV/c",hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p4),hInvMassK0s3d->GetYaxis()->GetBinLowEdge(p5+1)));
  tp4->SetNDC();
  tp4->SetTextFont(43);
  tp4->SetTextSize(26);
  tp4->Draw();
  plotFileName=Form("K0s-MassSpectra-VsPt.%s",outputForm.Data());
  ck0->SaveAs(plotFileName.Data());
  if(outputForm=="pdf") pdfFileNames+=Form("%s ",plotFileName.Data());

  // K0 pt resolution vs. pt
  const Int_t nPtBinsK0=7;
  Double_t ptbinlimsK0[nPtBinsK0+1]={0.,0.4,0.8,1.2,2.0,3.,4.,5.};
  TH1F* hSigmaK0AllR=new TH1F("hSigmaK0AllR"," ; p_{T} (GeV/c) ; #sigma_{K0} (MeV/c^{2})",nPtBinsK0,ptbinlimsK0);
  TH1F* hSigmaK0R4=new TH1F("hSigmaK0R4"," ; p_{T} (GeV/c) ; #sigma_{K0} (MeV/c^{2})",nPtBinsK0,ptbinlimsK0);
  TH1F* hMassK0AllR=new TH1F("hMassK0AllR"," ; p_{T} (GeV/c) ; #mu_{K0} (GeV/c^{2})",nPtBinsK0,ptbinlimsK0);
  TH1F* hYieldK0AllR=new TH1F("hYieldK0AllR"," ; p_{T} (GeV/c) ; N_{K0}/event",nPtBinsK0,ptbinlimsK0);
  TCanvas* ctmpk0=new TCanvas("ctmpk0","K0s vs. pt R<4",1600,900);
  ctmpk0->Divide(4,4);
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
    ctmpk0->cd(ipt+9);
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
    TFile* fouttree=new TFile("trending.root","recreate");
    trtree->Write();
    TDirectory* outdir=fouttree->mkdir(df->GetName());
    outdir->cd();
    l->Write(l->GetName(),1);
    fouttree->Close();
    delete fouttree;
  }

  if(outputForm=="pdf") gSystem->Exec(Form("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=PlotsESDTrackQA.pdf %s",pdfFileNames.Data()));
  
  printf("SUMMARY:\n");
  printf("Number of events used in the plots = %d\n",nSelectedEvents);
}

void FillMeanAndRms(TH2F* h2d, TGraphErrors* gMean, TGraphErrors* gRms){
  Int_t jpt=0;
  Int_t jptr=0;
  for(Int_t j=1; j<=h2d->GetNbinsX(); j++){
    Double_t pt=h2d->GetXaxis()->GetBinCenter(j);
    Double_t ept=0;
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
  if(hnumer->GetSumw2N()==0) hnumer->Sumw2();
  if(hdenom->GetSumw2N()==0) hdenom->Sumw2();
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
  if(hnumer->GetSumw2N()==0) hnumer->Sumw2();
  if(hdenom->GetSumw2N()==0) hdenom->Sumw2();
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
  h1->GetXaxis()->SetRangeUser(0.,maxPtHypoPlots);
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

Double_t fp2bkgk0(Double_t *x, Double_t *par){
  if (x[0] > 0.47 && x[0] < 0.53) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}
