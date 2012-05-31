#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TGraph.h>
#include <TExec.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TPaveStats.h>
#include "AliCDBEntry.h"
#include "AliITSresponseSDD.h"
#endif

/*  $Id: ShowResponseSDD.C 44072 2010-10-04 15:48:32Z prino $    */

// Macro to plot the calibration parameters 
// (time zero, ADCtokeV and Vdrift correction)
// from the OCDB file created from SDD offline calibration 
// (OCDB/ITS/Calib/RespSDD)

void ShowResponseSDD(TString filename="$ALICE_ROOT/OCDB/ITS/Calib/RespSDD/Run0_999999999_v0_s0.root"){
  if(filename.Contains("alien")){
    TGrid::Connect("alien:");
  }
  TFile* fr=TFile::Open(filename.Data());
  AliCDBEntry* e=(AliCDBEntry*)fr->Get("AliCDBEntry");
  e->PrintMetaData();
  AliITSresponseSDD* r=(AliITSresponseSDD*)e->GetObject();
  TH1F* hTimeZero=new TH1F("hTimeZero","",260,239.5,499.5);
  TH1F* hVdriftCorrLeft=new TH1F("hVdriftCorrLeft","",260,239.5,499.5);
  TH1F* hVdriftCorrRight=new TH1F("hVdriftCorrRight","",260,239.5,499.5);
  TH1F* hdistVdriftCorrLeft=new TH1F("hdistVdriftCorrLeft","",100,-0.2,0.2);
  TH1F* hdistVdriftCorrRight=new TH1F("hdistVdriftCorrRight","",100,-0.2,0.2);
  TH1F* hADCtokeV=new TH1F("hADCtokeV","",260,239.5,499.5);
  TH1F* hADCvsTime=new TH1F("hADCvsTime","",260,239.5,499.5);
  Float_t averTz=0.;
  Float_t averCv0=0.;
  Float_t averCv1=0.;
  Float_t averAk=0.;
  Float_t averAt=0.;
  for(Int_t iMod=240; iMod<500; iMod++){
    Float_t tz=r->GetTimeZero(iMod);
    Float_t cv0=r->GetDeltaVDrift(iMod,kFALSE);
    Float_t cv1=r->GetDeltaVDrift(iMod,kTRUE);
    Float_t ak=r->GetADCtokeV(iMod);
    Float_t at=r->GetADCvsDriftTime(iMod);
    hTimeZero->SetBinContent(iMod-240+1,tz);    
    hVdriftCorrLeft->SetBinContent(iMod-240+1,cv0);
    hVdriftCorrRight->SetBinContent(iMod-240+1,cv1);
    hdistVdriftCorrLeft->Fill(cv0);
    hdistVdriftCorrRight->Fill(cv1);
    hADCtokeV->SetBinContent(iMod-240+1,ak);
    hADCvsTime->SetBinContent(iMod-240+1,at);
    averTz+=tz;
    averCv0+=cv0;
    averCv1+=cv1;
    averAk+=ak;
    averAt+=at;
  }
  averTz/=260.;
  averCv0/=260.;
  averCv1/=260.;
  averAk/=260.;
  averAt/=260.;

  hTimeZero->SetMarkerStyle(20);
  hADCtokeV->SetMarkerStyle(20);
  hADCvsTime->SetMarkerStyle(20);
  hVdriftCorrLeft->SetMarkerStyle(22);
  hVdriftCorrLeft->SetMarkerColor(2);
  hVdriftCorrRight->SetMarkerStyle(24);
  hVdriftCorrRight->SetMarkerColor(4);
  hTimeZero->SetMaximum(hTimeZero->GetMaximum()*1.2);
  hADCtokeV->SetMaximum(hADCtokeV->GetMaximum()*1.2);
  hADCvsTime->SetMaximum(hADCvsTime->GetMaximum()*1.2);
  Float_t maxV=TMath::Max(hVdriftCorrLeft->GetMaximum(),hVdriftCorrRight->GetMaximum());
  Float_t minV=TMath::Min(hVdriftCorrLeft->GetMinimum(),hVdriftCorrRight->GetMinimum());
  Float_t scale=TMath::Max(TMath::Abs(maxV),TMath::Abs(minV));
  if(scale<0.02){
    hVdriftCorrLeft->SetMinimum(-0.05);
    hVdriftCorrLeft->SetMaximum(0.05);
  }else{
    hVdriftCorrLeft->SetMaximum(1.2*scale);
    hVdriftCorrLeft->SetMinimum(-1.2*scale);
  }

  //  gStyle->SetOptStat(0);

  printf("Charge vs. Time correction factor = %f\n",r->GetChargevsTime());

  TCanvas* c1=new TCanvas("c1","Time Zero");
  hTimeZero->SetStats(0);
  hTimeZero->Draw("P");
  hTimeZero->GetXaxis()->SetTitle("Module Id");
  hTimeZero->GetYaxis()->SetTitle("Time Zero (ns)");
  TLatex* tt=new TLatex(0.2,0.8,Form("Average Tzero = %.2f",averTz));
  tt->SetNDC();
  tt->Draw();
  c1->Modified();

  TCanvas* c2=new TCanvas("c2","Vdrift Corr",800,900);
  c2->Divide(1,2);
  c2->cd(1);
  hVdriftCorrLeft->SetStats(0);
  hVdriftCorrLeft->Draw("P");
  hVdriftCorrRight->Draw("PSAME");
  hVdriftCorrLeft->GetXaxis()->SetTitle("Module Id");
  hVdriftCorrLeft->GetYaxis()->SetTitle("Vdrift correction (#mum/ns)");
  TLatex* tc0=new TLatex(0.2,0.8,Form("Average Vdrift corr Left = %.3f",averCv0));
  tc0->SetNDC();
  tc0->SetTextColor(2);
  tc0->Draw();
  TLatex* tc1=new TLatex(0.2,0.72,Form("Average Vdrift corr Right = %.3f",averCv1));
  tc1->SetNDC();
  tc1->SetTextColor(4);
  tc1->Draw();
  c2->cd(2);
  hdistVdriftCorrLeft->SetLineColor(2);
  hdistVdriftCorrLeft->GetXaxis()->SetTitle("Vdrift correction (#mum/ns)");
  hdistVdriftCorrLeft->Draw(); 
  c2->Update();
  TPaveStats *stL=(TPaveStats*)hdistVdriftCorrLeft->GetListOfFunctions()->FindObject("stats");
  stL->SetY1NDC(0.71);
  stL->SetY2NDC(0.9);
  stL->SetTextColor(2);
  hdistVdriftCorrRight->SetLineColor(4);
  hdistVdriftCorrRight->Draw("SAMES");
  c2->Update();
  TPaveStats *stR=(TPaveStats*)hdistVdriftCorrRight->GetListOfFunctions()->FindObject("stats");
  stR->SetY1NDC(0.51);
  stR->SetY2NDC(0.7);
  stR->SetTextColor(4);

  TCanvas* c3=new TCanvas("c3","ADC calib",800,900);
  c3->Divide(1,2);
  c3->cd(1);
  hADCvsTime->SetStats(0);
  hADCvsTime->Draw("P");
  hADCvsTime->GetXaxis()->SetTitle("Module Id");
  hADCvsTime->GetYaxis()->SetTitle("ADC vs. Drift Time Slope (ADC/ns)");
  TLatex* ts=new TLatex(0.2,0.8,Form("Average ADCvsTime = %.3f",averAt));
  ts->SetNDC();
  ts->Draw(); 
  c3->cd(2);
  hADCtokeV->SetStats(0);
  hADCtokeV->Draw("P");
  hADCtokeV->GetXaxis()->SetTitle("Module Id");
  hADCtokeV->GetYaxis()->SetTitle("ADC to keV conv. factor");
  TLatex* ta=new TLatex(0.2,0.8,Form("Average ADCtokeV = %.3f",averAk));
  ta->SetNDC();
  ta->Draw(); 

}
