#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TObjArray.h>
#include "AliCDBEntry.h"
#include "AliITSCalibrationSDD.h"
#endif

/*  $Id$    */

// Macro to plot the calibration parameters from the OCDB files 
// created from PEDESTAL and PULSER runs vs. Time
// Origin: F. Prino (prino@to.infn.it)

void PlotCalibSDDVsTime(Int_t year=2010, Int_t firstRun=77677, 
			Int_t lastRun=999999999){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetTitleOffset(1.4,"Y");  


  TGrid::Connect("alien:",0,0,"t");
  gSystem->Exec(Form("gbbox find \"/alice/data/%d/OCDB/ITS/Calib/CalibSDD\" \"Run*.root\" > runCalibAlien.txt",year));
  FILE* listruns=fopen("runCalibAlien.txt","r");

  TH1F* hbase=new TH1F("hbase","",60,0.5,120.5);
  TH1F* hnoise=new TH1F("hnoise","",100,0.,7.);
  TH1F* hgain=new TH1F("hgain","",100,0.,4.);
  TH1F* hchstatus=new TH1F("hchstatus","",2,-0.5,1.5);
  TH1F* hchstatus3=new TH1F("hchstatus3","",2,-0.5,1.5);
  TH1F* hchstatus4=new TH1F("hchstatus4","",2,-0.5,1.5);
  TGraphErrors* gbasevstim=new TGraphErrors(0);
  TGraphErrors* gnoisevstim=new TGraphErrors(0);
  TGraphErrors* ggainvstim=new TGraphErrors(0);
  TGraphErrors* gstatvstim=new TGraphErrors(0);
  TGraphErrors* gfracvstim=new TGraphErrors(0);
  TGraphErrors* gfrac3vstim=new TGraphErrors(0);
  TGraphErrors* gfrac4vstim=new TGraphErrors(0);
  gbasevstim->SetName("gbasevstim");
  gnoisevstim->SetName("gnoisevstim");
  ggainvstim->SetName("ggainvstim");
  gstatvstim->SetName("gstatvstim");
  gfracvstim->SetName("gfracvstim");
  gfrac3vstim->SetName("gfrac3vstim");
  gfrac4vstim->SetName("gfrac4vstim");
  gbasevstim->SetTitle("Baseline vs. run");
  gnoisevstim->SetTitle("Noise vs. run");
  ggainvstim->SetTitle("Gain vs. run");
  gstatvstim->SetTitle("Good Anodes vs. run");
  gfracvstim->SetTitle("Fraction of Good Anodes vs. run");
  gfrac3vstim->SetTitle("Fraction of Good Anodes vs. run");
  gfrac4vstim->SetTitle("Fraction of Good Anodes vs. run");


  Char_t filnam[200],filnamalien[200];
  Int_t iPoint=0;
  Int_t nrun,nrun2,nv,ns;

  while(!feof(listruns)){
    hbase->Reset();
    hnoise->Reset();
    hgain->Reset();
    hchstatus->Reset();
    hchstatus3->Reset();
    hchstatus4->Reset();
    fscanf(listruns,"%s\n",filnam);    
    Char_t directory[100];
    sprintf(directory,"/alice/data/%d",year);
    if(!strstr(filnam,directory)) continue;
    sscanf(filnam,"/alice/data/%d/OCDB/ITS/Calib/CalibSDD/Run%d_%d_v%d_s%d.root",&year,&nrun,&nrun2,&nv,&ns);
    if(year==2009 && (nrun<85639 && nrun2> 85639)) continue; // protection for files with swapped ladders 4-5 of layer 3 
    if(year==2009 && (nrun>100000 && nv< 184)) continue; // protection for files with swapped ladder 0-1 of layer 4
    if(year==2010 && (nrun>=114603 && nv< 98)) continue; // protection for files without treatment of masked hybrids 
    if(nrun<firstRun) continue;
    if(nrun>lastRun) continue;
    sprintf(filnamalien,"alien://%s",filnam);
    printf("Open file: %s\n",filnam);
    TFile *f= TFile::Open(filnamalien);  
    AliCDBEntry *ent=(AliCDBEntry*)f->Get("AliCDBEntry");
    TObjArray *calSDD = (TObjArray *)ent->GetObject();
    printf("Run %d Entries in array=%d \n",nrun,calSDD->GetEntriesFast());


    AliITSCalibrationSDD *cal;
    for(Int_t iMod=0; iMod<260;iMod++){
      cal=(AliITSCalibrationSDD*)calSDD->At(iMod);
      if(cal==0) continue;
      for(Int_t iAn=0; iAn<512; iAn++){
	Int_t ic=cal->GetChip(iAn);
	Float_t base=cal->GetBaseline(iAn);
	Float_t noise=cal->GetNoiseAfterElectronics(iAn);
	Float_t gain=cal->GetChannelGain(iAn);
	if(cal->IsBadChannel(iAn)){
	  hchstatus->Fill(0);
	  if(iMod<84) hchstatus3->Fill(0);
	  else hchstatus4->Fill(0);
	}
	if(!cal->IsBadChannel(iAn) && !cal->IsChipBad(ic) && !cal->IsBad() ){
	  hbase->Fill(base);
	  hchstatus->Fill(1);
	  if(iMod<84) hchstatus3->Fill(1);
	  else hchstatus4->Fill(1);
	  hnoise->Fill(noise);
	  hgain->Fill(gain);
	}
      } 
    }
    printf("Run %d <Base> = %f <Noise> =%f Entries = %d\n",nrun,hbase->GetMean(),hnoise->GetMean(),(Int_t)hbase->GetEntries());
    if((Int_t)hbase->GetEntries()==0) continue;
    gbasevstim->SetPoint(iPoint,(Double_t)nrun,hbase->GetMean());
    gbasevstim->SetPointError(iPoint,0.,hbase->GetRMS());
    gnoisevstim->SetPoint(iPoint,(Double_t)nrun,hnoise->GetMean());
    gnoisevstim->SetPointError(iPoint,0.,hnoise->GetRMS());
    ggainvstim->SetPoint(iPoint,(Double_t)nrun,hgain->GetMean());
    ggainvstim->SetPointError(iPoint,0.,hgain->GetRMS());
    gstatvstim->SetPoint(iPoint,(Double_t)nrun,hchstatus->GetBinContent(2));
    gfracvstim->SetPoint(iPoint,(Double_t)nrun,hchstatus->GetBinContent(2)/260./512.);
    gfrac3vstim->SetPoint(iPoint,(Double_t)nrun,hchstatus3->GetBinContent(2)/84./512.);
    gfrac4vstim->SetPoint(iPoint,(Double_t)nrun,hchstatus4->GetBinContent(2)/176./512.);
    iPoint++;
    f->Close();
  }

  TFile *ofil=new TFile(Form("Calib%dVsTime.root",year),"recreate");
  gbasevstim->Write();
  gnoisevstim->Write();
  ggainvstim->Write();
  gstatvstim->Write();
  ofil->Close();

  TCanvas* cbase=new TCanvas("cbase","Baselines");
  gbasevstim->SetFillColor(kOrange-2);
  gbasevstim->SetMarkerStyle(20);
  gbasevstim->Draw("AP3");
  gbasevstim->Draw("PLXSAME");
  gbasevstim->SetMinimum(0.);
  gbasevstim->SetMaximum(70.);  
  gbasevstim->GetXaxis()->SetTitle("Run number");
  gbasevstim->GetYaxis()->SetTitle("<Baseline> (ADC counts)");
  cbase->SaveAs(Form("BaseRun%d.gif",year));

  TCanvas* cnoise=new TCanvas("cnoise","Noise");
  gnoisevstim->SetFillColor(kOrange-2);
  gnoisevstim->SetMarkerStyle(20);
  gnoisevstim->Draw("AP3");
  gnoisevstim->Draw("PLXSAME");
  gnoisevstim->SetMinimum(0.);
  gnoisevstim->SetMaximum(4.);
  gnoisevstim->GetXaxis()->SetTitle("Run number");
  gnoisevstim->GetYaxis()->SetTitle("<Noise> (ADC counts)");
  cnoise->SaveAs(Form("NoiseRun%d.gif",year));

  TCanvas* cgain=new TCanvas("cgain","Gain");
  ggainvstim->SetFillColor(kOrange-2);
  ggainvstim->SetMarkerStyle(20);
  ggainvstim->Draw("AP3");
  ggainvstim->Draw("PLXSAME");
  ggainvstim->SetMinimum(0.);
  ggainvstim->SetMaximum(4.);
  ggainvstim->GetXaxis()->SetTitle("Run number");
  ggainvstim->GetYaxis()->SetTitle("<Gain> (ADC/DAC)");
  cgain->SaveAs(Form("GainRun%d.gif",year));

  TCanvas* cstatus=new TCanvas("cstatus","Good channels");
  gstatvstim->SetFillColor(kOrange-2);
  gstatvstim->SetMarkerStyle(20);
  gstatvstim->Draw("AP3");
  gstatvstim->Draw("PLXSAME");
  gstatvstim->SetMinimum(100000.);
  gstatvstim->SetMaximum(133000.);
  gstatvstim->GetXaxis()->SetTitle("Run number");
  gstatvstim->GetYaxis()->SetTitle("Number of good anodes in acquisition");
  cstatus->SaveAs(Form("GoodAnodesRun%d.gif",year));

  TCanvas* cfrac=new TCanvas("cfrac","Fraction of Good");
  gfracvstim->SetMarkerStyle(20);
  gfrac3vstim->SetMarkerStyle(22);
  gfrac3vstim->SetMarkerColor(2);
  gfrac3vstim->SetLineColor(2);
  gfrac4vstim->SetMarkerStyle(23);
  gfrac4vstim->SetMarkerColor(4);
  gfrac4vstim->SetLineColor(4);
  gfracvstim->Draw("APL");
  gfrac3vstim->Draw("PLSAME");
  gfrac4vstim->Draw("PLSAME");
  gfracvstim->SetMinimum(0.7);
  gfracvstim->SetMaximum(1.05);
  gfracvstim->GetXaxis()->SetTitle("Run number");
  gfracvstim->GetYaxis()->SetTitle("Fraction of good anodes in acquisition");
  TLegend* leg=new TLegend(0.2,0.15,0.45,0.35);
  leg->SetFillColor(0);
  TLegendEntry* entr=leg->AddEntry(gfrac3vstim,"Layer 3","P");
  entr->SetTextColor(2);
  entr=leg->AddEntry(gfrac4vstim,"Layer 4","P");
  entr->SetTextColor(4);
  entr=leg->AddEntry(gfracvstim,"All","P");
  entr->SetTextColor(1);
  leg->Draw();
  cfrac->SaveAs(Form("FractionGoodRun%d.gif",year));

}
