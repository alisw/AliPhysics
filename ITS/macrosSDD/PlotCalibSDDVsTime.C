#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
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

/*  $Id: PlotCalibSDDVsTime.C 41568 2010-06-03 09:08:39Z prino $    */

// Macro to plot the calibration parameters from the OCDB files 
// created from PEDESTAL and PULSER runs vs. Time
// Origin: F. Prino (prino@to.infn.it)

void PlotCalibSDDVsTime(Int_t year=2016, Int_t firstRun=259000,
			Int_t lastRun=999999999,
			Int_t selectedMod=-1){

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
  TH1F* hlowzsthr3=new TH1F("hlowzsthr3","",256,-0.5,255.5);
  TH1F* hlowzsthr4=new TH1F("hlowzsthr4","",256,-0.5,255.5);
  TH1F* hhizsthr3=new TH1F("hhizsthr3","",256,-0.5,255.5);
  TH1F* hhizsthr4=new TH1F("hhizsthr4","",256,-0.5,255.5);

  TH2F* hlowzsthr3VsRun=new TH2F("hlowzsthr3VsRun"," ; Run Number ; Low ZS Threshold",400,-0.5,399.5,256,-0.5,255.5);
  TH2F* hlowzsthr4VsRun=new TH2F("hlowzsthr4VsRun"," ; Run Number ; Low ZS Threshold",400,-0.5,399.5,256,-0.5,255.5);
  TH2F* hhizsthr3VsRun=new TH2F("hhizsthr3VsRun"," ; Run Number ; High ZS Threshold",400,-0.5,399.5,256,-0.5,255.5);
  TH2F* hhizsthr4VsRun=new TH2F("hhizsthr4VsRun"," ; Run Number ; High ZS Threshold",400,-0.5,399.5,256,-0.5,255.5);

  TGraphErrors* gbasevstim=new TGraphErrors(0);
  TGraphErrors* gnoisevstim=new TGraphErrors(0);
  TGraphErrors* ggainvstim=new TGraphErrors(0);
  TGraphErrors* gstatvstim=new TGraphErrors(0);
  TGraphErrors* gfracvstim=new TGraphErrors(0);
  TGraphErrors* gfrac3vstim=new TGraphErrors(0);
  TGraphErrors* gfrac4vstim=new TGraphErrors(0);
  TGraphErrors* gzslow3vstim=new TGraphErrors(0);
  TGraphErrors* gzslow4vstim=new TGraphErrors(0);
  TGraphErrors* gzshi3vstim=new TGraphErrors(0);
  TGraphErrors* gzshi4vstim=new TGraphErrors(0);

  gbasevstim->SetName("gbasevstim");
  gnoisevstim->SetName("gnoisevstim");
  ggainvstim->SetName("ggainvstim");
  gstatvstim->SetName("gstatvstim");
  gfracvstim->SetName("gfracvstim");
  gfrac3vstim->SetName("gfrac3vstim");
  gfrac4vstim->SetName("gfrac4vstim");
  gzslow3vstim->SetName("gzslow3vstim");
  gzslow4vstim->SetName("gzslow4vstim");
  gzshi3vstim->SetName("gzshi3vstim");
  gzshi4vstim->SetName("gzshi4vstim");

  gbasevstim->SetTitle("Baseline vs. run");
  gnoisevstim->SetTitle("Noise vs. run");
  ggainvstim->SetTitle("Gain vs. run");
  gstatvstim->SetTitle("Good Anodes vs. run");
  gfracvstim->SetTitle("Fraction of Good Anodes vs. run");
  gfrac3vstim->SetTitle("Fraction of Good Anodes vs. run");
  gfrac4vstim->SetTitle("Fraction of Good Anodes vs. run");
  gzslow3vstim->SetTitle("Low ZS thr vs run");
  gzslow4vstim->SetTitle("Low ZS thr vs run");
  gzshi3vstim->SetTitle("High ZS thr vs run");
  gzshi4vstim->SetTitle("High ZS thr vs run");


  Char_t filnam[200],filnamalien[200];
  Int_t iPoint=0;
  Int_t nrun,nrun2,nv,ns;
  Double_t minzs=256.;
  Double_t maxzs=0.;

  while(!feof(listruns)){

    hbase->Reset();
    hnoise->Reset();
    hgain->Reset();
    hchstatus->Reset();
    hchstatus3->Reset();
    hchstatus4->Reset();
    hlowzsthr3->Reset();
    hlowzsthr4->Reset();
    hhizsthr3->Reset();
    hhizsthr4->Reset();

    fscanf(listruns,"%s\n",filnam);    
    Char_t directory[100];
    sprintf(directory,"/alice/data/%d",year);
    if(!strstr(filnam,directory)) continue;
    sscanf(filnam,"/alice/data/%d/OCDB/ITS/Calib/CalibSDD/Run%d_%d_v%d_s%d.root",&year,&nrun,&nrun2,&nv,&ns);
    if(year==2009 && (nrun<85639 && nrun2> 85639)) continue; // protection for files with swapped ladders 4-5 of layer 3 
    if(year==2009 && (nrun>100000 && nv< 184)) continue; // protection for files with swapped ladder 0-1 of layer 4
    if(year==2010 && (nrun>=114603 && nv< 98)) continue; // protection for files without treatment of masked hybrids 
    if(year==2011 && (nrun>=145349 && nrun<=148978) && nrun2> 148978) continue; // protection for files affected by problem in second DA
    if(year==2011 && nrun==156856) continue;
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
      if(selectedMod>=240 && (iMod+240)!=selectedMod) continue;
      cal=(AliITSCalibrationSDD*)calSDD->At(iMod);
      if(cal==0) continue;
      Float_t zshi0=cal->GetZSHighThreshold(0);
      Float_t zshi1=cal->GetZSHighThreshold(1);
      Float_t zslo0=cal->GetZSLowThreshold(0);
      Float_t zslo1=cal->GetZSLowThreshold(1);
      if(zshi0>maxzs) maxzs=zshi0;
      if(zshi1>maxzs) maxzs=zshi1;
      if(zslo0<minzs) minzs=zslo0;
      if(zslo1<minzs) minzs=zslo1;
      if(!cal->IsBad()){
	if(iMod<84){
	  hlowzsthr3->Fill(zslo0);
	  hlowzsthr3->Fill(zslo1);
	  hhizsthr3->Fill(zshi0);
	  hhizsthr3->Fill(zshi1);
	  hlowzsthr3VsRun->Fill(iPoint,zslo0);
	  hlowzsthr3VsRun->Fill(iPoint,zslo1);
	  hhizsthr3VsRun->Fill(iPoint,zshi0);
	  hhizsthr3VsRun->Fill(iPoint,zshi1);
	}else{
	  hlowzsthr4->Fill(zslo0);
	  hlowzsthr4->Fill(zslo1);
	  hhizsthr4->Fill(zshi0);
	  hhizsthr4->Fill(cal->GetZSHighThreshold(1));
	  hlowzsthr4VsRun->Fill(iPoint,zslo0);
	  hlowzsthr4VsRun->Fill(iPoint,zslo1);
	  hhizsthr4VsRun->Fill(iPoint,zshi0);
	  hhizsthr4VsRun->Fill(iPoint,zshi1);
	}
      }
      hlowzsthr3VsRun->GetXaxis()->SetBinLabel(iPoint+1,Form("%d",nrun));
      hhizsthr3VsRun->GetXaxis()->SetBinLabel(iPoint+1,Form("%d",nrun));
      hlowzsthr4VsRun->GetXaxis()->SetBinLabel(iPoint+1,Form("%d",nrun));
      hhizsthr4VsRun->GetXaxis()->SetBinLabel(iPoint+1,Form("%d",nrun));
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
    if(selectedMod==-1 && (Int_t)hbase->GetEntries()==0) continue;
    gbasevstim->SetPoint(iPoint,(Double_t)nrun,hbase->GetMean());
    gbasevstim->SetPointError(iPoint,0.,hbase->GetRMS());
    gnoisevstim->SetPoint(iPoint,(Double_t)nrun,hnoise->GetMean());
    gnoisevstim->SetPointError(iPoint,0.,hnoise->GetRMS());
    ggainvstim->SetPoint(iPoint,(Double_t)nrun,hgain->GetMean());
    ggainvstim->SetPointError(iPoint,0.,hgain->GetRMS());
    gstatvstim->SetPoint(iPoint,(Double_t)nrun,hchstatus->GetBinContent(2));
    Float_t normMod=260.;
    if(selectedMod!=-1) normMod=1.;
    gfracvstim->SetPoint(iPoint,(Double_t)nrun,hchstatus->GetBinContent(2)/normMod/512.);
    gfrac3vstim->SetPoint(iPoint,(Double_t)nrun,hchstatus3->GetBinContent(2)/84./512.);
    gfrac4vstim->SetPoint(iPoint,(Double_t)nrun,hchstatus4->GetBinContent(2)/176./512.);
    gzslow3vstim->SetPoint(iPoint,(Double_t)nrun,hlowzsthr3->GetMean());
    //  gzslow3vstim->SetPointError(iPoint,0.,hlowzsthr3->GetRMS());
    gzslow4vstim->SetPoint(iPoint,(Double_t)nrun,hlowzsthr4->GetMean());
    //  gzslow4vstim->SetPointError(iPoint,0.,hlowzsthr4->GetRMS());
    gzshi3vstim->SetPoint(iPoint,(Double_t)nrun,hhizsthr3->GetMean());
    //  gzshi3vstim->SetPointError(iPoint,0.,hhizsthr3->GetRMS());
    gzshi4vstim->SetPoint(iPoint,(Double_t)nrun,hhizsthr4->GetMean());
    // gzshi4vstim->SetPointError(iPoint,0.,hhizsthr4->GetRMS());    
    iPoint++;
    f->Close();
  }
  hlowzsthr3VsRun->GetXaxis()->SetRange(1,iPoint);
  hlowzsthr4VsRun->GetXaxis()->SetRange(1,iPoint);
  hhizsthr3VsRun->GetXaxis()->SetRange(1,iPoint);
  hhizsthr4VsRun->GetXaxis()->SetRange(1,iPoint);
  hlowzsthr3VsRun->GetYaxis()->SetRangeUser(minzs-5,maxzs+5);
  hlowzsthr4VsRun->GetYaxis()->SetRangeUser(minzs-5,maxzs+5);
  hhizsthr3VsRun->GetYaxis()->SetRangeUser(minzs-5,maxzs+5);
  hhizsthr4VsRun->GetYaxis()->SetRangeUser(minzs-5,maxzs+5);

  TFile *ofil=new TFile(Form("Calib%dVsTime.root",year),"recreate");
  gbasevstim->Write();
  gnoisevstim->Write();
  ggainvstim->Write();
  gstatvstim->Write();
  gfracvstim->Write();
  gfrac3vstim->Write();
  gfrac4vstim->Write();
  hlowzsthr3VsRun->Write();
  hlowzsthr4VsRun->Write();
  hhizsthr3VsRun->Write();
  hhizsthr4VsRun->Write();
  gzslow3vstim->Write();
  gzslow4vstim->Write();
  gzshi3vstim->Write();
  gzshi4vstim->Write();
  ofil->Close();

  TCanvas* cbase=new TCanvas("cbase","Baselines");
  gPad->SetTickx();
  gPad->SetTicky();
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
  gPad->SetTickx();
  gPad->SetTicky();
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
  gPad->SetTickx();
  gPad->SetTicky();
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
  gPad->SetTickx();
  gPad->SetTicky();
  gstatvstim->SetFillColor(kOrange-2);
  gstatvstim->SetMarkerStyle(20);
  gstatvstim->Draw("AP3");
  gstatvstim->Draw("PLXSAME");
  if(selectedMod==-1){
    gstatvstim->SetMinimum(100000.);
    gstatvstim->SetMaximum(133000.);
  }else{
    gstatvstim->SetMinimum(0.);
    gstatvstim->SetMaximum(512.);
  }
  gstatvstim->GetXaxis()->SetTitle("Run number");
  if(selectedMod==-1){
    gstatvstim->GetYaxis()->SetTitle("Number of good anodes in acquisition");
  }else{
    gstatvstim->GetYaxis()->SetTitle(Form("Number of good anodes in od %d",selectedMod));
  }
  cstatus->SaveAs(Form("GoodAnodesRun%d.gif",year));

  TCanvas* cfrac=new TCanvas("cfrac","Fraction of Good");
  gPad->SetTickx();
  gPad->SetTicky();
  gfracvstim->SetMarkerStyle(20);
  gfrac3vstim->SetMarkerStyle(22);
  gfrac3vstim->SetMarkerColor(2);
  gfrac3vstim->SetLineColor(2);
  gfrac4vstim->SetMarkerStyle(23);
  gfrac4vstim->SetMarkerColor(4);
  gfrac4vstim->SetLineColor(4);
  gfracvstim->Draw("APL");
  gfracvstim->SetMinimum(0.7);
  gfracvstim->SetMaximum(1.05);
  gfracvstim->GetXaxis()->SetTitle("Run number");
  if(selectedMod==-1){
    gfracvstim->GetYaxis()->SetTitle("Fraction of good anodes in acquisition");
    gfrac3vstim->Draw("PLSAME");
    gfrac4vstim->Draw("PLSAME");
  
    TLegend* leg=new TLegend(0.2,0.15,0.45,0.35);
    leg->SetFillColor(0);
    TLegendEntry* entr=leg->AddEntry(gfrac3vstim,"Layer 3","P");
    entr->SetTextColor(2);
    entr=leg->AddEntry(gfrac4vstim,"Layer 4","P");
    entr->SetTextColor(4);
    entr=leg->AddEntry(gfracvstim,"All","P");
    entr->SetTextColor(1);
    leg->Draw();
  }else{
    gfracvstim->GetYaxis()->SetTitle(Form("Fraction of good anodes in mod %d",selectedMod));
  }
  cfrac->SaveAs(Form("FractionGoodRun%d.gif",year));

  TCanvas* czsd=new TCanvas("czsd","Zero Supp Details",1500,800);
  czsd->Divide(2,2);
  czsd->cd(1);
  gPad->SetLogz();
  hlowzsthr3VsRun->Draw("colz");
  TLatex* tlay3=new TLatex(0.18,0.82,"Layer 3");
  tlay3->SetNDC();
  tlay3->SetTextFont(63);
  tlay3->SetTextSize(24);
  tlay3->Draw();
  czsd->cd(2);
  gPad->SetLogz();
  hlowzsthr4VsRun->Draw("colz");
  TLatex* tlay4=new TLatex(0.18,0.82,"Layer 4");
  tlay4->SetNDC();
  tlay4->SetTextFont(63);
  tlay4->SetTextSize(24);
  tlay4->Draw();
  czsd->cd(3);
  gPad->SetLogz();
  hhizsthr3VsRun->Draw("colz");
  tlay3->Draw();
  czsd->cd(4);
  gPad->SetLogz();
  hhizsthr4VsRun->Draw("colz");
  tlay4->Draw();
  czsd->SaveAs(Form("ZeroSuppDistThrRun%d.gif",year));

  TCanvas* czs=new TCanvas("czs","Zero Supp");
  gPad->SetTickx();
  gPad->SetTicky();
  gzslow3vstim->SetMarkerStyle(22);
  gzslow3vstim->SetMarkerColor(2);
  gzslow3vstim->SetLineColor(2);
  gzslow4vstim->SetMarkerStyle(23);
  gzslow4vstim->SetMarkerColor(4);
  gzslow4vstim->SetLineColor(4);
  gzshi3vstim->SetMarkerStyle(26);
  gzshi3vstim->SetMarkerColor(kRed+1);
  gzshi3vstim->SetLineColor(kRed+1);
  gzshi4vstim->SetMarkerStyle(32);
  gzshi4vstim->SetMarkerColor(kBlue+1);
  gzshi4vstim->SetLineColor(kBlue+1);

  gzslow3vstim->SetMinimum(20.);
  gzslow3vstim->SetMaximum(35.);
  gzslow3vstim->Draw("APL");
  gzslow3vstim->GetXaxis()->SetTitle("Run number");
  gzslow3vstim->GetYaxis()->SetTitle("<Zero suppression threshold>");
  gzslow4vstim->Draw("SAMEPL");
  gzshi3vstim->Draw("SAMEPL");
  gzshi4vstim->Draw("SAMEPL");
  if(selectedMod==-1){
    TLegend* legz=new TLegend(0.2,0.15,0.45,0.35);
    legz->SetFillColor(0);
    legz->AddEntry(gzslow3vstim,"Low Thr - Layer 3","P");
    legz->AddEntry(gzslow4vstim,"Low Thr - Layer 4","P");
    legz->AddEntry(gzshi3vstim,"High Thr - Layer 3","P");
    legz->AddEntry(gzshi4vstim,"High Thr - Layer 4","P");
    legz->Draw();
  }
  czs->SaveAs(Form("ZeroSuppMeanThrRun%d.gif",year));
}
