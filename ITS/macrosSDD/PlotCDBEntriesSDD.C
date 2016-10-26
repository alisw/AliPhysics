#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TF1.h>
#include <TGraph.h>
#include <TExec.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TLatex.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TPaveStats.h>
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliITSgeomTGeo.h"
#include "AliITSDDLModuleMapSDD.h"
#include "AliITSresponseSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSDriftSpeedArraySDD.h"
#include "AliITSDriftSpeedSDD.h"
#endif

void MakePalette();
void PlotCalib(AliCDBEntry *pedpul, Bool_t optVerbose);
void PlotDriftSpeed(AliCDBEntry *inject, Bool_t optVerbose);
void PlotOfflineCalib(AliCDBEntry *resp);

void MakePalette(){
  Int_t palette[3]={kGray,2,3};  
  gStyle->SetPalette(3,palette);
}


void PlotCDBEntriesSDD(Int_t nrun=238145, Bool_t optVerbose=kFALSE){
  TGrid::Connect("alien:");
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(nrun);
  AliCDBEntry *pedpul = man->Get("ITS/Calib/CalibSDD");
  AliCDBEntry *inject = man->Get("ITS/Calib/DriftSpeedSDD");
  AliCDBEntry *resp = man->Get("ITS/Calib/RespSDD");
  PlotCalib(pedpul,optVerbose);
  PlotDriftSpeed(inject,optVerbose);
  PlotOfflineCalib(resp);
}


void PlotCalib(AliCDBEntry *pedpul, Bool_t optVerbose){

  printf("====== PARAMETERS FROM PEDESTAL+PULSER RUN ======\n");

  TH2I* hlay3=new TH2I("hlay3","Layer 3",12,-0.5,5.5,14,-0.5,13.5);
  hlay3->GetXaxis()->SetTitle("Detector");
  hlay3->GetYaxis()->SetTitle("Ladder");
  hlay3->GetXaxis()->SetTickLength(0);
  hlay3->GetYaxis()->SetTickLength(0);
  hlay3->SetStats(0);
  hlay3->SetMinimum(-1);
  TH2I* hlay4=new TH2I("hlay4","Layer 4",16,-0.5,7.5,22,-0.5,21.5);
  hlay4->GetXaxis()->SetTitle("Detector");
  hlay4->GetYaxis()->SetTitle("Ladder");
  hlay4->GetXaxis()->SetTickLength(0);
  hlay4->GetYaxis()->SetTickLength(0);
  hlay4->GetYaxis()->SetTitle("Ladder");
  hlay4->SetStats(0);
  hlay4->SetMinimum(-1);
  TH2I* hdeadlay3=new TH2I("hdlay3","Layer 3",6,-0.5,5.5,14,-0.5,13.5);
  hdeadlay3->GetXaxis()->SetTitle("Detector");
  hdeadlay3->GetYaxis()->SetTitle("Ladder");
  hdeadlay3->GetXaxis()->SetTickLength(0);
  hdeadlay3->GetYaxis()->SetTickLength(0);
  hdeadlay3->SetStats(0);
  hdeadlay3->SetMinimum(-1.);
  TH2I* hdeadlay4=new TH2I("hdlay4","Layer 4",8,-0.5,7.5,22,-0.5,21.5);
  hdeadlay4->GetXaxis()->SetTitle("Detector");
  hdeadlay4->GetYaxis()->SetTitle("Ladder");
  hdeadlay4->GetXaxis()->SetTickLength(0);
  hdeadlay4->GetYaxis()->SetTickLength(0);
  hdeadlay4->GetYaxis()->SetTitle("Ladder");
  hdeadlay4->SetStats(0);
  hdeadlay4->SetMinimum(-1.);

  TObjArray *calSDD = (TObjArray *)pedpul->GetObject();
  TH1F* hmodstatus=new TH1F("hmodstatus","",260,0.5,260.5);
  TH1F* hnbadch=new TH1F("hnbadch","",260,0.5,260.5);
  TH1F* hbase=new TH1F("hbase","",60,0.5,120.5);
  TH2F* hbasemod=new TH2F("hbasemod","",260,239.5,499.5,50,0.,100.);
  TH1F* hnoise=new TH1F("hnoise","",100,0.,7.);
  TH2F* hnoisemod=new TH2F("hnoisemod","",260,239.5,499.5,50,0.,10.);
  TH1F* hgain=new TH1F("hgain","",100,0.,4.);
  TH2F* hgainmod=new TH2F("hgainmod","",260,239.5,499.5,50,0.,4.);
  TH1F* hchstatus=new TH1F("hchstatus","",2,-0.5,1.5);

  AliITSDDLModuleMapSDD* dmap=new AliITSDDLModuleMapSDD();
  dmap->SetJun09Map();
  TH2I *hddlcarlos=new TH2I("hddlcarlos","",24,-0.5,11.5,24,-0.5,23.5);
  hddlcarlos->GetXaxis()->SetTitle("Module");
  hddlcarlos->GetYaxis()->SetTitle("DDL");
  hddlcarlos->GetXaxis()->SetTickLength(0);
  hddlcarlos->GetYaxis()->SetTickLength(0);
  hddlcarlos->SetStats(0);
  hddlcarlos->SetMinimum(-1.);

  AliITSCalibrationSDD *cal;
  Int_t badModCounter3=0;
  Int_t badModCounter4=0;
  Int_t badHybridCounter3=0;
  Int_t badHybridCounter4=0;
  Int_t badAnodeCounter3=0;
  Int_t badAnodeCounter4=0;
  Int_t badAnodeCounterGoodMod3=0;
  Int_t badAnodeCounterGoodMod4=0;
  Int_t badAnodeCounterGoodHybrid3=0;
  Int_t badAnodeCounterGoodHybrid4=0;
  Int_t badAnodeCounterGoodModAndChip3=0;
  Int_t badAnodeCounterGoodModAndChip4=0;
  Int_t badChipCounter3=0;
  Int_t badChipCounter4=0;
  for(Int_t i=0; i<260; i++){
    cal=(AliITSCalibrationSDD*)calSDD->At(i);
    if(cal==0) continue;
    if(optVerbose) printf("Module %d (%d)   status = ",i,i+240);
    Int_t lay,lad,det;
    AliITSgeomTGeo::GetModuleId(i+240,lay,lad,det);
    Int_t index=1+(det-1)*2;
    Int_t ddl,carlos;
    dmap->FindInDDLMap(i+240,ddl,carlos);
    Int_t index2=1+carlos*2;
    ddl+=1;
    if(cal->IsBad()){ 
      if(optVerbose) printf("BAD\t");
      hddlcarlos->SetBinContent(index2,ddl,0);
      hddlcarlos->SetBinContent(index2+1,ddl,0);
      if(lay==3){ 
	badModCounter3++;
	badHybridCounter3+=2;
	hlay3->SetBinContent(index,lad,0);
	hlay3->SetBinContent(index+1,lad,0);
      }else if(lay==4){ 
	badModCounter4++;
	badHybridCounter4+=2;
	hlay4->SetBinContent(index,lad,0);
	hlay4->SetBinContent(index+1,lad,0);
      }
      hmodstatus->SetBinContent(i+1,0);
    }else{ 
      if(optVerbose) printf("OK\t");
      hmodstatus->SetBinContent(i+1,1);
      if(lay==3){ 
	badAnodeCounterGoodMod3+=cal->GetDeadChannels();
	if(cal->IsChipBad(0) && cal->IsChipBad(1) && cal->IsChipBad(2) && cal->IsChipBad(3)){
	  hlay3->SetBinContent(index,lad,0);
	  hddlcarlos->SetBinContent(index2,ddl,0);
	  badHybridCounter3++;
	}else{
	  hlay3->SetBinContent(index,lad,1);
	  hddlcarlos->SetBinContent(index2,ddl,1);
	  for(Int_t iAn=0; iAn<256; iAn++){
	    if(cal->IsBadChannel(iAn)) badAnodeCounterGoodHybrid3++;
	  }
	}
	if(cal->IsChipBad(4) && cal->IsChipBad(5) && cal->IsChipBad(6) && cal->IsChipBad(7)){
	  hlay3->SetBinContent(index+1,lad,0);
	  hddlcarlos->SetBinContent(index2+1,ddl,0);
	  badHybridCounter3++;
	}else{
	  hlay3->SetBinContent(index+1,lad,1);
	  hddlcarlos->SetBinContent(index2+1,ddl,1);
	  for(Int_t iAn=256; iAn<512; iAn++){
	    if(cal->IsBadChannel(iAn)) badAnodeCounterGoodHybrid3++;
	  }
	}
      }else{ 
	badAnodeCounterGoodMod4+=cal->GetDeadChannels();
	if(cal->IsChipBad(0) && cal->IsChipBad(1) && cal->IsChipBad(2) && cal->IsChipBad(3)){
	  hlay4->SetBinContent(index,lad,0);
	  hddlcarlos->SetBinContent(index2,ddl,0);
	  badHybridCounter4++;
	}else{
	  hlay4->SetBinContent(index,lad,1);
	  hddlcarlos->SetBinContent(index2,ddl,1);
	  for(Int_t iAn=0; iAn<256; iAn++){
	    if(cal->IsBadChannel(iAn)) badAnodeCounterGoodHybrid4++;
	  }
	}
	if(cal->IsChipBad(4) && cal->IsChipBad(5) && cal->IsChipBad(6) && cal->IsChipBad(7)){
	  hlay4->SetBinContent(index+1,lad,0);
	  hddlcarlos->SetBinContent(index2+1,ddl,0);
	  badHybridCounter4++;
	}else{
	  hlay4->SetBinContent(index+1,lad,1);
	  hddlcarlos->SetBinContent(index2+1,ddl,1);
	  for(Int_t iAn=256; iAn<512; iAn++){
	    if(cal->IsBadChannel(iAn)) badAnodeCounterGoodHybrid4++;
	  }
	}
      }
    }
    if(optVerbose) printf("   Chip Status (0=OK, 1=BAD): ");  
    for(Int_t ic=0; ic<8;ic++){ 
      if(optVerbose) printf("%d ",cal->IsChipBad(ic));
      if(cal->IsChipBad(ic) && !cal->IsBad()){ 
	if(i<84) badChipCounter3++;
	else badChipCounter4++;
      }
    }
    if(optVerbose){
      printf(" # bad anodes = %d  ",cal->GetDeadChannels());
      if(cal->IsAMAt20MHz()) printf("      20 MHz sampling");
      else printf("      40 MHz sampling");
      printf(" Threshold L %d %d H %d %d\n",cal->GetZSLowThreshold(0),cal->GetZSLowThreshold(1),cal->GetZSHighThreshold(0),cal->GetZSHighThreshold(1));
    }
    if(i<84) badAnodeCounter3+=cal->GetDeadChannels();
    else badAnodeCounter4+=cal->GetDeadChannels();
    hnbadch->SetBinContent(i+1,cal->GetDeadChannels());
    if(lay==3) hdeadlay3->SetBinContent(det,lad,cal->GetDeadChannels());
    if(lay==4) hdeadlay4->SetBinContent(det,lad,cal->GetDeadChannels());
    for(Int_t iAn=0; iAn<512; iAn++){
      Int_t ic=cal->GetChip(iAn);
      if(!cal->IsChipBad(ic) && !cal->IsBad() && cal->IsBadChannel(iAn)){ 
	if(i<84) badAnodeCounterGoodModAndChip3++;
	else badAnodeCounterGoodModAndChip4++;
      }
      Float_t base=cal->GetBaseline(iAn);
      Float_t noise=cal->GetNoiseAfterElectronics(iAn);
      Float_t gain=cal->GetChannelGain(iAn);
      if(cal->IsBadChannel(iAn)) hchstatus->Fill(0);
      if(!cal->IsBadChannel(iAn) && !cal->IsChipBad(ic) && !cal->IsBad() ){
	hbase->Fill(base);
	hbasemod->Fill(i+240,base);
	hchstatus->Fill(1);
	hnoise->Fill(noise);
	hnoisemod->Fill(i+240,noise);
	hgain->Fill(gain);
	hgainmod->Fill(i+240,gain);
      }
    }
  }
  Int_t totbad3=badModCounter3*512+badChipCounter3*64+badAnodeCounterGoodModAndChip3;
  Int_t tot3=6*14*512;
  Float_t fracbad3=(Float_t)totbad3/(Float_t)tot3;
  Float_t fracgm3=(Float_t)(84.-badModCounter3)/84.;
  Float_t fracgm4=(Float_t)(176.-badModCounter4)/176.;
  Float_t fracgh3=(Float_t)(84.*2.-badHybridCounter3)/84./2.;
  Float_t fracgh4=(Float_t)(176.*2-badHybridCounter4)/176./2.;
  Float_t fraccgm3=0.;
  if(badModCounter3!=84){
    fraccgm3=1.-(Float_t)(badAnodeCounterGoodModAndChip3+badChipCounter3*64)/(512.*(Float_t)(84.-badModCounter3));
  }
  Float_t fraccgm4=0.;
  if(badModCounter4!=176){
    fraccgm4=1.-(Float_t)(badAnodeCounterGoodModAndChip4+badChipCounter4*64)/(512.*(Float_t)(176.-badModCounter4));
  }
  Float_t fraccgh3=0.;
  if(badHybridCounter3!=(84*2)){
    fraccgh3=1.-(Float_t)badAnodeCounterGoodHybrid3/(256.*(84.*2.-badHybridCounter3));
  }
  Float_t fraccgh4=0.;
  if(badHybridCounter4!=(176*2)){
    fraccgh4=1.-(Float_t)badAnodeCounterGoodHybrid4/(256.*(176.*2.-badHybridCounter4));
  }
  Int_t totbad4=badModCounter4*512+badChipCounter4*64+badAnodeCounterGoodModAndChip4;
  Int_t tot4=8*22*512;
  Float_t fracbad4=(Float_t)totbad4/(Float_t)tot4;
  Float_t fractot=(Float_t)(totbad3+totbad4)/(Float_t)(tot3+tot4);

  printf("---- Layer 3 ----\n");
  printf("# of bad modules                      = %d\n",badModCounter3);
  printf("# of bad chips in good modules        = %d\n",badChipCounter3);
  printf("# of bad hybrids                      = %d\n",badHybridCounter3);
  printf("# of bad anodes in good hybrids       = %d\n",badAnodeCounterGoodHybrid3);
  printf("Fraction of Good modules=%f\n",fracgm3);
  printf("Fraction of Good hybrids=%f\n",fracgh3);
  printf("Fraction of good anodes in good modules = %f\n",fraccgm3);
  printf("Fraction of good anodes in good hybrids = %f\n",fraccgh3);
  printf("Fraction of bads (anodes+chips+mod)     = %f\n",fracbad3);
  printf("---- Layer 4 ----\n");
  printf("# of bad modules                      = %d\n",badModCounter4);
  printf("# of bad chips in good modules        = %d\n",badChipCounter4);
  printf("# of bad hybrids                      = %d\n",badHybridCounter4);
  printf("# of bad anodes in good hybrids       = %d\n",badAnodeCounterGoodHybrid4);
  printf("Fraction of Good modules=%f\n",fracgm4);
  printf("Fraction of Good hybrids=%f\n",fracgh4);
  printf("Fraction of good anodes in good modules = %f\n",fraccgm4);
  printf("Fraction of good anodes in good hybrids = %f\n",fraccgh4);
  printf("Fraction of bads (anodes+chips+mod)     = %f\n",fracbad4);
  printf("---- Total   ----\n");
  printf("# of bad modules                      = %d\n",badModCounter3+badModCounter4);
  printf("# of bad chips in good modules        = %d\n",badChipCounter3+badChipCounter4);
  printf("# of bad anodes in good modules+chips = %d\n",badAnodeCounterGoodModAndChip3+badAnodeCounterGoodModAndChip4);
  printf("Fraction of bads (anodes+chips+mod)   = %f\n",fractot);
  printf("---------------------------------------------------\n");
  

  TLine* lin=new TLine(0,0,0,23);  
  TExec *ex1 = new TExec("ex1","MakePalette();");
  TExec *ex2 = new TExec("ex2","gStyle->SetPalette(1);");

  TCanvas* clay=new TCanvas("clay","Layer status",900,600);
  clay->Divide(2,1);
  clay->cd(1);
  hlay3->Draw("col");
  ex1->Draw();
  hlay3->DrawCopy("col same");
  for(Int_t i=0;i<6;i++){
    lin->SetY1(-0.5);
    lin->SetY2(13.5);
    lin->SetX1(i+0.5);
    lin->SetX2(i+0.5);
    lin->DrawClone();
  }
  for(Int_t i=0;i<14;i++){
    lin->SetX1(-0.5);
    lin->SetX2(5.5);
    lin->SetY1(i+0.5);
    lin->SetY2(i+0.5);
    lin->DrawClone();
  }
  clay->cd(2);
  hlay4->DrawCopy("col");
  for(Int_t i=0;i<8;i++){
    lin->SetY1(-0.5);
    lin->SetY2(21.5);
    lin->SetX1(i+0.5);
    lin->SetX2(i+0.5);
    lin->DrawClone();
  }
  for(Int_t i=0;i<22;i++){
    lin->SetX1(-0.5);
    lin->SetX2(7.5);
    lin->SetY1(i+0.5);
    lin->SetY2(i+0.5);
    lin->DrawClone();
  }


  TCanvas* cddl=new TCanvas("cddl","DDL status",800,800);
  hddlcarlos->Draw("col");
  ex1->Draw();
  hddlcarlos->DrawCopy("col same");
  for(Int_t i=0;i<12;i++){
    lin->SetY1(-0.5);
    lin->SetY2(23.5);
    lin->SetX1(i+0.5);
    lin->SetX2(i+0.5);
    lin->DrawClone();
  }
  for(Int_t i=0;i<24;i++){
    lin->SetX1(-0.5);
    lin->SetX2(11.5);
    lin->SetY1(i+0.5);
    lin->SetY2(i+0.5);
    lin->DrawClone();
  }
 

  TCanvas *c0b=new TCanvas("c0b","Bad Channels",900,600);
  c0b->Divide(2,1);
  c0b->cd(1);
  hdeadlay3->DrawCopy("colz");
  ex2->Draw();
  hdeadlay3->DrawCopy("colz same");
  for(Int_t i=0;i<6;i++){
    lin->SetY1(-0.5);
    lin->SetY2(13.5);
    lin->SetX1(i+0.5);
    lin->SetX2(i+0.5);
    lin->DrawClone();
  }
  for(Int_t i=0;i<14;i++){
    lin->SetX1(-0.5);
    lin->SetX2(5.5);
    lin->SetY1(i+0.5);
    lin->SetY2(i+0.5);
    lin->DrawClone();
  }
  c0b->cd(2);
  hdeadlay4->DrawCopy("colz");
  ex2->Draw();
  hdeadlay4->DrawCopy("colz same");
  for(Int_t i=0;i<8;i++){
    lin->SetY1(-0.5);
    lin->SetY2(21.5);
    lin->SetX1(i+0.5);
    lin->SetX2(i+0.5);
    lin->DrawClone();
  }
  for(Int_t i=0;i<22;i++){
    lin->SetX1(-0.5);
    lin->SetX2(7.5);
    lin->SetY1(i+0.5);
    lin->SetY2(i+0.5);
    lin->DrawClone();
  }


  
  TCanvas *c1=new TCanvas("c1","Anode calibration",800,800);
  c1->Divide(2,2);
  c1->cd(1);
  hbase->Draw();
  hbase->GetXaxis()->SetTitle("Baseline after equalization");
  hbase->GetXaxis()->CenterTitle();  
  c1->cd(2);
  hnoise->Draw(); 
  hnoise->GetXaxis()->SetTitle("Noise");
  hnoise->GetXaxis()->CenterTitle();  
  c1->cd(3);
  hgain->Draw();
  hgain->GetXaxis()->SetTitle("Gain");
  hgain->GetXaxis()->CenterTitle();  
  c1->cd(4);
  hchstatus->Draw();
  hchstatus->GetXaxis()->SetTitle("Anode status (0=bad, 1=OK)");
  hchstatus->GetXaxis()->CenterTitle();

  TCanvas *c1m=new TCanvas("c1m","Calib. vs. mod",1000,800);
  c1m->Divide(2,2);
  c1m->cd(1);
  gPad->SetRightMargin(0.14);
  hbasemod->SetStats(0);
  hbasemod->Draw("colz"); 
  hbasemod->GetXaxis()->SetTitle("Module Number");
  hbasemod->GetYaxis()->SetTitle("Baseline");
  c1m->cd(2);
  gPad->SetRightMargin(0.14);
  hnoisemod->SetStats(0);
  hnoisemod->Draw("colz"); 
  hnoisemod->GetXaxis()->SetTitle("Module Number");
  hnoisemod->GetYaxis()->SetTitle("Noise");
  c1m->cd(3);
  gPad->SetRightMargin(0.14);
  hgainmod->SetStats(0);
  hgainmod->Draw("colz");
  hgainmod->GetXaxis()->SetTitle("Module Number");
  hgainmod->GetYaxis()->SetTitle("Gain");
  c1m->cd(4);
  hnbadch->Scale(1/512.);
  hnbadch->SetMarkerStyle(20);
  hnbadch->SetMarkerSize(0.8);
  hnbadch->SetStats(0);
  hnbadch->Draw("P");
  hnbadch->GetXaxis()->SetTitle("Module number");   
  hnbadch->GetYaxis()->SetTitle("Fraction of bad anodes");
  return;
}


void PlotDriftSpeed(AliCDBEntry *inject, Bool_t optVerbose){

  printf("====== PARAMETERS FROM INJECTOR RUN ======\n");

  TObjArray *drspSDD = (TObjArray *)inject->GetObject();
  AliITSDriftSpeedArraySDD *vdriftarr0;
  AliITSDriftSpeedArraySDD *vdriftarr1;

  TH2F* hlay3=new TH2F("hinjlay3","Injector Status Layer 3",12,-0.5,5.5,14,-0.5,13.5);
  hlay3->GetXaxis()->SetTitle("Detector");
  hlay3->GetYaxis()->SetTitle("Ladder");
  hlay3->GetXaxis()->SetTickLength(0);
  hlay3->GetYaxis()->SetTickLength(0);
  hlay3->SetStats(0);
  hlay3->SetMinimum(-0.01);
  hlay3->SetMaximum(7.);
  TH2F* hlay4=new TH2F("hinjlay4","Injector Status Layer 4",16,-0.5,7.5,22,-0.5,21.5);
  hlay4->GetXaxis()->SetTitle("Detector");
  hlay4->GetYaxis()->SetTitle("Ladder");
  hlay4->GetXaxis()->SetTickLength(0);
  hlay4->GetYaxis()->SetTickLength(0);
  hlay4->GetYaxis()->SetTitle("Ladder");
  hlay4->SetStats(0);
  hlay4->SetMinimum(-0.01);
  hlay4->SetMaximum(7.);

  TGraph *tempvsmod0=new TGraph(0);
  TGraph *tempvsmod1=new TGraph(0);
  TGraph *vvsmod0=new TGraph(0);
  TGraph *vvsmod1=new TGraph(0);
  TGraph *poldegvsmod0=new TGraph(0); 
  TGraph *poldegvsmod1=new TGraph(0); 
  TGraph *anmaxvsmod0=new TGraph(0); 
  TGraph *anmaxvsmod1=new TGraph(0); 
  TGraph *dvcevsmod0=new TGraph(0);
  TGraph *dvcevsmod1=new TGraph(0);
  TGraph *dveevsmod0=new TGraph(0);
  TGraph *dveevsmod1=new TGraph(0);

  char tit0[100];
  Int_t nrun=0;
  sprintf(tit0,"Temperature vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Temperature vs. mod. number - Run %d",nrun);
  tempvsmod0->SetTitle(tit0);
  tempvsmod1->SetTitle(tit0);

  sprintf(tit0,"Drift Speed vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Drift Speed vs. mod. number - Run %d",nrun);
  vvsmod0->SetTitle(tit0);
  vvsmod1->SetTitle(tit0);

  sprintf(tit0,"Degree of poly fit vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Degree of poly fit vs. mod. number - Run %d",nrun);
  poldegvsmod0->SetTitle(tit0);
  poldegvsmod1->SetTitle(tit0);

  sprintf(tit0,"Anode with max. vdrift vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Anode with max. vdrift vs. mod. number - Run %d",nrun);
  anmaxvsmod0->SetTitle(tit0);
  anmaxvsmod1->SetTitle(tit0);

  sprintf(tit0,"Delta Vdrift 128-0 vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Delta Vdrift 128-0 vs. mod. number - Run %d",nrun);
  dvcevsmod0->SetTitle(tit0);
  dvcevsmod1->SetTitle(tit0);

  sprintf(tit0,"Delta Vdrift 256-0 vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Delta Vdrift 256-0 vs. mod. number - Run %d",nrun);
  dveevsmod0->SetTitle(tit0);
  dveevsmod1->SetTitle(tit0);

  TF1* fPoly=new TF1("fPoly","pol3",0.,256.);
  Char_t tit[100];
  Int_t cntpad = 0;
  Int_t iGoodInj=0;
  Int_t iRescaledSpeed=0;
  Int_t iAverSpeed=0;
  TLatex* tleft=new TLatex(0.2,0.82,"Side 0");
  tleft->SetNDC();
  tleft->SetTextColor(1);
  TLatex* tright=new TLatex(0.2,0.75,"Side 1");
  tright->SetNDC();
  tright->SetTextColor(2);

  for(Int_t i=0; i<260; i++){
    Int_t iMod=i+240;
    Int_t lay,lad,det;
    AliITSgeomTGeo::GetModuleId(iMod,lay,lad,det);
    Int_t i0=2*i;
    Int_t i1=1+2*i;
    vdriftarr0=(AliITSDriftSpeedArraySDD*)drspSDD->At(i0);
    vdriftarr1=(AliITSDriftSpeedArraySDD*)drspSDD->At(i1);
    AliITSDriftSpeedSDD* vdrift0=0x0;
    if(vdriftarr0) vdrift0=vdriftarr0->GetDriftSpeedObject(0);
    AliITSDriftSpeedSDD* vdrift1=0x0;
    if(vdriftarr1) vdrift1=vdriftarr1->GetDriftSpeedObject(0);

    Int_t statusInj0=vdriftarr0->GetInjectorStatus();
    Int_t statusInj1=vdriftarr1->GetInjectorStatus();
    if(statusInj0>1) iGoodInj++;
    else if(statusInj0==1) iRescaledSpeed++;
    else iAverSpeed++;
    if(statusInj1>1) iGoodInj++;
    else if(statusInj1==1) iRescaledSpeed++;
    else iAverSpeed++;

    if(optVerbose) printf(" Mod. %d \tStatusLR=%X %X \t TimeStamp=%d \t v(an 128l)= %f",iMod,statusInj0,statusInj1,vdrift0->GetEventTimestamp(),vdriftarr0->GetDriftSpeed(0,128));
    if(optVerbose) printf("        \t v(an 128r)= %f  Degree=%d %d\n",vdriftarr1->GetDriftSpeed(0,128),vdrift0->GetDegreeofPoly(),vdrift1->GetDegreeofPoly());

    Int_t n7=(statusInj0&(0x1F<<25))>>25;
    Int_t n6=(statusInj0&(0x1F<<20))>>20;
    Int_t n5=(statusInj0&(0x1F<<15))>>15;
    Int_t n4=(statusInj0&(0x1F<<5))>>10;
    Int_t n3=(statusInj0&(0x1F<<5))>>5;
    Int_t n2=statusInj0&0x1F;
    Float_t aveStatus0=(7.*n7+6.*n6+5.*n5+4.*n4+3.*n3+2.*n2)/32.;
    n7=(statusInj1&(0x1F<<25))>>25;
    n6=(statusInj1&(0x1F<<20))>>20;
    n5=(statusInj1&(0x1F<<15))>>15;
    n4=(statusInj1&(0x1F<<5))>>10;
    n3=(statusInj1&(0x1F<<5))>>5;
    n2=statusInj1&0x1F;
    Float_t aveStatus1=(7.*n7+6.*n6+5.*n5+4.*n4+3.*n3+2.*n2)/32.;

    Int_t index=1+(det-1)*2;
    if(lay==3){ 
      hlay3->SetBinContent(index,lad,aveStatus0);
      hlay3->SetBinContent(index+1,lad,aveStatus1);
    }
    if(lay==4){ 
      hlay4->SetBinContent(index,lad,aveStatus0);
      hlay4->SetBinContent(index+1,lad,aveStatus1);
    }



    Float_t vel0=0;
    Float_t pd0=0;
    if(vdrift0){ 
      vel0=vdrift0->GetDriftSpeedAtAnode(128);
      pd0=vdrift0->GetDegreeofPoly();
    }
    Float_t vel1=0;
    Float_t pd1=0;
    if(vdrift1){ 
      vel1=vdrift1->GetDriftSpeedAtAnode(128);
      pd1=vdrift1->GetDegreeofPoly();
    }
    Float_t Edrift=(1800-45)/291/0.012;  
    Float_t mob0=vel0*1.E5/Edrift;  
    Float_t temper0=293.15*TMath::Power((mob0/1350.),-1/2.4); 
    Float_t mob1=vel1*1.E5/Edrift;  
    Float_t temper1=293.15*TMath::Power((mob1/1350.),-1/2.4); 
    tempvsmod0->SetPoint(tempvsmod0->GetN(),(Float_t)iMod,temper0);
    tempvsmod1->SetPoint(tempvsmod1->GetN(),(Float_t)iMod,temper1);
    vvsmod0->SetPoint(vvsmod0->GetN(),(Float_t)iMod,vel0);
    vvsmod1->SetPoint(vvsmod1->GetN(),(Float_t)iMod,vel1);
    poldegvsmod0->SetPoint(poldegvsmod0->GetN(),(Float_t)iMod,pd0);
    poldegvsmod1->SetPoint(poldegvsmod1->GetN(),(Float_t)iMod,pd1);

    for(Int_t ipar=0; ipar<=vdrift0->GetDegreeofPoly(); ipar++){
      fPoly->SetParameter(ipar,vdrift0->GetDriftSpeedParameter(ipar));
    }
    if(vdrift0->GetDegreeofPoly()<3){
      for(Int_t ipar=vdrift0->GetDegreeofPoly()+1; ipar<=3; ipar++) fPoly->SetParameter(ipar,0.);
    }

    anmaxvsmod0->SetPoint(anmaxvsmod0->GetN(),(Float_t)iMod,fPoly->GetMaximumX(0.,256.));
    dvcevsmod0->SetPoint(dvcevsmod0->GetN(),(Float_t)iMod,fPoly->Eval(128)-fPoly->Eval(0));
    dveevsmod0->SetPoint(dveevsmod0->GetN(),(Float_t)iMod,fPoly->Eval(256)-fPoly->Eval(0));
    
    for(Int_t ipar=0; ipar<=vdrift1->GetDegreeofPoly(); ipar++){
      fPoly->SetParameter(ipar,vdrift1->GetDriftSpeedParameter(ipar));
    }
    if(vdrift1->GetDegreeofPoly()<3){
      for(Int_t ipar=vdrift1->GetDegreeofPoly()+1; ipar<=3; ipar++) fPoly->SetParameter(ipar,0.);
    }
    anmaxvsmod1->SetPoint(anmaxvsmod1->GetN(),(Float_t)iMod,fPoly->GetMaximumX(0.,256.));
    dvcevsmod1->SetPoint(dvcevsmod1->GetN(),(Float_t)iMod,fPoly->Eval(128)-fPoly->Eval(0));
    dveevsmod1->SetPoint(dveevsmod1->GetN(),(Float_t)iMod,fPoly->Eval(256)-fPoly->Eval(0));
    //    getchar();
  }

  printf("Number of half-modules with drift speed from injectors               = %d\n",iGoodInj);
  printf("Number of half-modules with drift speed rescaled from golden module = %d\n",iRescaledSpeed);
  printf("Number of half-modules with average drift speed                      = %d\n",iAverSpeed);

  gStyle->SetPalette(59);

  TCanvas* cinjst=new TCanvas("cinjst","Injector Status",900,600);
  cinjst->Divide(2,1);
  cinjst->cd(1);
  hlay3->Draw("colz");
  TLine** linv3=new TLine*[5];
  for(Int_t i=0;i<5;i++){
    linv3[i]=new TLine(i+0.5,-0.5,i+0.5,13.5);
    linv3[i]->SetLineColor(kGray+1);
    linv3[i]->Draw();
  }
  TLine** linh3=new TLine*[13];
  for(Int_t i=0;i<13;i++){
    linh3[i]=new TLine(-0.5,i+0.5,5.5,i+0.5);
    linh3[i]->SetLineColor(kGray+1);
    linh3[i]->Draw();
  }
  cinjst->cd(2);
  hlay4->Draw("colz");
  TLine** linv4=new TLine*[7];
  for(Int_t i=0;i<7;i++){
    linv4[i]=new TLine(i+0.5,-0.5,i+0.5,21.5);
    linv4[i]->SetLineColor(kGray+1);
    linv4[i]->Draw();
  }
  TLine** linh4=new TLine*[21];
  for(Int_t i=0;i<21;i++){
    linh4[i]=new TLine(-0.5,i+0.5,7.5,i+0.5);
    linh4[i]->SetLineColor(kGray+1);
    linh4[i]->Draw();
  }
  cinjst->Modified();


  TCanvas* c2;
  c2=new TCanvas("c2","Vdrift vs. mod",1000,700);
  vvsmod0->SetMarkerStyle(20);
  vvsmod0->Draw("AP");
  vvsmod0->GetXaxis()->SetTitle("Module Number");
  vvsmod0->GetYaxis()->SetTitle("Vdrift (#mum/ns)");
  vvsmod1->SetMarkerStyle(21);
  vvsmod1->SetMarkerColor(2);
  vvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();

  TCanvas* c2t;
  c2t=new TCanvas("c2t","Temper vs. mod",1000,700);
  tempvsmod0->SetMarkerStyle(20);
  tempvsmod0->Draw("AP");
  tempvsmod0->GetXaxis()->SetTitle("Module Number");
  tempvsmod0->GetYaxis()->SetTitle("Temperature (K)");
  tempvsmod1->SetMarkerStyle(21);
  tempvsmod1->SetMarkerColor(2);
  tempvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();

  TCanvas* c3;
  c3=new TCanvas("c3","Params vs. mod",900,900);
  c3->Divide(2,2);
  
  c3->cd(1);
  gPad->SetLeftMargin(0.14);
  poldegvsmod0->SetMarkerStyle(20);
  poldegvsmod0->Draw("AP");
  poldegvsmod0->GetXaxis()->SetTitle("Module Number");
  poldegvsmod0->GetYaxis()->SetTitle("Degree of Polynomial fit");
  poldegvsmod0->GetYaxis()->SetTitleOffset(1.4);
  poldegvsmod1->SetMarkerStyle(21);
  poldegvsmod1->SetMarkerColor(2);
  poldegvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  c3->cd(2);
  gPad->SetLeftMargin(0.14);
  anmaxvsmod0->SetMarkerStyle(20);
  anmaxvsmod0->Draw("AP");
  anmaxvsmod0->GetXaxis()->SetTitle("Module Number");
  anmaxvsmod0->GetYaxis()->SetTitle("Anode with max. drift speed");
  anmaxvsmod0->GetYaxis()->SetTitleOffset(1.4);
  anmaxvsmod1->SetMarkerStyle(21);
  anmaxvsmod1->SetMarkerColor(2);
  anmaxvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  c3->cd(3);
  gPad->SetLeftMargin(0.14);
  dvcevsmod0->SetMarkerStyle(20);
  dvcevsmod0->Draw("AP");
  dvcevsmod0->GetXaxis()->SetTitle("Module Number");
  dvcevsmod0->GetYaxis()->SetTitle("vdrift(anode128)-vdrift(anode0)");
  dvcevsmod0->GetYaxis()->SetTitleOffset(1.4);
  dvcevsmod1->SetMarkerStyle(21);
  dvcevsmod1->SetMarkerColor(2);
  dvcevsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  c3->cd(4);
  gPad->SetLeftMargin(0.14);
  dveevsmod0->SetMarkerStyle(20);
  dveevsmod0->Draw("AP");
  dveevsmod0->GetYaxis()->SetTitleOffset(1.4);
  dveevsmod0->GetXaxis()->SetTitle("Module Number");
  dveevsmod0->GetYaxis()->SetTitle("vdrift(anode256)-vdrift(anode0)");
  dveevsmod1->SetMarkerStyle(21);
  dveevsmod1->SetMarkerColor(2);
  dveevsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();

  
}


void PlotOfflineCalib(AliCDBEntry *resp){
  printf("====== PARAMETERS FROM OFFLINE CALIBRATION ======\n");

  AliITSresponseSDD* r=(AliITSresponseSDD*)resp->GetObject();
  
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

