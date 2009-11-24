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
#include <TLine.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include "AliCDBEntry.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSgeomTGeo.h"
#endif

// Macro to plot the calibration parameters from the OCDB file 
// created from PEDESTAL and PULSER runs (OCDB/ITS/Calib/CalibSDD)
// Two methods ShowCalibrationSDD:
//  - the first takes the name of the file to be displayed
//  - the second builds the alien path+name from run number and file version
//
// Origin: F. Prino (prino@to.infn.it)

void MakePalette(){
  Int_t palette[3]={kGray,2,3};  
  gStyle->SetPalette(3,palette);
}

void ShowCalibrationSDD(Char_t *filnam="$ALICE_ROOT/OCDB/ITS/Calib/CalibSDD/Run0_9999999_v0_s0.root", Int_t iMod=0){


  TFile *f=TFile::Open(filnam);
  AliCDBEntry *ent=(AliCDBEntry*)f->Get("AliCDBEntry");
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
  TH2I* hdeadlay3=new TH2I("hlay3","Layer 3",6,-0.5,5.5,14,-0.5,13.5);
  hdeadlay3->GetXaxis()->SetTitle("Detector");
  hdeadlay3->GetYaxis()->SetTitle("Ladder");
  hdeadlay3->GetXaxis()->SetTickLength(0);
  hdeadlay3->GetYaxis()->SetTickLength(0);
  hdeadlay3->SetStats(0);
  hdeadlay3->SetMinimum(-1.);
  TH2I* hdeadlay4=new TH2I("hlay4","Layer 4",8,-0.5,7.5,22,-0.5,21.5);
  hdeadlay4->GetXaxis()->SetTitle("Detector");
  hdeadlay4->GetYaxis()->SetTitle("Ladder");
  hdeadlay4->GetXaxis()->SetTickLength(0);
  hdeadlay4->GetYaxis()->SetTickLength(0);
  hdeadlay4->GetYaxis()->SetTitle("Ladder");
  hdeadlay4->SetStats(0);
  hdeadlay4->SetMinimum(-1.);

  TObjArray *calSDD = (TObjArray *)ent->GetObject();
  printf("Entries in array=%d\n",calSDD->GetEntriesFast());
  TH1F* hmodstatus=new TH1F("hmodstatus","",260,0.5,260.5);
  TH1F* hnbadch=new TH1F("hnbadch","",260,0.5,260.5);
  TH1F* hbase=new TH1F("hbase","",60,0.5,120.5);
  TH1F* hnoise=new TH1F("hnoise","",100,0.,7.);
  TH1F* hgain=new TH1F("hgain","",100,0.,4.);
  TH1F* hchstatus=new TH1F("hchstatus","",2,-0.5,1.5);
  AliITSCalibrationSDD *cal;
  Int_t badModCounter3=0;
  Int_t badModCounter4=0;
  Int_t badAnodeCounter3=0;
  Int_t badAnodeCounter4=0;
  Int_t badAnodeCounterGoodMod3=0;
  Int_t badAnodeCounterGoodMod4=0;
  Int_t badAnodeCounterGoodModAndChip3=0;
  Int_t badAnodeCounterGoodModAndChip4=0;
  Int_t badChipCounter3=0;
  Int_t badChipCounter4=0;
  for(Int_t i=0; i<260; i++){
    cal=(AliITSCalibrationSDD*)calSDD->At(i);
    if(cal==0) continue;
    printf("Module %d (%d)   status = ",i,i+240);
    Int_t lay,lad,det;
    AliITSgeomTGeo::GetModuleId(i+240,lay,lad,det);
    Int_t index=1+(det-1)*2;
    if(cal->IsBad()){ 
      printf("BAD\t");
      if(lay==3){ 
	badModCounter3++;
	hlay3->SetBinContent(index,lad,0);
	hlay3->SetBinContent(index+1,lad,0);
      }else if(lay==4){ 
	badModCounter4++;
	hlay4->SetBinContent(index,lad,0);
	hlay4->SetBinContent(index+1,lad,0);
      }
      hmodstatus->SetBinContent(i+1,0);
    }else{ 
      printf("OK\t");
      hmodstatus->SetBinContent(i+1,1);
      if(lay==3){ 
	badAnodeCounterGoodMod3+=cal->GetDeadChannels();
	if(cal->IsChipBad(0) && cal->IsChipBad(1) && cal->IsChipBad(2) && cal->IsChipBad(3)){
	  hlay3->SetBinContent(index,lad,0);
	}else{
	  hlay3->SetBinContent(index,lad,1);
	}
	if(cal->IsChipBad(4) && cal->IsChipBad(5) && cal->IsChipBad(6) && cal->IsChipBad(7)){
	  hlay3->SetBinContent(index+1,lad,0);
	}else{
	  hlay3->SetBinContent(index+1,lad,1);
	}
      }else{ 
	badAnodeCounterGoodMod4+=cal->GetDeadChannels();
	if(cal->IsChipBad(0) && cal->IsChipBad(1) && cal->IsChipBad(2) && cal->IsChipBad(3)){
	  hlay4->SetBinContent(index,lad,0);
	}else{
	  hlay4->SetBinContent(index,lad,1);
	}
	if(cal->IsChipBad(4) && cal->IsChipBad(5) && cal->IsChipBad(6) && cal->IsChipBad(7)){
	  hlay4->SetBinContent(index+1,lad,0);
	}else{
	  hlay4->SetBinContent(index+1,lad,1);
	}
      }
     }
    printf("   Chip Status (0=OK, 1=BAD): ");  
    for(Int_t ic=0; ic<8;ic++){ 
      printf("%d ",cal->IsChipBad(ic));
      if(cal->IsChipBad(ic) && !cal->IsBad()){ 
	if(i<84) badChipCounter3++;
	else badChipCounter4++;
      }
    }
    printf(" # bad anodes = %d  ",cal->GetDeadChannels());
    if(cal->IsAMAt20MHz()) printf("      20 MHz sampling");
    else printf("      40 MHz sampling");
    printf("\n");
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
	hchstatus->Fill(1);
	hnoise->Fill(noise);
	hgain->Fill(gain);
      }
    }
  }
  Int_t totbad3=badModCounter3*512+badChipCounter3*64+badAnodeCounterGoodModAndChip3;
  Int_t tot3=6*14*512;
  Float_t fracbad3=(Float_t)totbad3/(Float_t)tot3;
  Int_t totbad4=badModCounter4*512+badChipCounter4*64+badAnodeCounterGoodModAndChip4;
  Int_t tot4=8*22*512;
  Float_t fracbad4=(Float_t)totbad4/(Float_t)tot4;
  Float_t fractot=(Float_t)(totbad3+totbad4)/(Float_t)(tot3+tot4);
  printf("----------------------Summary----------------------\n");
  printf("---- Layer 3 ----\n");
  printf("# of bad modules                      = %d\n",badModCounter3);
  printf("# of bad chips in good modules        = %d\n",badChipCounter3);
  printf("# of bad anodes in good modules+chips = %d\n",badAnodeCounterGoodModAndChip3);
  printf("Fraction of bads (anodes+chips+mod)   = %f\n",fracbad3);
  printf("---- Layer 4 ----\n");
  printf("# of bad modules                      = %d\n",badModCounter4);
  printf("# of bad chips in good modules        = %d\n",badChipCounter4);
  printf("# of bad anodes in good modules+chips = %d\n",badAnodeCounterGoodModAndChip4);
  printf("Fraction of bads (anodes+chips+mod)   = %f\n",fracbad4);
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

  TCanvas *c0=new TCanvas("c0","Module status",800,800);
  c0->Divide(1,2);
  c0->cd(1);
  hmodstatus->Draw();
  hmodstatus->GetXaxis()->SetTitle("Module number");
  hmodstatus->GetYaxis()->SetTitle("Module status (1=OK, 0=BAD)");
  c0->cd(2);
  hnbadch->Draw();
  hnbadch->GetXaxis()->SetTitle("Module number");   
  hnbadch->GetYaxis()->SetTitle("Number of bad anodes");




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










  // Plot quantities for specified module

  cal=(AliITSCalibrationSDD*)calSDD->At(iMod);
  if(cal==0) return;
  printf("-----------------------------------\n");
  printf("Module %d    status = ",iMod);
  if(cal->IsBad()) printf("BAD\n");
  else printf("OK\n");
  printf("   Chip Status (0=OK, 1=BAD): ");  
  for(Int_t ic=0; ic<8;ic++) printf("%d ",cal->IsChipBad(ic));
  printf("\n");
  printf("   Number of bad anodes =%d\n",cal->GetDeadChannels());
  printf("-----------------------------------\n");
  Int_t ipt=0;
  TGraph *gbad=new TGraph(0);
  gbad->SetTitle("Bad Channels");
  TGraph *gbase=new TGraph(0);
  gbase->SetTitle("Baselines");
  TGraph *gnoi=new TGraph(0);
  gnoi->SetTitle("Noise");
  TGraph *ggain=new TGraph(0);
  ggain->SetTitle("Gain");
  for(Int_t iAn=0; iAn<512; iAn++){
    Float_t bad=1;
    if(cal->IsBadChannel(iAn)) bad=0;
    Float_t base=cal->GetBaseline(iAn);
    Float_t noise=cal->GetNoiseAfterElectronics(iAn);
    Float_t gain=cal->GetChannelGain(iAn);
    gbad->SetPoint(ipt,(Float_t)iAn,bad);
    gbase->SetPoint(ipt,(Float_t)iAn,base);
    ggain->SetPoint(ipt,(Float_t)iAn,gain);
    gnoi->SetPoint(ipt,(Float_t)iAn,noise);
    ipt++;
  }
  Char_t ctit[100];
  sprintf(ctit,"Module %d",iMod);

  TCanvas *c2=new TCanvas("c2",ctit,1200,800);
  c2->Divide(2,2);
  
  c2->cd(1);
  gbase->SetMarkerStyle(7);
  gbase->Draw("AP");
  gbase->GetXaxis()->SetTitle("Anode Number");
  gbase->GetYaxis()->SetTitle("Baseline after equalization");  
  c2->cd(2);
  gnoi->SetMarkerStyle(7);
  gnoi->Draw("AP");
  gnoi->GetXaxis()->SetTitle("Anode Number");
  gnoi->GetYaxis()->SetTitle("Noise");  
  c2->cd(3);
  ggain->SetMarkerStyle(7);
  ggain->Draw("AP");
  ggain->GetXaxis()->SetTitle("Anode Number");
  ggain->GetYaxis()->SetTitle("Gain");
  c2->cd(4);
  gbad->SetMarkerStyle(7);
  gbad->Draw("AP");
  gbad->SetMinimum(-0.1);
  gbad->GetXaxis()->SetTitle("Anode Number");
  gbad->GetYaxis()->SetTitle("Anode Status (1=OK, 0=bad)");
}

void ShowCalibrationSDD(Int_t nrun, Int_t year=2009, Int_t nmod=0){
  TGrid::Connect("alien:",0,0,"t");
  TString cmd=Form("gbbox find \"/alice/data/%d/OCDB/ITS/Calib/CalibSDD\" \"Run%d*.root\" > run.txt",year,nrun);
  gSystem->Exec(cmd.Data());
  Char_t filnam[200],filnamalien[200];
  FILE* runtxt=fopen("run.txt","r");
  fscanf(runtxt,"%s\n",filnam);    
  if(!strstr(filnam,"/alice/data/")){
    printf("Bad run number\n");
    gSystem->Exec("rm run.txt");
    return;
  }  
  sprintf(filnamalien,"alien://%s",filnam);
  
  printf("Open file: %s\n",filnamalien);
  ShowCalibrationSDD(filnamalien,nmod);
  fclose(runtxt);
  gSystem->Exec("rm run.txt");
}
