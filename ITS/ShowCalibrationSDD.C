#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include "AliCDBEntry.h"
#include "AliITSCalibrationSDD.h"
#endif

// Macro to plot the calibration parameters from the OCDB file 
// created from PEDESTAL and PULSER runs (OCDB/ITS/Calib/CalibSDD)
// Two methods ShowCalibrationSDD:
//  - the first takes the name of the file to be displayed
//  - the second builds the alien path+name from run number and file version
//
// Origin: F. Prino (prino@to.infn.it)

void ShowCalibrationSDD(Int_t iMod=0, Char_t *filnam="$ALICE_ROOT/ITS/Calib/CalibSDD/Run0_9999999_v0_s0.root"){


  TFile *f=TFile::Open(filnam);
  AliCDBEntry *ent=(AliCDBEntry*)f->Get("AliCDBEntry");
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
    if(cal->IsBad()){ 
      printf("BAD\t");
      if(i<84) badModCounter3++;
      else badModCounter4++;
      hmodstatus->SetBinContent(i+1,0);
    }
    else{ 
      printf("OK\t");
      hmodstatus->SetBinContent(i+1,1);
      if(i<84) badAnodeCounterGoodMod3+=cal->GetDeadChannels();
      else badAnodeCounterGoodMod4+=cal->GetDeadChannels();
    }
    printf("   Chip Status (0=OK, 1=BAD): ");  
    for(Int_t ic=0; ic<8;ic++){ 
      printf("%d ",cal->IsChipBad(ic));
      if(cal->IsChipBad(ic) && !cal->IsBad()){ 
	if(i<84) badChipCounter3++;
	else badChipCounter4++;
      }
    }
    if(cal->IsAMAt20MHz()) printf("      20 MHz sampling");
    else printf("      40 MHz sampling");
    printf("\n");
    if(i<84) badAnodeCounter3+=cal->GetDeadChannels();
    else badAnodeCounter4+=cal->GetDeadChannels();
    hnbadch->SetBinContent(i+1,cal->GetDeadChannels());
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
      if(!cal->IsBadChannel(iAn)){
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

void ShowCalibrationSDD(Int_t nrun, Int_t nv, Char_t* dir="LHC08d", Int_t nmod=0){
  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/2008/%s/OCDB/ITS/Calib/CalibSDD/Run%d_999999999_v%d_s0.root",dir,nrun,nv);
  printf("Open file: %s\n",filnam);
  ShowCalibrationSDD(nmod,filnam);
}
