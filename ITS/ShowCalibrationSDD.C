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

void ShowCalibrationSDD(Char_t *filnam="$ALICE_ROOT/ITS/Calib/CalibSDD/Run0_9999999_v0_s0.root", Int_t iMod=0){


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
  for(Int_t i=0; i<260; i++){
    cal=(AliITSCalibrationSDD*)calSDD->At(i);
    if(cal==0) continue;
    printf("Module %d (%d)   status = ",i,i+240);
    if(cal->IsBad()) printf("BAD\t");
    else printf("OK\t");
    printf("   Chip Status (0=OK, 1=BAD): ");  
    for(Int_t ic=0; ic<8;ic++) printf("%d ",cal->IsChipBad(ic));
    if(cal->IsAMAt20MHz()) printf("      20 MHz sampling");
    else printf("      40 MHz sampling");
    printf("\n");
    if(cal->IsBad()) hmodstatus->SetBinContent(i+1,0);
    else hmodstatus->SetBinContent(i+1,1);
    hnbadch->SetBinContent(i+1,cal->GetDeadChannels());
    for(Int_t iAn=0; iAn<512; iAn++){
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

  cal=(AliITSCalibrationSDD*)calSDD->At(iMod);
  if(cal==0) return;
  printf("-----------------------------------\n");
  printf("Module %d    status = ",iMod);
  if(cal->IsBad()) printf("BAD\n");
  else printf("OK\n");
  printf("   Chip Status (0=OK, 1=BAD): ");  
  for(Int_t ic=0; ic<8;ic++) printf("%d ",cal->IsChipBad(ic));
  printf("\n");
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

void ShowCalibrationSDD(Int_t nrun, Int_t nv,Int_t nmod=0){
  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/2008/LHC08c/OCDB/ITS/Calib/CalibSDD/Run%d_999999999_v%d_s0.root",nrun,nv);
  printf("Open file: %s\n",filnam);
  ShowCalibrationSDD(filnam,nmod);
}
