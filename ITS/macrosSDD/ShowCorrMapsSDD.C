#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TGrid.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include "AliCDBEntry.h"
#include "AliITSCorrMapSDD.h"
#include "AliITSMapSDD.h"
#endif


// Macro to plot the SDD correction maps 
// from the OCDB file (OCDB/ITS/Calib/MapsTimeSDD)
//
// Origin: F. Prino (prino@to.infn.it)

void ShowCorrMapsSDD(TString filname="alien:///alice/data/2010/OCDB/ITS/Calib/MapsTimeSDD/Run117224_999999999_v3_s0.root"){

  if(filname.Contains("alien")){
    TGrid::Connect("alien:");
  }
  TFile *fil =TFile::Open(filname.Data());
    
  AliCDBEntry* ent=(AliCDBEntry*)fil->Get("AliCDBEntry");
  TObjArray *maptSDD = (TObjArray *)ent->GetObject();
  printf("Entries in array=%d\n",maptSDD->GetEntriesFast());
  TString psnm0 = "mapsSDD.ps[";
  TString psnm1 = "mapsSDD.ps";
  TString psnm2 = "mapsSDD.ps]";
  Bool_t oldMapFormat=kFALSE;
  TObject* objmap=(TObject*)maptSDD->At(0);
  TString cname(objmap->ClassName());
  if(cname.CompareTo("AliITSMapSDD")==0){ 
    oldMapFormat=kTRUE;
    printf("SDD Maps converted to new format\n");
  }
 

  TCanvas * c3=new TCanvas("c3","Layer 3",1200,900);
  c3->Print(psnm0.Data());
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  Int_t cntpad=0;
  for(Int_t i=0; i<84; i++){
    if(i%6==0){
      c3->cd();
      c3->Modified();
      c3->Update();
      if(i) c3->Print(psnm1.Data());
      c3->Clear();
      c3->Divide(4,3);
      cntpad=0;
    }
    Int_t index0=i*2;
    Int_t index1=i*2+1;
    AliITSCorrMapSDD* map0 = 0;
    AliITSCorrMapSDD* map1 = 0;
    if(oldMapFormat){ 
      AliITSMapSDD* oldmap0=(AliITSMapSDD*)maptSDD->At(index0);
      AliITSMapSDD* oldmap1=(AliITSMapSDD*)maptSDD->At(index1);
      map0=oldmap0->ConvertToNewFormat();
      map1=oldmap1->ConvertToNewFormat();
    }else{
      map0=(AliITSCorrMapSDD*)maptSDD->At(index0);
      map1=(AliITSCorrMapSDD*)maptSDD->At(index1);
    }

    printf("Module %s Entries in map %dx%d\n",map0->GetName(),map0->GetNBinsAnode(),map0->GetNBinsDrift());
    if(map0->GetNBinsAnode()==1){
      TH1F* hp0=map0->GetMapProfile();
      TH1F* hp1=map1->GetMapProfile();
      if(hp0->GetMinimum()>-0.000001) hp0->SetMinimum(-10);
      if(hp0->GetMaximum()<0.000001) hp0->SetMaximum(10);
      if(hp1->GetMinimum()>-0.000001) hp1->SetMinimum(-10);
      if(hp1->GetMaximum()<0.000001) hp1->SetMaximum(10);
      hp0->SetTitle(Form("Module %d - Left",i+240));
      hp1->SetTitle(Form("Module %d - Right",i+240));
      hp0->GetXaxis()->SetTitle("Drift distance (mm)");
      hp0->GetYaxis()->SetTitle("Correction (#mum)");
      hp1->GetXaxis()->SetTitle("Drift distance (mm)");
      hp1->GetYaxis()->SetTitle("Correction (#mum)");
      c3->cd(++cntpad);
      hp0->Draw();
      c3->cd(++cntpad);
      hp1->Draw();
    }else{
      TH2F* hmap0=map0->GetMapHisto();
      TH2F* hmap1=map1->GetMapHisto();
      c3->cd(++cntpad);
      hmap0->Draw("colz");
      c3->cd(++cntpad);
      hmap1->Draw("colz");
    }
  }
  c3->cd();
  c3->Modified();
  c3->Update();
  c3->Print(psnm1.Data());

  TCanvas * c4=new TCanvas("c4","Layer 4",1200,900);
  c4->Divide(4,4);
  gStyle->SetPalette(1);
  for(Int_t i=84; i<260; i++){
    Int_t j=i-84;
    if(j%8==0){
      c4->cd();
      c4->Modified();
      c4->Update();
      if(j) c4->Print(psnm1.Data());
      c4->Clear();
      c4->Divide(4,4);
      cntpad=0;
    }
    Int_t index0=i*2;
    Int_t index1=i*2+1;
    AliITSCorrMapSDD* map0 = 0;
    AliITSCorrMapSDD* map1 = 0;
    if(oldMapFormat){ 
      AliITSMapSDD* oldmap0=(AliITSMapSDD*)maptSDD->At(index0);
      AliITSMapSDD* oldmap1=(AliITSMapSDD*)maptSDD->At(index1);
      map0=oldmap0->ConvertToNewFormat();
      map1=oldmap1->ConvertToNewFormat();
    }else{
      map0=(AliITSCorrMapSDD*)maptSDD->At(index0);
      map1=(AliITSCorrMapSDD*)maptSDD->At(index1);
    }
    printf("Module %d Entries in map %dx%d\n",i,map0->GetNBinsAnode(),map0->GetNBinsDrift());
    if(map0->GetNBinsAnode()==1){
      TH1F* hp0=map0->GetMapProfile();
      TH1F* hp1=map1->GetMapProfile();
      if(hp0->GetMinimum()>-0.000001) hp0->SetMinimum(-10);
      if(hp0->GetMaximum()<0.000001) hp0->SetMaximum(10);
      if(hp1->GetMinimum()>-0.000001) hp1->SetMinimum(-10);
      if(hp1->GetMaximum()<0.000001) hp1->SetMaximum(10);
      hp0->SetTitle(Form("Module %d - Left",i+240));
      hp1->SetTitle(Form("Module %d - Right",i+240));
      hp0->GetXaxis()->SetTitle("Drift distance (mm)");
      hp0->GetYaxis()->SetTitle("Correction (#mum)");
      hp1->GetXaxis()->SetTitle("Drift distance (mm)");
      hp1->GetYaxis()->SetTitle("Correction (#mum)");
      c4->cd(++cntpad);
      hp0->Draw();
      c4->cd(++cntpad);
      hp1->Draw();
    }else{
      TH2F* hmap0=map0->GetMapHisto();
      TH2F* hmap1=map1->GetMapHisto();
      c4->cd(++cntpad);
      hmap0->Draw("colz");
      c4->cd(++cntpad);
      hmap1->Draw("colz");
    }
  }
  c4->cd();
  c4->Modified();
  c4->Update();
  c4->Print(psnm1.Data());
  c4->Print(psnm2.Data());

}
