#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH2F.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGrid.h>
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSDD.h"
#endif

// Macro for a z-phi event display of the SDD Raw Data
// Origin: F. Prino,   prino@to.infn.it

void DisplaySDDRawData(TString filename, Int_t firstEv=0, Int_t lastEv=5){


  Bool_t writtenoutput=kFALSE;
  AliITSDDLModuleMapSDD* ddlmap=new AliITSDDLModuleMapSDD();
  ddlmap->SetJun09Map();

  TH2F* hzphi3=new TH2F("hzphi3","Layer 3",1536,-0.5,1535.5,3584,-0.5,3584.5);
  TH2F* hzphi4=new TH2F("hzphi4","Layer 4",2048,-0.5,2047.5,5632,-0.5,5631.5);

  TLine** lA3=new TLine*[5];
  for(Int_t ilin=0;ilin<5;ilin++){
    lA3[ilin]=new TLine((ilin+1)*256,0,(ilin+1)*256,3584.5);
    lA3[ilin]->SetLineColor(kGray);
    lA3[ilin]->SetLineStyle(2);
  }
  TLine** lT3=new TLine*[13];
  for(Int_t ilin=0;ilin<13;ilin++){
    lT3[ilin]=new TLine(0,(ilin+1)*256,1535.5,(ilin+1)*256);
    lT3[ilin]->SetLineColor(kGray);
    lT3[ilin]->SetLineStyle(2);
  }

  TLine** lA4=new TLine*[7];
  for(Int_t ilin=0;ilin<7;ilin++){
    lA4[ilin]=new TLine((ilin+1)*256,0,(ilin+1)*256,5631.5);
    lA4[ilin]->SetLineColor(kGray);
    lA4[ilin]->SetLineStyle(2);
  }
  TLine** lT4=new TLine*[21];
  for(Int_t ilin=0;ilin<21;ilin++){
    lT4[ilin]=new TLine(0,(ilin+1)*256,2047.5,(ilin+1)*256);
    lT4[ilin]->SetLineColor(kGray);
    lT4[ilin]->SetLineStyle(2);
  }

  hzphi3->SetStats(0);
  hzphi4->SetStats(0);

  Int_t iev=firstEv;
  AliRawReader *rd; 
  if(filename.Contains(".root")){
    rd=new AliRawReaderRoot(filename.Data(),iev);
  }else{
    rd=new AliRawReaderDate(filename.Data(),iev);
  }

  TStopwatch *evtime=new TStopwatch();
  TCanvas* c0 = new TCanvas("cd0","c0",800,800);
  gStyle->SetPalette(1);
  do{
    c0->Clear();				
    c0->Divide(1,2,0.001,0.001);

    evtime->Start();
    printf("Event # %d\n",iev);
    rd->Reset();
    hzphi3->Reset();
    hzphi4->Reset();

    UChar_t cdhAttr=AliITSRawStreamSDD::ReadBlockAttributes(rd);
    UInt_t amSamplFreq=AliITSRawStreamSDD::ReadAMSamplFreqFromCDH(cdhAttr);
    AliITSRawStream* s=AliITSRawStreamSDD::CreateRawStreamSDD(rd,cdhAttr);
    if(!writtenoutput){
      printf("Use %s raw stream, sampling frequency %d MHz\n",s->ClassName(),amSamplFreq);
      writtenoutput=kTRUE;

    }

    while(s->Next()){
      
      if(s->IsCompletedModule()==kFALSE && s->IsCompletedDDL()==kFALSE){
	Int_t lay,lad,det;
	Int_t modID=ddlmap->GetModuleNumber(rd->GetDDLID(),s->GetCarlosId());
	AliITSgeomTGeo::GetModuleId(modID,lay,lad,det);
	Int_t iz=s->GetCoord1()+256*(det-1);
	Int_t iphi=s->GetCoord2()+256*(lad-1)+128*s->GetChannel();
	if(lay==3){
	  hzphi3->SetBinContent(iz+1,iphi+1,s->GetSignal());
	}else if(lay==4){
	  hzphi4->SetBinContent(iz+1,iphi+1,s->GetSignal());
	}
      }
    }
    evtime->Stop();
    printf("**** Event=%d \n",iev);
    evtime->Print("u");
    evtime->Reset();
    iev++;
    
    c0->cd(1);
    hzphi3->Draw("colz");
    for(Int_t ilin=0;ilin<5;ilin++) lA3[ilin]->Draw("same");
    for(Int_t ilin=0;ilin<13;ilin++) lT3[ilin]->Draw("same");
    hzphi3->GetXaxis()->SetTitle("Z (anode)");
    hzphi3->GetYaxis()->SetTitle("PHI (time bin)");
      
    c0->cd(2);
    hzphi4->Draw("colz");
    for(Int_t ilin=0;ilin<7;ilin++) lA4[ilin]->Draw("same");
    for(Int_t ilin=0;ilin<21;ilin++) lT4[ilin]->Draw("same");
    hzphi4->GetXaxis()->SetTitle("Z (anode)");
    hzphi4->GetYaxis()->SetTitle("PHI (time bin)");
    c0->Update();
  }while(rd->NextEvent()&&iev<=lastEv);

}

void DisplaySDDRawData(Int_t nrun, Int_t n2, Int_t chunk=10, Int_t year=2009, Char_t* dir="LHC09b", 
		       Int_t firstEv=21*3, 
		       Int_t lastEv=21*3+1){  
  TGrid::Connect("alien:",0,0,"t");
  TString filnam(Form("alien:///alice/data/%d/%s/%09d/raw/%02d%09d%03d.%02d.root",year,dir,nrun,year-2000,nrun,n2,chunk));
  printf("Open file %s\n",filnam.Data());
  DisplaySDDRawData(filnam,firstEv,lastEv);
}

