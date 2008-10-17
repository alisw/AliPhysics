#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH2F.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSDD.h"
#endif

// Macro to display the SDD Raw Data for 1 DDL
// Origin: F. Prino,   prino@to.infn.it

void PlotSDDRawData(Char_t datafil[100], Int_t nDDL, Int_t firstEv=0, Int_t lastEv=5){

  const Int_t nHybrids=24;

  TH2F** histo = new TH2F*[nHybrids];
  Char_t nome[20];
  for(Int_t i=0;i<nHybrids;i++){
    sprintf(nome,"histo%d",i);
    histo[i]=new TH2F(nome,"",256,-0.5,255.5,256,-0.5,255.5);
    histo[i]->SetStats(0);
  }

  Int_t iev=firstEv;
  AliRawReader *rd; 
  if(strstr(datafil,".root")!=0){
    rd=new AliRawReaderRoot(datafil,iev);
  }else{
    rd=new AliRawReaderDate(datafil,iev);
  }
  TStopwatch *evtime=new TStopwatch();
  TCanvas* c0 = new TCanvas("cd0","c0",900,900);
  gStyle->SetPalette(1);
  Int_t idev;
  do{
    c0->Clear();				
    c0->Divide(4,6,0.001,0.001);

    evtime->Start();
    printf("Event # %d\n",iev);
    rd->Reset();
    for(Int_t i=0;i<nHybrids;i++) histo[i]->Reset();
    AliITSRawStreamSDD s(rd);
    Int_t iCountNext=0;    
    while(s.Next()){
      iCountNext++;
      if(s.IsCompletedModule()==kFALSE && s.IsCompletedDDL()==kFALSE){
 	Int_t i=s.GetCarlosId()*2+s.GetChannel();
	if(rd->GetDDLID()==nDDL) histo[i]->Fill(s.GetCoord2(),s.GetCoord1(),s.GetSignal());
      }
    }
    idev=s.GetEventId();
    evtime->Stop();
    printf("**** Event=%d  ID=%d\n",iev,idev);
    evtime->Print("u");
    evtime->Reset();
    iev++;
    
    for(Int_t i=0;i<nHybrids;i++){
      c0->cd(i+1);
      histo[i]->DrawCopy("colz");
    }
    c0->Update();
  }while(rd->NextEvent()&&iev<=lastEv);

}

