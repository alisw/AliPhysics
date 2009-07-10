#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH2F.h>
#include <TGrid.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSDDCompressed.h"
#endif

// Macro to display the SDD Raw Data for 1 DDL
// Origin: F. Prino,   prino@to.infn.it

void PlotSDDRawData(Char_t datafil[100], Int_t nDDL, Int_t firstEv=18, Int_t lastEv=20){

  const Int_t nHybrids=24;
  Bool_t writtenoutput=kFALSE;
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

  do{
    c0->Clear();				
    c0->Divide(4,6,0.001,0.001);

    evtime->Start();
    printf("Event # %d\n",iev);
    rd->Reset();
    for(Int_t i=0;i<nHybrids;i++) histo[i]->Reset();
    UChar_t cdhAttr=AliITSRawStreamSDD::ReadBlockAttributes(rd);
    UInt_t amSamplFreq=AliITSRawStreamSDD::ReadAMSamplFreqFromCDH(cdhAttr);
    AliITSRawStream* s=AliITSRawStreamSDD::CreateRawStreamSDD(rd,cdhAttr);
    if(!writtenoutput){
      printf("Use %s raw stream, sampling frequency %d MHz\n",s->ClassName(),amSamplFreq);
      writtenoutput=kTRUE;
    }

    Int_t iCountNext=0;    
    while(s->Next()){
      iCountNext++;
      if(s->IsCompletedModule()==kFALSE && s->IsCompletedDDL()==kFALSE){
 	Int_t i=s->GetCarlosId()*2+s->GetChannel();
	if(rd->GetDDLID()==nDDL) histo[i]->Fill(s->GetCoord2(),s->GetCoord1(),s->GetSignal());
      }
    }
    evtime->Stop();
    printf("**** Event=%d  \n",iev);
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

void PlotSDDRawData(Int_t nrun, Int_t n2, Int_t year=2009, Char_t* dir="LHC09b_SDD",
		    Int_t nDDL=0, 
		    Int_t firstEv=18, 
		    Int_t lastEv=20){

  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/%d/%s/%09d/raw/%02d%09d%03d.10.root",year,dir,nrun,year-2000,nrun,n2);
  printf("Open file %s\n",filnam);
  PlotSDDRawData(filnam,nDDL,firstEv,lastEv);
}

