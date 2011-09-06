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

void PlotSDDRawData(Char_t datafil[100], 
		    Int_t nDDL, 
		    Int_t firstEv=18, 
		    Int_t lastEv=20){

  // Main function

  const Int_t nHybrids=24;
  Bool_t writtenoutput=kFALSE;
  TH2F** histo = new TH2F*[nHybrids];
  for(Int_t i=0;i<nHybrids;i++){
    histo[i]=new TH2F(Form("histo%d",i),"",256,-0.5,255.5,256,-0.5,255.5);
    histo[i]->SetStats(0);
  }

  Bool_t isSingleMod=kFALSE;
  Int_t npx=4;
  Int_t npy=6;
  Int_t xsiz=700;
  Int_t ysiz=700;
  Int_t nHybrToPlot=24;
  Int_t iMod=-1;
  Int_t nCarlos;
  if(nDDL>=240 && nDDL<500){
    iMod=nDDL;
    AliITSDDLModuleMapSDD* dmap=new AliITSDDLModuleMapSDD();
    dmap->SetJun09Map();
    dmap->FindInDDLMap(iMod,nDDL,nCarlos);
    histo[nCarlos*2]->SetTitle(Form("Module %d Side 0",iMod));
    histo[nCarlos*2+1]->SetTitle(Form("Module %d Side 1",iMod));
    isSingleMod=kTRUE;
    npx=2;
    npy=1;
    xsiz=900;
    ysiz=450;
    nHybrToPlot=2;
  }

  Int_t iev=firstEv;
  AliRawReader *rd; 
  if(strstr(datafil,".root")!=0){
    rd=new AliRawReaderRoot(datafil,iev);
  }else{
    rd=new AliRawReaderDate(datafil,iev);
  }
  TStopwatch *evtime=new TStopwatch();
  TCanvas* c0 = new TCanvas("cd0","c0",xsiz,ysiz);
  gStyle->SetPalette(1);

  do{
    c0->Clear();				
    c0->Divide(npx,npy,0.001,0.001);

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
    
    for(Int_t i=0;i<nHybrToPlot;i++){
      c0->cd(i+1);
      if(isSingleMod){
	histo[nCarlos*2+i]->DrawCopy("colz");
      }else{
	histo[i]->DrawCopy("colz");
      }
    }
    c0->Update();
    //    if(histo[nCarlos*2]->GetMaximum()>1) getchar();
  }while(rd->NextEvent()&&iev<=lastEv);

}

void PlotSDDRawData(Int_t nrun, 
		    Int_t n2, 
		    Int_t year=2011, 
		    Char_t* dir="LHC11d_SDD",
		    Int_t nDDL=0, 
		    Int_t firstEv=18, 
		    Int_t lastEv=20){

  // Get file directly from alien

  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/%d/%s/%09d/raw/%02d%09d%03d.10.root",year,dir,nrun,year-2000,nrun,n2);
  printf("Open file %s\n",filnam);
  PlotSDDRawData(filnam,nDDL,firstEv,lastEv);
}

void PlotSDDRawData(Char_t datafil[100], 
		    Int_t nLay, 
		    Int_t nLad, 
		    Int_t nDet,
		    Int_t firstEv, 
		    Int_t lastEv){

  // plot raw data for single module starting from 
  // Layer, Ladder and detector numbers (counted from 1)
  Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(nLay,nLad,nDet);
  PlotSDDRawData(datafil,modIndex,firstEv,lastEv);

}
