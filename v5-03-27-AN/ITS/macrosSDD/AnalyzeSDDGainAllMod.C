#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH2F.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TFile.h>
#include <TGrid.h>
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSOnlineSDDBase.h"
#include "AliITSOnlineSDDCMN.h"
#include "AliITSOnlineSDDTP.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSDDCompressed.h"
#endif

// Macro for the analysis of PULSER runs (equivalent to ITSSDDGAINda.cxx)
// Two functions named AnalyzeSDDGainAllModules: 
// The first is for analyzing a local raw data file and takes as agrument the file name.
// The second is for running on ALIEN
// All DDLs are analyzed, the argument nDDL selects the DDL to be plotted
// Origin: F. Prino (prino@to.infn.it)


void AnalyzeSDDGainAllMod(Char_t *datafil, 
			  Int_t adcfreq=20, 
			  Int_t nDDL=0, 
			  Int_t firstEv=18, 
			  Int_t lastEv=22, 
			  Float_t pascalDAC=100){

  const Int_t kTotDDL=24;
  const Int_t kModPerDDL=12;
  const Int_t kSides=2;
  Bool_t writtenoutput=kFALSE;


  TH2F** histo = new TH2F*[kTotDDL*kModPerDDL*kSides];
  AliITSOnlineSDDTP **anal=new AliITSOnlineSDDTP*[kTotDDL*kModPerDDL*kSides];
  Bool_t isFilled[kTotDDL*kModPerDDL*kSides];

  Char_t hisnam[20];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	anal[index]=new AliITSOnlineSDDTP(iddl,imod,isid,pascalDAC);
	if(adcfreq==40) anal[index]->SetLastGoodTB(254);
	else anal[index]->SetLastGoodTB(126);
	sprintf(hisnam,"h%02dc%02ds%d",iddl,imod,isid);
	histo[index]=new TH2F(hisnam,"",256,-0.5,255.5,256,-0.5,255.5);
	isFilled[index]=0;
      }
    }
  }

  TCanvas* c0 = new TCanvas("c0","Ev Display",900,900);
  gStyle->SetPalette(1);
  Char_t text[50];

  Int_t iev=firstEv;
  AliRawReader *rd; 
  if(strstr(datafil,".root")!=0){
    rd=new AliRawReaderRoot(datafil,iev);
  }else{
    rd=new AliRawReaderDate(datafil,iev);
  }
  TLatex *t0=new TLatex();
  t0->SetNDC();
  t0->SetTextSize(0.06);
  t0->SetTextColor(4);

  do{
    c0->Clear();
    c0->Divide(4,6,0.001,0.001);
    printf("Event # %d\n",iev);
    rd->Reset();
    for(Int_t iddl=0; iddl<kTotDDL;iddl++){
      for(Int_t imod=0; imod<kModPerDDL;imod++){
	for(Int_t isid=0;isid<kSides;isid++){
	  Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	  histo[index]->Reset();
	}
      }
    }

    UChar_t cdhAttr=AliITSRawStreamSDD::ReadBlockAttributes(rd);
    UInt_t amSamplFreq=AliITSRawStreamSDD::ReadAMSamplFreqFromCDH(cdhAttr);
    AliITSRawStream* s=AliITSRawStreamSDD::CreateRawStreamSDD(rd,cdhAttr);
    if(!writtenoutput){
      printf("Use %s raw stream, sampling frequency %d MHz\n",s->ClassName(),amSamplFreq);
      writtenoutput=kTRUE;
    }
    while(s->Next()){
      Int_t iDDL=rd->GetDDLID();
      Int_t iCarlos=s->GetCarlosId();
      if(s->IsCompletedModule()) continue;
      if(s->IsCompletedDDL()) continue;
      if(iDDL>=0 && iDDL<kTotDDL){ 
	Int_t index=kSides*(kModPerDDL*iDDL+iCarlos)+s->GetChannel(); 
	histo[index]->Fill(s->GetCoord2(),s->GetCoord1(),s->GetSignal());
	isFilled[index]=1;
      }
    }
    delete s;
    iev++;
    for(Int_t iddl=0; iddl<kTotDDL;iddl++){
      for(Int_t imod=0; imod<kModPerDDL;imod++){
	for(Int_t isid=0;isid<kSides;isid++){
	  Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	  anal[index]->AddEvent(histo[index]);
	  if(iddl==nDDL){
	    Int_t index2=kSides*imod+isid;
	    c0->cd(index2+1);
	    histo[index]->DrawCopy("colz");
	    sprintf(text,"DDL %d channel %d Side %d",nDDL,imod,isid);
	    t0->DrawLatex(0.15,0.92,text);
	    c0->Update();
	  }
	}
      }
    }
    printf(" --- OK\n");
  }while(rd->NextEvent()&&iev<=lastEv);

  TH1F *htotgain=new TH1F("htotgain","",100,0.2,4.2);
  TH1F *htotpeakpos=new TH1F("htotpeakpos","",256,-0.5,255.5);
  TH1F *hstatus=new TH1F("hstatus","",2,-0.5,1.5);

  TFile *outfil=new TFile("SDDgain-results.root","recreate");
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	if(isFilled[index]){
	  anal[index]->ValidateAnodes();
	  anal[index]->WriteToASCII();
	  anal[index]->WriteToROOT(outfil);
	  for(Int_t ian=0; ian<256;ian++){
	    Float_t gain=anal[index]->GetChannelGain(ian);
	    Float_t ppos=anal[index]->GetTimeBinTPPeak(ian);
	    Int_t anstatus=anal[index]->IsAnodeGood(ian);
	    hstatus->Fill(anstatus);
	    htotgain->Fill(gain);
	    htotpeakpos->Fill(ppos);
	  }
	}
      }
    }
  }
  outfil->Close();

  // Draw Statistics of baselines and noise
  TCanvas *call=new TCanvas("call","General stats",700,700);
  call->Divide(2,2);
  call->cd(1);
  htotpeakpos->Draw();
  htotpeakpos->GetXaxis()->SetTitle("TP peak position (Time Bin)");
  htotpeakpos->GetXaxis()->SetTitleSize(0.07);
  htotpeakpos->GetXaxis()->SetTitleOffset(0.6);
  call->cd(2);
  htotgain->Draw();
  htotgain->GetXaxis()->SetTitle("Gain (ADC/DAC)");
  htotgain->GetXaxis()->SetTitleSize(0.07);
  htotgain->GetXaxis()->SetTitleOffset(0.6);
  call->cd(3);
  hstatus->Draw();
  hstatus->GetXaxis()->SetTitle("Anode Status (0=bad 1=good)");
  hstatus->GetXaxis()->SetTitleSize(0.07);
  hstatus->GetXaxis()->SetTitleOffset(0.6);
  call->Update();
  call->SaveAs("GenStatsPulser.gif");

  // Draw baselines and noisegain and TestPulse time bin for all modules

  TH1F** hgain = new TH1F*[kModPerDDL*kSides];
  TH1F** htptb = new TH1F*[kModPerDDL*kSides];

  TCanvas *c1=new TCanvas("c1","DDL: TP position",900,900);
  c1->SetBottomMargin(0.14);
  c1->Divide(4,6,0.001,0.001);
  TCanvas *c2=new TCanvas("c2","DDL: gain",900,900);
  c2->SetBottomMargin(0.14);
  c2->Divide(4,6,0.001,0.001);

  for(Int_t imod=0; imod<kModPerDDL;imod++){
    for(Int_t isid=0;isid<kSides;isid++){
      Int_t index1=kSides*(kModPerDDL*nDDL+imod)+isid;
      Int_t index2=kSides*imod+isid;
      sprintf(text,"DDL %d channel %d Side %d",nDDL,imod,isid);

      TLatex *t3=new TLatex(0.15,0.92,text);
      t3->SetNDC();
      t3->SetTextSize(0.06);
      t3->SetTextColor(4);
      sprintf(hisnam,"hgain%ds%d",imod,isid);
      hgain[index2]=new TH1F(hisnam,"",256,-0.5,255.5);
      sprintf(hisnam,"htptb%ds%d",imod,isid);
      htptb[index2]=new TH1F(hisnam,"",256,-0.5,255.5);
      for(Int_t ian=0;ian<256;ian++){
	hgain[index2]->SetBinContent(ian+1,anal[index1]->GetChannelGain(ian));
	htptb[index2]->SetBinContent(ian+1,anal[index1]->GetTimeBinTPPeak(ian));
      }

      c1->cd(index2+1);
      htptb[index2]->Draw();
    //    htptb[imod]->SetMinimum(0);
    //    htptb[imod]->SetMaximum(75);
      htptb[index2]->GetXaxis()->SetTitle("Anode");
      htptb[index2]->GetYaxis()->SetTitle("TP position (Time Bin)");
      htptb[index2]->GetXaxis()->SetTitleSize(0.07);
      htptb[index2]->GetYaxis()->SetTitleSize(0.07);
      htptb[index2]->GetXaxis()->SetTitleOffset(0.6);
      htptb[index2]->GetYaxis()->SetTitleOffset(0.7);
      t3->Draw();
      c1->Update();


      c2->cd(index2+1); 
      hgain[index2]->SetMinimum(0.);
      hgain[index2]->SetMaximum(4.);
      hgain[index2]->Draw();
      hgain[index2]->GetXaxis()->SetTitle("Anode");
      hgain[index2]->GetYaxis()->SetTitle("Gain");
      hgain[index2]->GetXaxis()->SetTitleSize(0.07);
      hgain[index2]->GetYaxis()->SetTitleSize(0.07);
      hgain[index2]->GetXaxis()->SetTitleOffset(0.6);
      hgain[index2]->GetYaxis()->SetTitleOffset(0.7);
      hgain[index2]->SetStats(0);
      t3->Draw();
      c2->Update();
    }
  }

  c1->SaveAs("TPtimebin.gif");
  c2->SaveAs("Gain.gif");

}

void AnalyzeSDDGainAllMod(Int_t nrun, Int_t n2, Int_t year=2009, Char_t* dir="LHC09b_SDD",
			  Int_t adcfreq=20, 
			  Int_t nDDL=0, 
			  Int_t firstEv=18, 
			  Int_t lastEv=22, 
			  Float_t pascalDAC=100){


  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/%d/%s/%09d/raw/%02d%09d%03d.10.root",year,dir,nrun,year-2000,nrun,n2);
  printf("Open file %s\n",filnam);
  AnalyzeSDDGainAllMod(filnam,adcfreq,nDDL,firstEv,lastEv,pascalDAC);
}
