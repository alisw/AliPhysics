#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH2F.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TFile.h>
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSOnlineSDDBase.h"
#include "AliITSOnlineSDDCMN.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSDDCompressed.h"
#include "TPaveStats.h"
#endif

// Macro for the analysis of PEDESTAL runs (equivalent to ITSSDDBASda.cxx)
// Two functions named AnalyzeSDDNoiseAllModules: 
// The first is for analyzing a local raw data file and takes as agrument the file name.
// The second is for running on ALIEN
// All DDLs are analyzed, the argument nDDL selects the DDL to be plotted
// Origin: F. Prino (prino@to.infn.it)

void AnalyzeSDDNoiseAllMod(Char_t *datafil, 
			   Int_t adcfreq=20, 
			   Int_t nDDL=0, 
			   Int_t firstEv=18, 
			   Int_t lastEv=20){ 


  const Int_t kTotDDL=24;
  const Int_t kModPerDDL=12;
  const Int_t kSides=2;
  Bool_t writtenoutput=kFALSE;

  AliITSOnlineSDDBase **base=new AliITSOnlineSDDBase*[kTotDDL*kModPerDDL*kSides];
  TH2F **histo=new TH2F*[kTotDDL*kModPerDDL*kSides];

  Char_t hisnam[20];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	base[index]=new AliITSOnlineSDDBase(iddl,imod,isid);	
	if(adcfreq==40) base[index]->SetLastGoodTB(254);
	else base[index]->SetLastGoodTB(126);
	sprintf(hisnam,"h%02dc%02ds%d",iddl,imod,isid);
	histo[index]=new TH2F(hisnam,"",256,-0.5,255.5,256,-0.5,255.5);
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
    printf("Event # %d ",iev);
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
      }
    }
    delete s;
    iev++;
    for(Int_t iddl=0; iddl<kTotDDL;iddl++){
      for(Int_t imod=0; imod<kModPerDDL;imod++){
	for(Int_t isid=0;isid<kSides;isid++){
	  Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	  base[index]->AddEvent(histo[index]);
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

  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	base[index]->ValidateAnodes();
	base[index]->WriteToASCII(); // fondamentale!!!!!!!!!
	delete base[index];
      }
    }
  }
  delete rd;
  delete [] base;
  
  printf("Start second analysis for Common Mode correction\n");
  AliITSOnlineSDDCMN **corr=new AliITSOnlineSDDCMN*[kTotDDL*kModPerDDL*kSides];
  Bool_t isFilled[kTotDDL*kModPerDDL*kSides];

  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	corr[index]=new AliITSOnlineSDDCMN(iddl,imod,isid);
	if(adcfreq==40) corr[index]->SetLastGoodTB(254);
	else corr[index]->SetLastGoodTB(126);
	isFilled[index]=0;
      }
    }
  }

  iev=firstEv;
  AliRawReader *rd2; 
  if(strstr(datafil,".root")!=0){
    rd2=new AliRawReaderRoot(datafil,iev);
  }else{
    rd2=new AliRawReaderDate(datafil,iev);
  }
  do{
    c0->Clear();
    c0->Divide(4,6,0.001,0.001);
    printf("Event # %d ",iev);
    rd2->Reset();
    for(Int_t iddl=0; iddl<kTotDDL;iddl++){
      for(Int_t imod=0; imod<kModPerDDL;imod++){
	for(Int_t isid=0;isid<kSides;isid++){
	  Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	  histo[index]->Reset();
	}
      }
    }
    
    UChar_t cdhAttr=AliITSRawStreamSDD::ReadBlockAttributes(rd2);
    UInt_t amSamplFreq=AliITSRawStreamSDD::ReadAMSamplFreqFromCDH(cdhAttr);
    AliITSRawStream* s=AliITSRawStreamSDD::CreateRawStreamSDD(rd2,cdhAttr);
    if(!writtenoutput){
      printf("Use %s raw stream, sampling frequency %d MHz\n",s->ClassName(),amSamplFreq);
      writtenoutput=kTRUE;
    }
    while(s->Next()){
      Int_t iDDL=rd2->GetDDLID();
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
	  if(isFilled[index]) corr[index]->AddEvent(histo[index]);
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
  }while(rd2->NextEvent()&&iev<=lastEv);

  TH1F *htotbas=new TH1F("htotbas","",100,0.,150.);
  TH1F *htotbaseq=new TH1F("htotbaseq","",100,0.,150.);
  TH1F *htotnoise=new TH1F("htotnoise","",100,0.,10.);
  TH1F *htotnoisecorr=new TH1F("htotnoisecorr","",100,0.,10.);
  TH1F *hstatus=new TH1F("hstatus","",2,-0.5,1.5);

  TFile *outfil=new TFile("SDDbase-results.root","recreate");
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	if(isFilled[index]){
	  corr[index]->ValidateAnodes();
	  corr[index]->WriteToASCII();
	  corr[index]->WriteToROOT(outfil);
	  for(Int_t ian=0; ian<256;ian++){
	    Float_t basl=corr[index]->GetAnodeBaseline(ian);
	    Float_t basleq=corr[index]->GetAnodeEqualizedBaseline(ian);
	    Float_t noi=corr[index]->GetAnodeRawNoise(ian);
	    Float_t cornoi=corr[index]->GetAnodeCorrNoise(ian);
	    Int_t anstatus=corr[index]->IsAnodeGood(ian);
	    hstatus->Fill(anstatus);
	    htotbas->Fill(basl);
	    htotbaseq->Fill(basleq);
	    htotnoise->Fill(noi);
	    htotnoisecorr->Fill(cornoi);
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
  htotbas->Draw();
  htotbas->GetXaxis()->SetTitle("Baselines");
  htotbas->GetXaxis()->SetTitleSize(0.07);
  htotbas->GetXaxis()->SetTitleOffset(0.6);
  call->cd(2);
  htotbaseq->Draw();
  htotbaseq->GetXaxis()->SetTitle("Baselines after equalization");
  htotbaseq->GetXaxis()->SetTitleSize(0.07);
  htotbaseq->GetXaxis()->SetTitleOffset(0.6);
  call->cd(3);
  htotnoisecorr->SetLineColor(2);
  htotnoisecorr->Draw();
  call->Update();
  TPaveStats *st1=(TPaveStats*)htotnoisecorr->GetListOfFunctions()->FindObject("stats");
  st1->SetY1NDC(0.51);
  st1->SetY2NDC(0.7);
  htotnoisecorr->GetXaxis()->SetTitle("Noise");
  htotnoisecorr->GetXaxis()->SetTitleSize(0.07);
  htotnoisecorr->GetXaxis()->SetTitleOffset(0.6);
  htotnoise->Draw("SAMES");
  call->Update();
  TPaveStats *st2=(TPaveStats*)htotnoise->GetListOfFunctions()->FindObject("stats");
  st2->SetY1NDC(0.71);
  st2->SetY2NDC(0.9);
  
  call->cd(4);
  hstatus->Draw();
  hstatus->GetXaxis()->SetTitle("Anode Status (0=bad 1=good)");
  hstatus->GetXaxis()->SetTitleSize(0.07);
  hstatus->GetXaxis()->SetTitleOffset(0.6);
  call->Update();
  call->SaveAs("GenStatsPedestal.gif");

  // Draw baselines and noise for all modules of the selected DDL

  TH1F** hbas = new TH1F*[kSides*kModPerDDL];
  TH1F** hrawn = new TH1F*[kSides*kModPerDDL];
  TH1F** hcorrn = new TH1F*[kSides*kModPerDDL];
  TH1F** hdbas = new TH1F*[kSides*kModPerDDL];
  TH1F** hdrawn = new TH1F*[kSides*kModPerDDL];
  TH1F** hdcorrn = new TH1F*[kSides*kModPerDDL];
  TCanvas *c1=new TCanvas("c1","DDL: Baselines vs anode",900,900);
  c1->SetBottomMargin(0.14);
  c1->Divide(4,6,0.001,0.001);
  TCanvas *c2=new TCanvas("c2","DDL: Noise vs anode",900,900);
  c2->SetBottomMargin(0.14);
  c2->Divide(4,6,0.001,0.001);
  TCanvas *c3=new TCanvas("c3","DDL: Baselines distr",900,900);
  c3->SetBottomMargin(0.14);
  c3->Divide(4,6,0.001,0.001);
  TCanvas *c4=new TCanvas("c4","DDL: Noise Distr",900,900);
  c4->SetBottomMargin(0.14);
  c4->Divide(4,6,0.001,0.001);
  TLatex *t1=new TLatex(0.15,0.2,"Raw Noise");
  t1->SetNDC();
  t1->SetTextSize(0.05);
  TLatex *t2=new TLatex(0.4,0.2,"Corrected Noise");
  t2->SetNDC();
  t2->SetTextSize(0.05);
  t2->SetTextColor(2);
  TLatex *t3=new TLatex();
  t3->SetNDC();
  t3->SetTextSize(0.06);
  t3->SetTextColor(4);

  for(Int_t imod=0; imod<kModPerDDL;imod++){
    for(Int_t isid=0;isid<kSides;isid++){
      Int_t index1=kSides*(kModPerDDL*nDDL+imod)+isid;
      Int_t index2=kSides*imod+isid;
      sprintf(text,"DDL %d channel %d Side %d",nDDL,imod,isid);
      hbas[index2]=corr[index1]->GetBaselineAnodeHisto();
      hrawn[index2]=corr[index1]->GetRawNoiseAnodeHisto();
      hcorrn[index2]=corr[index1]->GetCorrNoiseAnodeHisto();
      hdbas[index2]=corr[index1]->GetBaselineHisto();
      hdrawn[index2]=corr[index1]->GetRawNoiseHisto();
      hdcorrn[index2]=corr[index1]->GetCorrNoiseHisto();
      c1->cd(index2+1);
      hbas[index2]->Draw();
      hbas[index2]->SetMinimum(0);
      hbas[index2]->SetMaximum(75);
      hbas[index2]->GetXaxis()->SetTitle("Anode");
      hbas[index2]->GetYaxis()->SetTitle("Baseline");
      hbas[index2]->GetXaxis()->SetTitleSize(0.07);
      hbas[index2]->GetYaxis()->SetTitleSize(0.07);
      hbas[index2]->GetXaxis()->SetTitleOffset(0.6);
      hbas[index2]->GetYaxis()->SetTitleOffset(0.7);
      t3->DrawLatex(0.15,0.92,text);
      c1->Update();


      c2->cd(index2+1); 
      hrawn[index2]->SetMinimum(1.);
      hrawn[index2]->SetMaximum(6.);
      hrawn[index2]->Draw();
      hrawn[index2]->GetXaxis()->SetTitle("Anode");
      hrawn[index2]->GetYaxis()->SetTitle("Noise");
      hrawn[index2]->GetXaxis()->SetTitleSize(0.07);
      hrawn[index2]->GetYaxis()->SetTitleSize(0.07);
      hrawn[index2]->GetXaxis()->SetTitleOffset(0.6);
      hrawn[index2]->GetYaxis()->SetTitleOffset(0.7);
      gStyle->SetOptStat(0);
      hrawn[index2]->SetStats(0);
      hcorrn[index2]->SetLineColor(2);
      hcorrn[index2]->Draw("SAME");
      t1->Draw();
      t2->Draw();
      t3->DrawLatex(0.15,0.92,text);
      c2->Update();

      c3->cd(index2+1);
      hdbas[index2]->Draw();
      hdbas[index2]->GetXaxis()->SetTitle("Baseline");
      hdbas[index2]->GetXaxis()->SetTitleSize(0.07);
      hdbas[index2]->GetXaxis()->SetTitleOffset(0.6);
      t3->DrawLatex(0.15,0.92,text);
      c3->Update();

      c4->cd(index2+1); 
      hdrawn[index2]->Draw();
      hdrawn[index2]->GetXaxis()->SetTitle("Noise");
      hdrawn[index2]->GetXaxis()->SetTitleSize(0.07);
      hdrawn[index2]->GetXaxis()->SetTitleOffset(0.6);
      hdcorrn[index2]->SetLineColor(2);
      hdcorrn[index2]->Draw("SAME");
      t1->Draw();
      t2->Draw();
      t3->DrawLatex(0.15,0.92,text);
      c4->Update();
    }
  }

  c1->SaveAs("Baselines.gif");
  c2->SaveAs("Noise.gif");
  c3->SaveAs("BaselinesDist.gif");
  c4->SaveAs("NoiseDist.gif");
  
  Char_t delfil[100];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	sprintf(delfil,"rm SDDbase_step1_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	gSystem->Exec(delfil);
      }
    }
  }

}

void AnalyzeSDDNoiseAllMod(Int_t nrun, Int_t n2, Int_t year=2009, Char_t* dir="LHC09b_SDD",
			   Int_t adcfreq=20, 
			   Int_t nDDL=0, 
			   Int_t firstEv=18, 
			   Int_t lastEv=20){

  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/%d/%s/%09d/raw/%02d%09d%03d.10.root",year,dir,nrun,year-2000,nrun,n2);
  printf("Open file %s\n",filnam);
  AnalyzeSDDNoiseAllMod(filnam,adcfreq,nDDL,firstEv,lastEv);
}
