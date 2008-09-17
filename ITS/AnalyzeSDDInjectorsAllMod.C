#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TFile.h>
#include <TGrid.h>
#include <TF1.h>
#include <TLine.h>
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSOnlineSDDInjectors.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSDDLModuleMapSDD.h"
#endif

// Macro for the analysis of PULSER runs (equivalent to ITSSDDINJda.cxx)
// Two functions named AnalyzeSDDInjectorsAllModules: 
// The first is for analyzing a local raw data file and takes as agrument the file name.
// The second is for running on ALIEN
// All DDLs are analyzed, the argument nDDL selects the DDL to be plotted
// Origin: F. Prino (prino@to.infn.it)


void AnalyzeSDDInjectorsAllMod(Char_t *datafil, Int_t nDDL, Int_t firstEv=10, Int_t lastEv=15){

  const Int_t kTotDDL=24;
  const Int_t kModPerDDL=12;
  const Int_t kSides=2;

  AliITSDDLModuleMapSDD* dmap=new AliITSDDLModuleMapSDD();
  dmap->SetJun08Map();

  TH2F** histo = new TH2F*[kTotDDL*kModPerDDL*kSides];
  Int_t nWrittenEv[kTotDDL*kModPerDDL*kSides];
  TGraphErrors** gvel = new TGraphErrors*[kTotDDL*kModPerDDL*kSides];
  AliITSOnlineSDDInjectors **anal=new AliITSOnlineSDDInjectors*[kTotDDL*kModPerDDL*kSides];

  Char_t hisnam[20];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	sprintf(hisnam,"h%02dc%02ds%d",iddl,imod,isid);
	histo[index]=new TH2F(hisnam,"",256,-0.5,255.5,256,-0.5,255.5);
	anal[index]=new AliITSOnlineSDDInjectors(iddl,imod,isid);
/* Uncomment these lines for analysis of runs with 40 MHz sapling */
// 	anal[index]->SetInjLineRange(0,20,50);
// 	anal[index]->SetInjLineRange(1,90,160);
// 	anal[index]->SetInjLineRange(2,170,240);
// 	anal[index]->SetTimeStep(25.);
/* END of lines to be uncommented */
	nWrittenEv[index]=0;
      }
    }
  }
  TGraph *gvvsmod0=new TGraph(0);
  TGraph *gvvsmod1=new TGraph(0);

  TCanvas* c0 = new TCanvas("c0","",900,900);
  gStyle->SetPalette(1);
  TCanvas* c1 = new TCanvas("c1","",900,900);
  Char_t text[50];

  Int_t iev=firstEv;
  AliRawReader *rd; 
  if(strstr(datafil,".root")!=0){
    rd=new AliRawReaderRoot(datafil,iev);
  }else{
    rd=new AliRawReaderDate(datafil,iev);
  }

  Char_t gname[15];
  TF1 *funz=new TF1("funz","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.,255.);
  TLatex *t0=new TLatex();
  t0->SetNDC();
  t0->SetTextSize(0.06);
  t0->SetTextColor(4);

  do{
    c0->Clear();
    c0->Divide(4,6,0.001,0.001);
    c1->Clear();
    c1->Divide(4,6,0.001,0.001);
    printf("Event # %d\n",iev);
    UInt_t timeSt=rd->GetTimestamp();
    rd->Reset();
    for(Int_t iddl=0; iddl<kTotDDL;iddl++){
      for(Int_t imod=0; imod<kModPerDDL;imod++){
	for(Int_t isid=0;isid<kSides;isid++){
	  Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	  histo[index]->Reset();
	}
      }
    }

    AliITSRawStreamSDD s(rd);
    while(s.Next()){
      Int_t iDDL=rd->GetDDLID();
      Int_t iCarlos=s.GetCarlosId();
      if(iDDL>=0 && iDDL<kTotDDL && s.IsCompletedModule()==kFALSE){ 
	Int_t index=kSides*(kModPerDDL*iDDL+iCarlos)+s.GetChannel(); 
	histo[index]->Fill(s.GetCoord2(),s.GetCoord1(),s.GetSignal());
      }
    }
    
    for(Int_t iddl=0; iddl<kTotDDL;iddl++){
      for(Int_t imod=0; imod<kModPerDDL;imod++){
	for(Int_t isid=0;isid<kSides;isid++){
	  Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	  anal[index]->AnalyzeEvent(histo[index]); 
	  anal[index]->WriteToASCII(iev,timeSt,nWrittenEv[index]);
	  nWrittenEv[index]++;
	  if(iev==firstEv && anal[index]->GetInjPadStatus(16)>=6){
	    Float_t vel=anal[index]->GetDriftSpeed(16);
	    Int_t iMod=dmap->GetModuleNumber(iddl,imod);
	    if(isid==0) gvvsmod0->SetPoint(gvvsmod0->GetN(),(Float_t)iMod,vel);
	    if(isid==1) gvvsmod1->SetPoint(gvvsmod1->GetN(),(Float_t)iMod,vel);
	  }
	  if(iddl==nDDL){
	    Int_t index2=kSides*imod+isid;
	    c0->cd(index2+1);
	    histo[index]->SetMaximum(100.);
	    histo[index]->DrawCopy("colz");
	    sprintf(text,"DDL %d channel %d Side %d",nDDL,imod,isid);
	    t0->DrawLatex(0.15,0.92,text);
	    c0->Update();
	    c1->cd(index2+1);
	    gvel[index]=anal[index]->GetDriftSpeedGraph();
	    gvel[index]->SetMarkerStyle(20);
	    gvel[index]->SetTitle("");
	    sprintf(gname,"gvel%dev%d",index,iev);
	    gvel[index]->SetName(gname);
	    //	    gvel[index]->SetMinimum(0);
	    //gvel[index]->SetMaximum(200);

	    gvel[index]->GetXaxis()->SetLimits(0,256);
	    gvel[index]->GetXaxis()->SetTitle("Anode");
	    gvel[index]->GetYaxis()->SetTitle("Drift vel.");
	    gvel[index]->GetXaxis()->SetTitleSize(0.07);
	    gvel[index]->GetYaxis()->SetTitleSize(0.07);
	    gvel[index]->GetXaxis()->SetTitleOffset(0.6);
	    gvel[index]->GetYaxis()->SetTitleOffset(0.6);
	    if(gvel[index]->GetN()>0) gvel[index]->Draw("AP");
	    Float_t *param=anal[index]->GetDriftSpeedFitParam();
	    funz->SetParameters(param[0],param[1],param[2],param[3]);
	    funz->SetLineColor(2);
	    funz->DrawCopy("LSAME");
	    t0->DrawLatex(0.15,0.92,text);
	    c1->Update();
	  }
	}
      }
    }
    iev++;
    printf(" --- OK\n");
  }while(rd->NextEvent()&&iev<=lastEv);

  TCanvas* c8=new TCanvas("c8");
  gvvsmod0->SetTitle("");
  gvvsmod1->SetTitle("");

  gvvsmod0->SetMarkerStyle(20);
  gvvsmod1->SetMarkerStyle(21);
  gvvsmod1->SetMarkerColor(2);
  gvvsmod0->Draw("AP");
  gvvsmod0->SetMinimum(6.2);
  gvvsmod0->SetMaximum(7.2);
  gvvsmod0->GetXaxis()->SetTitle("Module Number");
  gvvsmod0->GetYaxis()->SetTitle("Vdrift at injector pad 16");  
  gvvsmod1->Draw("PSAME");
  TLatex* tleft=new TLatex(0.7,0.82,"Side 0");
  tleft->SetNDC();
  tleft->SetTextColor(1);
  tleft->Draw();
  TLatex* tright=new TLatex(0.7,0.75,"Side 1");
  tright->SetNDC();
  tright->SetTextColor(2);
  tright->Draw();

  TLine *lin=new TLine(323,6.2,323,7.2);
  lin->SetLineColor(4);
  lin->Draw();
  c8->Update();
}

void AnalyzeSDDInjectorsAllMod(Int_t nrun, Int_t n2, Int_t nDDL=0, Int_t firstEv=10, Int_t lastEv=15){
  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/2008/LHC08c_SDD/%09d/raw/08%09d%03d.10.root",nrun,nrun,n2);
  printf("Open file %s\n",filnam);
  AnalyzeSDDInjectorsAllMod(filnam,nDDL,firstEv,lastEv);
}



