#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TFile.h>
#include <TMath.h>
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


void AnalyzeSDDInjectorsAllMod(Char_t *datafil, Int_t nDDL, Int_t firstEv=10, Int_t lastEv=15,Int_t jpad=16, Int_t statuscut=7){

  const Int_t kTotDDL=24;
  const Int_t kModPerDDL=12;
  const Int_t kSides=2;

  AliITSDDLModuleMapSDD* dmap=new AliITSDDLModuleMapSDD();
  dmap->SetJun08Map();

  TH2F** histo = new TH2F*[kTotDDL*kModPerDDL*kSides];
  Int_t nWrittenEv[kTotDDL*kModPerDDL*kSides];
  TGraphErrors** gvel = new TGraphErrors*[kTotDDL*kModPerDDL*kSides];
  AliITSOnlineSDDInjectors **anal=new AliITSOnlineSDDInjectors*[kTotDDL*kModPerDDL*kSides];
  TH1F** hvdriftl=new TH1F*[260];  
  TH1F** hvdriftr=new TH1F*[260];  
  Char_t hisnam[20];
  for(Int_t idet=0; idet<260;idet++){
    sprintf(hisnam,"vdriftl%03d",idet);
    hvdriftl[idet]=new TH1F(hisnam,"",500,5.5,8.0);
    sprintf(hisnam,"vdriftr%03d",idet);
    hvdriftr[idet]=new TH1F(hisnam,"",500,5.5,8.0);
  }
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
  TGraphErrors *gvvsmod0=new TGraphErrors(0);
  TGraphErrors *gvvsmod1=new TGraphErrors(0);
  TGraphErrors *gtvsmod0=new TGraphErrors(0);
  TGraphErrors *gtvsmod1=new TGraphErrors(0);
  Float_t gvmin=6.0, gvmax=7.5;
  Float_t gtmin=288., gtmax=308.;
  TH1F* hanst=new TH1F("hanst","",8,-0.5,7.5);
  TH1F* hpad7l=new TH1F("hpad7l","",33,-0.5,32.5);
  TH1F* hpad7r=new TH1F("hpad7r","",33,-0.5,32.5);

  TCanvas* c0 = new TCanvas("c0","Event display",900,900);
  gStyle->SetPalette(1);
  TCanvas* c1 = new TCanvas("c1","Drift Speed vs. anode",900,900);
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
  Int_t readEv=0;
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
	  Int_t iMod=dmap->GetModuleNumber(iddl,imod);
	  if(iMod!=-1){
	    for(Int_t ipad=0;ipad<33;ipad++){
	      Int_t st=anal[index]->GetInjPadStatus(ipad);
	      hanst->Fill(st);
	      if(anal[index]->GetInjPadStatus(ipad)>=statuscut){
		if(isid==0) hpad7l->Fill(ipad);
		if(isid==1) hpad7r->Fill(ipad);
	      }
	    }
	    if(anal[index]->GetInjPadStatus(jpad)>=statuscut){
	      Float_t vel=anal[index]->GetDriftSpeed(jpad);
	      if(isid==0) hvdriftl[iMod-240]->Fill(vel);
	      if(isid==1) hvdriftr[iMod-240]->Fill(vel);
	    }
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
    readEv++;
    printf(" --- OK\n");
  }while(rd->NextEvent()&&iev<=lastEv);
  printf("Total number of events = %d\n",readEv);
  Float_t nfac=1./(Float_t)readEv/33./520.;
  hanst->Scale(nfac);
  nfac=1./(Float_t)readEv;
  hpad7l->Scale(nfac);
  hpad7r->Scale(nfac);

  Int_t ipt0=0, ipt1=0;
  Float_t Edrift=(1800-45)/291/0.012;  
  TFile *outfil=new TFile("DriftSpeedHistos.root","recreate");
  for(Int_t iMod=0; iMod<260; iMod++){
    outfil->cd();
    hvdriftl[iMod]->Write();    
    hvdriftr[iMod]->Write();
    Float_t modid=iMod+240;
    if(hvdriftl[iMod]->GetEntries()>0){
      Float_t avevell=hvdriftl[iMod]->GetMean();
      Float_t rmsvell=hvdriftl[iMod]->GetRMS();
      if(avevell > 5.5 && avevell < 8.5){
	gvvsmod0->SetPoint(ipt0,modid,avevell);
	gvvsmod0->SetPointError(ipt0,0,rmsvell);
	Float_t mob=avevell*1.E5/Edrift;  
	Float_t temper=293.15*TMath::Power((mob/1350.),-1/2.4); 
	gtvsmod0->SetPoint(ipt0,modid,temper);
	++ipt0;
      }
    }
    if(hvdriftr[iMod]->GetEntries()>0){
      Float_t avevelr=hvdriftr[iMod]->GetMean();
      Float_t rmsvelr=hvdriftr[iMod]->GetRMS();
      if(avevelr > 5.5 && avevelr < 8.5){
	gvvsmod1->SetPoint(ipt1,modid,avevelr);
	gvvsmod1->SetPointError(ipt1,0,rmsvelr);
	Float_t mob=avevelr*1.E5/Edrift;
	Float_t temper=293.15*TMath::Power((mob/1350.),-1./2.4); 
	gtvsmod1->SetPoint(ipt1,modid,temper);
	++ipt1;
      }
    }
  }
  gvvsmod0->SetName("gvvsmod0");
  gvvsmod1->SetName("gvvsmod1");
  gtvsmod0->SetName("gtvsmod0");
  gtvsmod1->SetName("gtvsmod1");
  outfil->cd();
  gvvsmod0->Write();
  gvvsmod1->Write();
  gtvsmod0->Write();
  gtvsmod1->Write();
  outfil->Close();

  TCanvas* c8=new TCanvas("c8","Drift Speed vs. mod");
  gvvsmod0->SetTitle("");
  gvvsmod1->SetTitle("");
  gvvsmod0->SetMarkerStyle(20);
  gvvsmod1->SetMarkerStyle(21);
  gvvsmod1->SetMarkerColor(2);
  gvvsmod0->Draw("AP");
  gvvsmod0->SetMinimum(gvmin);
  gvvsmod0->SetMaximum(gvmax);
  gvvsmod0->GetXaxis()->SetTitle("Module Number");
  Char_t title[25];
  sprintf(title,"Vdrift at injector pad %d",jpad);
  gvvsmod0->GetYaxis()->SetTitle(title);  
  gvvsmod1->Draw("PSAME");
  TLatex* tleft=new TLatex(0.7,0.82,"Side 0");
  tleft->SetNDC();
  tleft->SetTextColor(1);
  tleft->Draw();
  TLatex* tright=new TLatex(0.7,0.75,"Side 1");
  tright->SetNDC();
  tright->SetTextColor(2);
  tright->Draw();

  TLine *lin=new TLine(323,gvmin,323,gvmax);
  lin->SetLineColor(4);
  lin->Draw();
  c8->Update();
  c8->SaveAs("VdriftVsMod.gif");

  TCanvas* c8t=new TCanvas("c8t","Temeprature vs. mod");
  gtvsmod0->SetTitle("");
  gtvsmod1->SetTitle("");
  gtvsmod0->SetMarkerStyle(20);
  gtvsmod1->SetMarkerStyle(21);
  gtvsmod1->SetMarkerColor(2);
  gtvsmod0->Draw("AP");
  gtvsmod0->SetMinimum(gtmin);
  gtvsmod0->SetMaximum(gtmax);
  gtvsmod0->GetXaxis()->SetTitle("Module Number");
  sprintf(title,"Estimated Temperature (K)");
  gtvsmod0->GetYaxis()->SetTitle(title);  
  gtvsmod1->Draw("PSAME");
  tleft->Draw();
  tright->Draw();
  TLine *lint=new TLine(323,gtmin,323,gtmax);
  lint->SetLineColor(4);
  lint->Draw();
  c8t->Update();
  c8t->SaveAs("TempVsMod.gif");

  TCanvas* c9=new TCanvas("c9","Injector status");
  hanst->SetStats(0);
  hanst->Draw();
  hanst->GetXaxis()->SetTitle("Injector pad status");
  hanst->GetXaxis()->CenterTitle();
  c9->SaveAs("InjStatus.gif");

//   TCanvas* c10=new TCanvas("c10","Pad status 7",1200,600);
//   hpad7l->SetStats(0);
//   hpad7r->SetStats(0);
//   c10->Divide(2,1);
//   c10->cd(1);
//   hpad7l->Draw(); 
//   hpad7l->GetXaxis()->SetTitle("Side Left -- Pad number");
//   hpad7l->GetXaxis()->CenterTitle();
//   hpad7l->GetYaxis()->SetTitle("Number of status 7");
//   c10->cd(2);
//   hpad7r->Draw();
//   hpad7r->GetXaxis()->SetTitle("Side Right -- Pad number");
//   hpad7r->GetXaxis()->CenterTitle();
//   hpad7r->GetYaxis()->SetTitle("Number of status 7");
//   printf("Side 0, maximum pad=%d\n",hpad7l->GetMaximumBin());
//   printf("Side 1, maximum pad=%d\n",hpad7r->GetMaximumBin());
  

}

void AnalyzeSDDInjectorsAllMod(Int_t nrun, Int_t n2, Char_t* dir="LHC08d_SDD", Int_t nDDL=0, Int_t firstEv=15, Int_t lastEv=15){
  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/2008/%s/%09d/raw/08%09d%03d.10.root",dir,nrun,nrun,n2);
  printf("Open file %s\n",filnam);
  AnalyzeSDDInjectorsAllMod(filnam,nDDL,firstEv,lastEv);
}



