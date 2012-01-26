/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//Draw the 2-particle correlation histograms (from output of AliAnalysisTaskDiHadron)
//Author: Jason Glyndwr Ulery, ulery@uni-frankfurt.de

#ifndef __CINT__
#include <TF1.h>
#include <TF2.h>
#include "TSystem.h"
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TRandom.h>
#include <fstream>
#include "TList.h"
#include "TArrayF.h"
#endif

void DrawDPhi(float TPt1, float TPt2, float APt1, float APt2, float Cent1, float Cent2, int ptWeighted=0, int notDelta=0, int LSign=0){//ptWeighted 0 no weight, 1 weighted, 2 <pt>; Lsign 0-all, 1-like-sign 2-unlike-sign
  gROOT->Reset();
  gROOT->LoadMacro("Util9.C");
  Style(1);
  //gDebug=1;
  int SaveRoot=1;
  int SavePS=0;
  int SavePDF=0;
  float ZYACent=1.5;//0  for 2-D zyam lowest bins -number wanted for 2-D strips//1.5
  float ZYAWidth=0.2; 
  int ReBinPhi=1;
  int ReBinEta=1;
  int DrawAsHisto=0;
  int noTitle=1;
  int SaveMacro=0;
  int SaveText=0;
  int DrawMC=0;
  int EffMethod=0;//0 no correction, 1 angle dependent, 2 use mix for triggered, 3 <1>, 4 fits, 5  fits with VariablePtLimit
  //int SaveEffFits=1;

 
  
  //char *Folder="2010-03-10/pass5";
  //char *Folder="2010-04-07/7_2010_04_05";
  //char *EffFolder="2010-03-10/all_5";

  //char *Folder="2010-04-12/pass5";
  //char *EffFolder="2010-04-12/900Pythia";

  // char *Folder="2010-04-26/7pass1";
  //char *Folder="2010-04-26/7Pythia";
  //char *Folder="2010-05-28/7Pythia2";//with new cut
  //char *Folder="2010-05-28/900Pythia";//with new cut
  //char *Folder="2010-05-28/7Pythia_LHC10b5";//newer pythia
  //har *Folder="2010-06-30/7Pythia_LHC10b5";//newer pythia
  /*
  char *Folder="2010-07-02/LHC10c_900pass2";
  //char *Folder="2010-07-02/LHC10c6_900Pythia";
  char *EffFolder="2010-07-02/LHC10c6_900Pythia";
  */
  //12&21
  //7-22 paper cuts, 7-27 papercuts and early corr, 7-12 old cuts, 8-4 more multiplicity cuts
  char *Folder="2010-08-17/LHC10b_7pass2";
  // char *Folder="2010-08-17/7Pythia_LHC10b5";
  //char *Folder="2010-08-17/LHC10c_900pass2";
  //char *Folder="2010-08-17/LHC10c6_900Pythia";

   char *EffFolder="2010-08-17/7Pythia_LHC10b5";
   //char *EffFolder="2010-08-17/LHC10c6_900Pythia";

  //  char *EffFolder="2010-04-26/7Pythia";
  // char *EffFolder="2010-05-28/7Pythia2";//with new cut
  // char *EffFolder="2010-05-28/900Pythia";//with new cut
  //char *EffFolder="2010-06-30/7Pythia_LHC10b5";//newer pythia
  //char *EffFolder="2010-05-28/7Pythia_LHC10b5";//newer pythia


  char *cPt[3]={"","Pt","MPt"};
  char *cDelta[2]={"Delta",""};
  char *csign[3]={"","LS","ULS"};
  char name[300];
  char inName[100];
  char effName[100];
  sprintf(inName,"%s/julery_DiHadron.root",Folder);
  sprintf(effName,"%s/julery_DiHadron.root",EffFolder);
  cout << inName << endl;
  cout << effName << endl;
  char *titArray[2]={"","_NoTitle"};
  char *histArray[2]={"","_Histo"};
  
  Float_t MPt[4];//<pt>, error <pt>, # of triggers
  Float_t MPt2[4];
  Float_t TrigSum=0, TrigSum2=0;
  
  TH1F *hPhiRaw=new TH1F("hPhiRaw","",1,0,1);
  TH1F *hPhiCorr=new TH1F("hPhiCorr","",1,0,1);
  TH1F *hPhiEff=new TH1F("hPhiEff","",1,0,1);
  TH1F *hPhiMC=new TH1F("hPhiMC","",1,0,1);
  TH1F *hPhiMixRaw=new TH1F("hPhiMixRaw","",1,0,1);
  TH1F *hPhiMixCorr=new TH1F("hPhiMixCorr","",1,0,1);
  TH1F *hPhiMixEff=new TH1F("hPhiMixEff","",1,0,1);
  TH1F *hPhiMixMC=new TH1F("hPhiMixMC","",1,0,1);

  TH1F *hEtaNRaw=new TH1F("hEtaNRaw","",1,0,1);
  TH1F *hEtaNCorr=new TH1F("hEtaNCorr","",1,0,1);
  TH1F *hEtaNEff=new TH1F("hEtaNEff","",1,0,1);
  TH1F *hEtaNMC=new TH1F("hEtaNMC","",1,0,1);
  TH1F *hEtaNMixRaw=new TH1F("hEtaNMixRaw","",1,0,1);
  TH1F *hEtaNMixCorr=new TH1F("hEtaNMixCorr","",1,0,1);
  TH1F *hEtaNMixEff=new TH1F("hEtaNMixEff","",1,0,1);
  TH1F *hEtaNMixMC=new TH1F("hEtaNMixMC","",1,0,1);

  TH1F *hEtaARaw=new TH1F("hEtaARaw","",1,0,1);
  TH1F *hEtaACorr=new TH1F("hEtaACorr","",1,0,1);
  TH1F *hEtaAEff=new TH1F("hEtaAEff","",1,0,1);
  TH1F *hEtaAMC=new TH1F("hEtaAMC","",1,0,1);
  TH1F *hEtaAMixRaw=new TH1F("hEtaAMixRaw","",1,0,1);
  TH1F *hEtaAMixCorr=new TH1F("hEtaAMixCorr","",1,0,1);
  TH1F *hEtaAMixEff=new TH1F("hEtaAMixEff","",1,0,1);
  TH1F *hEtaAMixMC=new TH1F("hEtaAMixMC","",1,0,1);

  TH2F *hPhiEtaRaw=new TH2F("hPhiEtaRaw","",1,0,1,1,0,1);
  TH2F *hPhiEtaCorr=new TH2F("hPhiEtaCorr","",1,0,1,1,0,1);
  TH2F *hPhiEtaEff=new TH2F("hPhiEtaEff","",1,0,1,1,0,1);
  TH2F *hPhiEtaMC=new TH2F("hPhiEtaMC","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixRaw=new TH2F("hPhiEtaMixRaw","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixCorr=new TH2F("hPhiEtaMixCorr","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixEff=new TH2F("hPhiEtaMixEff","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixMC=new TH2F("hPhiEtaMixMC","",1,0,1,1,0,1);
  TH1F *hMult=new TH1F("hMult","",1,0,1);
  TH1F *hTMult=new TH1F("hTMult","",1,0,1);
  TH1F *hNTrig=new TH1F("hNTrig","Number of Triggers",1,-0.5,0.5);
  //TH1F *hNEvents=new TH1F("hNEvents","Number of Events & In Selected Cent",2,-0.5,1.5);

  TF1 *ZeroLine=new TF1("ZeroLine","[0]",-2,7);
  ZeroLine->SetParameter(0,0);
  ZeroLine->SetLineStyle(2);
   TF1 *fit1=new TF1("fit1","1/sqrt(2*3.1415926)*([0]/[1]*(exp(-0.5*pow(x/[1],2))+exp(-0.5*pow((x-6.29185)/[1],2)))+[2]/[3]*(exp(-0.5*pow((x-3.14159)/[3],2))+exp(-0.5*pow((x+3.14159)/[3],2))))+[4]");
  fit1->SetParameters(.1,.2,.1,.35,.1);
  fit1->SetParNames("Near-Yield", "Near-Width","Away-Yield","Away-Width","Bg-Level");
   TF1 *fit2=new TF1("fit2","1/sqrt(2*3.1415926)*([0]/[1]*(exp(-0.5*pow(x/[1],2))))+[2]");
  fit2->SetParameters(.1,.2,.1);
  fit2->SetParNames("Near-Yield", "Near-Width","Bg-Level");

 int VariablePtLimit=0;//if 1 and EffMethod==4 then a variable upper limit on the associated pt is used
 if(EffMethod==5)VariablePtLimit=1;
  
  TFile *inFile=new TFile(inName);
  TList *inList=inFile->Get("julery_DiHadron");
  //sprintf(xname,"fHistNEvents_C%d",xCent);
  TH1F *hNEvents=(TH1F*)inList->FindObject("fHistNEvents_C0");
  //cout << "Number of Events: " << hNEvents->GetBinContent(1) << endl;
  TFile *effFile=new TFile(effName);
  TList *effList=effFile->Get("julery_DiHadron");
  for(int Cent=Cent1;Cent<=Cent2;Cent++){
    int qqx=0;
    if(ptWeighted)qqx=1;
    for(int ipt=0;ipt<=qqx;ipt++){
    MakeProjections(TPt1,TPt2,APt1,APt2,Cent,inList,hPhiRaw,hPhiMixRaw,hEtaNRaw,hEtaNMixRaw,hEtaARaw,hEtaAMixRaw,hPhiEtaRaw,hPhiEtaMixRaw,MPt,ipt,0,notDelta,LSign);
    MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiEff,hPhiMixEff,hEtaNEff,hEtaNMixEff,hEtaAEff,hEtaAMixEff,hPhiEtaEff,hPhiEtaMixEff,MPt2,ipt,0,notDelta,LSign);
    MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiMC,hPhiMixMC,hEtaNMC,hEtaNMixMC,hEtaAMC,hEtaAMixMC,hPhiEtaMC,hPhiEtaMixMC,MPt2,ipt,1,notDelta,LSign);
    c100xx=new TCanvas("c100");
    hPhiMixEff->Draw();
    c100xx->SaveAs("Mix.gif");
    if(EffMethod<4){
      EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hPhiEff,hPhiMC,hPhiMixEff,hPhiMixMC,hEtaNEff,hEtaNMC,hEtaNMixEff,hEtaNMixMC,hEtaAEff,hEtaAMC,hEtaAMixEff,hEtaAMixMC,hPhiEtaEff,hPhiEtaMC,hPhiEtaMixEff,hPhiEtaMixMC,EffMethod);
    }
    else{
      EffFit(APt1,APt2,Cent,effList,hPhiEff,hPhiMixEff,hEtaNEff,hEtaNMixEff,hEtaAEff,hEtaAMixEff,hPhiEtaEff,hPhiEtaMixEff,LSign,VariablePtLimit);
    }
    //EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hPhiRaw,hPhiEff,hPhiMC,hPhiMixRaw,hPhiMixEff,hPhiMixMC,hEtaNRaw,hEtaNEff,hEtaNMC,hEtaNMixRaw,hEtaNMixEff,hEtaNMixMC,hEtaARaw,hEtaAEff,hEtaAMC,hEtaAMixRaw,hEtaAMixEff,hEtaAMixMC,hPhiEtaRaw,hPhiEtaEff,hPhiEtaMC,hPhiEtaMixRaw,hPhiEtaMixEff,hPhiEtaMixMC,EffMethod);
    
    hPhiCorr=(TH1F*)hPhiRaw->Clone();
    sprintf(name,"hPhiCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hPhiCorr->SetName(name);
    hPhiCorr->Divide(hPhiEff);
    
    hPhiMixCorr=(TH1F*)hPhiMixRaw->Clone();
    sprintf(name,"hPhiMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hPhiMixCorr->SetName(name);
    hPhiMixCorr->Divide(hPhiMixEff);

    hEtaNCorr=(TH1F*)hEtaNRaw->Clone();
    sprintf(name,"hEtaNCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hEtaNCorr->SetName(name);
    hEtaNCorr->Divide(hEtaNEff);
    
    hEtaNMixCorr=(TH1F*)hEtaNMixRaw->Clone();
    sprintf(name,"hEtaNMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hEtaNMixCorr->SetName(name);
    hEtaNMixCorr->Divide(hEtaNMixEff);

    hEtaACorr=(TH1F*)hEtaARaw->Clone();
    sprintf(name,"hEtaACorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hEtaACorr->SetName(name);
    hEtaACorr->Divide(hEtaAEff);
    
    hEtaAMixCorr=(TH1F*)hEtaAMixRaw->Clone();
    sprintf(name,"hEtaAMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hEtaAMixCorr->SetName(name);
    hEtaAMixCorr->Divide(hEtaAMixEff);
    
    hPhiEtaCorr=(TH2F*)hPhiEtaRaw->Clone();
    sprintf(name,"hPhiEtaCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hPhiEtaCorr->SetName(name);
    hPhiEtaCorr->Divide(hPhiEtaEff);
    
    hPhiEtaMixCorr=(TH2F*)hPhiEtaMixRaw->Clone();
    sprintf(name,"hPhiMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hPhiEtaMixCorr->SetName(name);
    hPhiEtaMixCorr->Divide(hPhiEtaMixEff);
    //GeoCorr(hEtaNCorr,hEtaNMixCorr,hEtaNMixMC,hEtaACorr,hEtaAMixCorr,hEtaAMixMC,hPhiEtaCorr,hPhiEtaMixCorr,hPhiEtaMixMC);
    //GeoCorr2(hEtaNCorr,hEtaNMixCorr,hEtaACorr,hEtaAMixCorr,hPhiEtaCorr,hPhiEtaMixCorr);
     MixedCorrect(hPhiCorr,hPhiMixCorr,hEtaNCorr,hEtaNMixCorr,hEtaACorr,hEtaAMixCorr,hPhiEtaCorr,hPhiEtaMixCorr);
    //Deal with efficiencies of 0
    for(int x=1;x<=hPhiCorr->GetNbinsX();x++){
      if(hPhiEff->GetBinContent(x)<1E-10)hPhiCorr->SetBinContent(x,0);
      if(hPhiMixEff->GetBinContent(x)<1E-10)hPhiMixCorr->SetBinContent(x,0);
    }
    for(int x=1;x<=hEtaNCorr->GetNbinsX();x++){
      if(hEtaNEff->GetBinContent(x)<1E-10)hEtaNCorr->SetBinContent(x,0);
      if(hEtaNMixEff->GetBinContent(x)<1E-10)hEtaNMixCorr->SetBinContent(x,0);
      if(hEtaAEff->GetBinContent(x)<1E-10)hEtaACorr->SetBinContent(x,0);
      if(hEtaAMixEff->GetBinContent(x)<1E-10)hEtaAMixCorr->SetBinContent(x,0);
    }
    for(int x=1;x<=hPhiEtaCorr->GetNbinsX();x++){ 
      for(int y=1;y<=hPhiCorr->GetNbinsX();y++){
	if(hPhiEtaEff->GetBinContent(x,y)<1E-10)hPhiEtaCorr->SetBinContent(x,y,0);
	if(hPhiEtaMixEff->GetBinContent(x,y)<1E-10)hPhiEtaMixCorr->SetBinContent(x,y,0);
      }
    }
    
  
    if(!notDelta&&!ipt){
      float a[3];
      if(ZYACent>0.0001)ZYA1(hPhiCorr,hPhiMixCorr,a,ZYACent,ZYAWidth);
      // cout << "a: "  << a[0] << "  Error: " << a[1] << "  Error Code: " << a[2] << "  <pT>: " << MPt[0]  << endl;
      else if(ZYACent==0){
	ZYAM2D(hPhiEtaCorr,hPhiEtaMixCorr,a,80,10);
      }
      else
	ZYAM2D2(hPhiEtaCorr,hPhiEtaMixCorr,a,fabs(ZYACent),ZYAWidth);
    }

    if(ptWeighted==2&&ipt==0){
      TH1F *TempPhiRaw=(TH1F*)hPhiRaw->Clone();
      TempPhiRaw->SetName("TempPhiRaw");
      TH1F *TempPhiMixRaw=(TH1F*)hPhiMixRaw->Clone();
      TempPhiMixRaw->SetName("TempPhiMixRaw");
      TH1F *TempPhiCorr=(TH1F*)hPhiCorr->Clone();
      TempPhiCorr->SetName("TempPhiCorr");
      TH1F *TempPhiMixCorr=(TH1F*)hPhiMixCorr->Clone();
      TempPhiMixCorr->SetName("TemPhiMixCorr");
      TH1F *TempPhiCorr2=(TH1F*)hPhiCorr->Clone();
      TempPhiCorr2->SetName("TempPhiCorr2");
      TH1F *TempPhiMixCorr2=(TH1F*)hPhiMixCorr->Clone();
      TempPhiMixCorr2->SetName("TemPhiMixCorr2");
      
      TH1F *TempEtaNRaw=(TH1F*)hEtaNRaw->Clone();
      TempEtaNRaw->SetName("TempEtaNRaw");
      TH1F *TempEtaNMixRaw=(TH1F*)hEtaNMixRaw->Clone();
      TempEtaNMixRaw->SetName("TempEtaNMixRaw");
      TH1F *TempEtaNCorr=(TH1F*)hEtaNCorr->Clone();
      TempEtaNCorr->SetName("TempEtaNCorr");
      TH1F *TempEtaNMixCorr=(TH1F*)hEtaNMixCorr->Clone();
      TempEtaNMixCorr->SetName("TemEtaNMixCorr");
      TH1F *TempEtaNCorr2=(TH1F*)hEtaNCorr->Clone();
      TempEtaNCorr2->SetName("TempEtaNCorr2");
      TH1F *TempEtaNMixCorr2=(TH1F*)hEtaNMixCorr->Clone();
      TempEtaNMixCorr2->SetName("TemEtaNMixCorr2");

      TH1F *TempEtaARaw=(TH1F*)hEtaARaw->Clone();
      TempEtaARaw->SetName("TempEtaARaw");
      TH1F *TempEtaAMixRaw=(TH1F*)hEtaAMixRaw->Clone();
      TempEtaAMixRaw->SetName("TempEtaAMixRaw");
      TH1F *TempEtaACorr=(TH1F*)hEtaACorr->Clone();
      TempEtaACorr->SetName("TempEtaACorr");
      TH1F *TempEtaAMixCorr=(TH1F*)hEtaAMixCorr->Clone();
      TempEtaAMixCorr->SetName("TemEtaAMixCorr");
      TH1F *TempEtaACorr2=(TH1F*)hEtaACorr->Clone();
      TempEtaACorr2->SetName("TempEtaAMixCorr2");
      TH1F *TempEtaAMixCorr2=(TH1F*)hEtaAMixCorr->Clone();
      TempEtaAMixCorr2->SetName("TemEtaAMixCorr2");

      TH1F *TempPhiEtaRaw=(TH1F*)hPhiEtaRaw->Clone();
      TempPhiEtaRaw->SetName("TempPhiEtaRaw");
      TH1F *TempPhiEtaMixRaw=(TH1F*)hPhiEtaMixRaw->Clone();
      TempPhiEtaMixRaw->SetName("TempPhiEtaMixRaw");
      TH1F *TempPhiEtaCorr=(TH1F*)hPhiEtaCorr->Clone();
      TempPhiEtaCorr->SetName("TempPhiEtaCorr");
      TH1F *TempPhiEtaMixCorr=(TH1F*)hPhiEtaMixCorr->Clone();
      TempPhiEtaMixCorr->SetName("TemPhiEtaMixCorr");
      TH1F *TempPhiEtaCorr2=(TH1F*)hPhiEtaCorr->Clone();
      TempPhiEtaCorr2->SetName("TempPhiEtaCorr2");
      TH1F *TempPhiEtaMixCorr2=(TH1F*)hPhiEtaMixCorr->Clone();
      TempPhiEtaMixCorr2->SetName("TempPhiEtaMixCorr2");

    }
    }//ptWeighted loop
    
    
    TH1F *hPhiCorr2=(TH1F*)hPhiCorr->Clone();
    hPhiCorr2->SetName("hPhiCorr2");
    sprintf(name,"Bg Sub. %s",hPhiCorr2->GetTitle());
    hPhiCorr2->SetTitle(name);
    if(noTitle)hPhiCorr2->SetTitle("");
    
    TH1F *hPhiMixCorr2=(TH1F*)hPhiMixCorr->Clone();
    hPhiMixCorr2->SetName("hPhiMixCorr2");
    sprintf(name,"Normalized %1.3f %s",a[0],hPhiMixCorr2->GetTitle());
    hPhiMixCorr2->SetTitle(name);
    if(noTitle)hPhiMixCorr2->SetTitle("");
    hPhiMixCorr2->SetMarkerStyle(25);
    hPhiMixCorr2->SetMarkerColor(4);
    hPhiMixCorr2->SetLineColor(4);

    TH1F *hEtaNCorr2=(TH1F*)hEtaNCorr->Clone();
    hEtaNCorr2->SetName("hEtaNCorr2");
    sprintf(name,"Bg Sub. %s",hEtaNCorr2->GetTitle());
    hEtaNCorr2->SetTitle(name);
    if(noTitle)hEtaNCorr2->SetTitle("");
    
    TH1F *hEtaNMixCorr2=(TH1F*)hEtaNMixCorr->Clone();
    hEtaNMixCorr2->SetName("hEtaNMixCorr2");
    sprintf(name,"Normalized %1.3f %s",a[0],hEtaNMixCorr2->GetTitle());
    hEtaNMixCorr2->SetTitle(name);
    if(noTitle)hEtaNMixCorr2->SetTitle("");
    hEtaNMixCorr2->SetMarkerStyle(25);
    hEtaNMixCorr2->SetMarkerColor(4);
    hEtaNMixCorr2->SetLineColor(4);

    TH1F *hEtaACorr2=(TH1F*)hEtaACorr->Clone();
    hEtaACorr2->SetName("hEtaACorr2");
    sprintf(name,"Bg Sub. %s",hEtaACorr2->GetTitle());
    hEtaACorr2->SetTitle(name);
    if(noTitle)hEtaACorr2->SetTitle("");
    
    TH1F *hEtaAMixCorr2=(TH1F*)hEtaAMixCorr->Clone();
    hEtaAMixCorr2->SetName("hEtaAMixCorr2");
    sprintf(name,"Normalized %1.3f %s",a[0],hEtaAMixCorr2->GetTitle());
    hEtaAMixCorr2->SetTitle(name);
    if(noTitle)hEtaAMixCorr2->SetTitle("");
    hEtaAMixCorr2->SetMarkerStyle(25);
    hEtaAMixCorr2->SetMarkerColor(4);
    hEtaAMixCorr2->SetLineColor(4);
    
    TH2F *hPhiEtaCorr2=(TH2F*)hPhiEtaCorr->Clone();
    hPhiEtaCorr2->SetName("hPhiEtaCorr2");
    sprintf(name,"BG Sub. %s",hPhiEtaCorr2->GetTitle());
    hPhiEtaCorr2->SetTitle(name);
    if(noTitle)hPhiEtaCorr2->SetTitle("");
    
    TH2F *hPhiEtaMixCorr2=(TH2F*)hPhiEtaMixCorr->Clone();
    hPhiEtaMixCorr2->SetName("hPhiEtaMixCorr2");
    sprintf(name,"Normalized %1.3f %s",a[0],hPhiEtaMixCorr2->GetTitle());
    hPhiEtaMixCorr2->SetTitle(name);
    if(noTitle)hPhiEtaMixCorr2->SetTitle("");
    if(!notDelta){
      hPhiMixCorr2->Scale(a[0]);
      hEtaNMixCorr2->Scale(a[0]);
      hEtaAMixCorr2->Scale(a[0]);
      hPhiEtaMixCorr2->Scale(a[0]);
      
      //	   hPhiCorr2->Add(hPhiMixCorr2,-1);
      Add1D(hPhiCorr2,hPhiMixCorr2,-1);
      Add1D(hEtaNCorr2,hEtaNMixCorr2,-1);
      Add1D(hEtaACorr2,hEtaAMixCorr2,-1);
      Add2D(hPhiEtaCorr2,hPhiEtaMixCorr2,-1); 
      if(ptWeighted==2){
	TempPhiMixCorr2->Scale(a[0]);
	TempEtaNMixCorr2->Scale(a[0]);
	TempEtaAMixCorr2->Scale(a[0]);
	TempPhiEtaMixCorr2->Scale(a[0]);
	Add1D(TempPhiCorr2,TempPhiMixCorr2,-1);
	Add1D(TempEtaNCorr2,TempEtaNMixCorr2,-1);
	Add1D(TempEtaACorr2,TempEtaAMixCorr2,-1);
	//cout << TempPhiEtaCorr2->GetBinContent(1,1) << " " << TempPhiEtaMixCorr2->GetBinContent(1,1) << endl;
	//Add2D(TempPhiEtaCorr2,TempPhiEtaMixCorr2,-1); 	
	TempPhiEtaCorr2->Add(TempPhiEtaMixCorr2,-1);


	/*
	hPhiRaw->Divide(TempPhiRaw);
	hPhiCorr->Divide(TempPhiCorr);

	hPhiMixRaw->Divide(TempPhiMixRaw);
	hPhiMixCorr->Divide(TempPhiMixCorr);
	hPhiMixCorr2->Divide(TemPhiMixCorr2);
	*/
       	hPhiCorr2->Divide(TempPhiCorr2);
	DivideMean1D(hPhiRaw,TempPhiRaw,APt1,APt2,MPt[2]);
	DivideMean1DCorr(hPhiCorr,TempPhiCorr,TempPhiRaw,APt1,APt2,MPt[2]);
	//DivideMean1DBgSub(hPhiCorr2,TempPhiCorr2,TempPhiRaw,TempPhiMixRaw,APt1,APt2,MPt[2],MPt[3]);
	DivideMean1D(hPhiMixRaw,TempPhiMixRaw,APt1,APt2,MPt[3]);
	DivideMean1DCorr(hPhiMixCorr,TempPhiMixCorr,TempPhiMixRaw,APt1,APt2,MPt[3]);
	DivideMean1DCorr(hPhiMixCorr2,TempPhiMixCorr2,TempPhiMixRaw,APt1,APt2,MPt[3]);

	/*
	hEtaNRaw->Divide(TempEtaNRaw);
	hEtaNCorr->Divide(TempEtaNCorr);

	hEtaNMixRaw->Divide(TempEtaNMixRaw);
	hEtaNMixCorr->Divide(TempEtaNMixCorr);
	hEtaNMixCorr2->Divide(TempEtaNMixCorr2);

	hEtaARaw->Divide(TempEtaARaw);
	hEtaACorr->Divide(TempEtaACorr);

	hEtaAMixRaw->Divide(TempEtaAMixRaw);
	hEtaAMixCorr->Divide(TempEtaAMixCorr);
	hEtaAMixCorr2->Divide(TempEtaAMixCorr2);
	*/
	hEtaNCorr2->Divide(TempEtaNCorr2);
	DivideMean1D(hEtaNRaw,TempEtaNRaw,APt1,APt2,MPt[2]);
	DivideMean1DCorr(hEtaNCorr,TempEtaNCorr,TempEtaNRaw,APt1,APt2,MPt[2]);
	//DivideMean1DBgSub(hEtaNCorr2,TempEtaNCorr2,TempEtaNRaw,TempEtaNMixRaw,APt1,APt2,MPt[2],MPt[3]);
	DivideMean1D(hEtaNMixRaw,TempEtaNMixRaw,APt1,APt2,MPt[3]);
	DivideMean1DCorr(hEtaNMixCorr,TempEtaNMixCorr,TempEtaNMixRaw,APt1,APt2,MPt[3]);
	DivideMean1DCorr(hEtaNMixCorr2,TempEtaNMixCorr2,TempEtaNMixRaw,APt1,APt2,MPt[3]);

	hEtaACorr2->Divide(TempEtaACorr2);
	DivideMean1D(hEtaARaw,TempEtaARaw,APt1,APt2,MPt[2]);
	DivideMean1DCorr(hEtaACorr,TempEtaACorr,TempEtaARaw,APt1,APt2,MPt[2]);
	//DivideMean1DBgSub(hEtaACorr2,TempEtaACorr2,TempEtaARaw,TempEtaAMixRaw,APt1,APt2,MPt[2],MPt[3]);
	DivideMean1D(hEtaAMixRaw,TempEtaAMixRaw,APt1,APt2,MPt[3]);
	DivideMean1DCorr(hEtaAMixCorr,TempEtaAMixCorr,TempEtaAMixRaw,APt1,APt2,MPt[3]);
	DivideMean1DCorr(hEtaAMixCorr2,TempEtaAMixCorr2,TempEtaAMixRaw,APt1,APt2,MPt[3]);
	/*
	//Its not liking my 2D functions for some reason
	DivideMean2D(hPhiEtaRaw,TempPhiEtaRaw,APt1,APt2,MPt[2]);
	DivideMean2D(hPhiEtaCorr,TempPhiEtaCorr,APt1,APt2,MPt[2]);
	DivideMean2D(hPhiEtaCorr2,TempPhiEtaCorr2,APt1,APt2,MPt[2]);
	DivideMean2D(hPhiEtaMixRaw,TempPhiEtaMixRaw,APt1,APt2,MPt[3]);
	DivideMean2D(hPhiEtaMixCorr,TempPhiEtaMixCorr,APt1,APt2,MPt[3]);
	DivideMean2D(hPhiEtaMixCorr2,TempPhiEtaMixCorr2,APt1,APt2,MPt[3]);	
	*/
	//Errors are messed up but its not liking my 2d divide
	hPhiEtaRaw->Divide(TempPhiEtaRaw);
	hPhiEtaCorr->Divide(TempPhiEtaCorr);
	hPhiEtaCorr2->Divide(TempPhiEtaCorr2);
	hPhiEtaMixRaw->Divide(TempPhiEtaMixRaw);
	hPhiEtaMixCorr->Divide(TempPhiEtaMixCorr);
	hPhiEtaMixCorr2->Divide(TempPhiEtaMixCorr2);
	

	hPhiRaw->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiCorr->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiCorr2->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiMixRaw->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiMixCorr->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiMixCorr2->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");

	hEtaNRaw->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaNCorr->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaNCorr2->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaNMixRaw->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaNMixCorr->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaNMixCorr2->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");

	hEtaARaw->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaACorr->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaACorr2->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaAMixRaw->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaAMixCorr->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");
	hEtaAMixCorr2->GetYaxis()->SetTitle("<p_{T}^{Assoc}>");

	hPhiEtaRaw->GetZaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiEtaCorr->GetZaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiEtaCorr2->GetZaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiEtaMixRaw->GetZaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiEtaMixCorr->GetZaxis()->SetTitle("<p_{T}^{Assoc}>");
	hPhiEtaMixCorr2->GetZaxis()->SetTitle("<p_{T}^{Assoc}>");
      }
    }	
    
    if(Cent==Cent1){
      TH1F *hPhiRawPlot=(TH1F*)hPhiRaw->Clone();
      hPhiRawPlot->SetName("hPhiRawPlot");
      TH1F *hPhiCorrPlot=(TH1F*)hPhiCorr->Clone();
      hPhiCorrPlot->SetName("hPhiCorrPlot");
      TH1F *hPhiCorr2Plot=(TH1F*)hPhiCorr2->Clone();
      hPhiCorr2Plot->SetName("hPhiCorr2Plot");
      TH1F *hPhiEffPlot=(TH1F*)hPhiEff->Clone();
      hPhiEffPlot->SetName("hPhiEffPlot");
      TH1F *hPhiMCPlot=(TH1F*)hPhiMC->Clone();
      hPhiMCPlot->SetName("hPhiMCPlot");
      TH1F *hPhiMixRawPlot=(TH1F*)hPhiMixRaw->Clone();
      hPhiMixRawPlot->SetName("hPhiMixRawPlot");
      TH1F *hPhiMixCorrPlot=(TH1F*)hPhiMixCorr->Clone();
      hPhiMixCorrPlot->SetName("hPhiMixCorrPlot");
      TH1F *hPhiMixCorr2Plot=(TH1F*)hPhiMixCorr2->Clone();
      hPhiMixCorr2Plot->SetName("hPhiMixCorr2Plot");
      TH1F *hPhiMixEffPlot=(TH1F*)hPhiMixEff->Clone();
      hPhiMixEffPlot->SetName("hPhiMixEffPlot");
      TH1F *hPhiMixMCPlot=(TH1F*)hPhiMixMC->Clone();
      hPhiMixMCPlot->SetName("hPHiMixMCPlot");
      

      TH1F *hEtaNRawPlot=(TH1F*)hEtaNRaw->Clone();
      hEtaNRawPlot->SetName("hEtaNRawPlot");
      TH1F *hEtaNCorrPlot=(TH1F*)hEtaNCorr->Clone();
      hEtaNCorrPlot->SetName("hEtaNCorrPlot");
      TH1F *hEtaNCorr2Plot=(TH1F*)hEtaNCorr2->Clone();
      hEtaNCorr2Plot->SetName("hEtaNCorr2Plot");
      TH1F *hEtaNEffPlot=(TH1F*)hEtaNEff->Clone();
      hEtaNEffPlot->SetName("hEtaNEffPlot");
      TH1F *hEtaNMixRawPlot=(TH1F*)hEtaNMixRaw->Clone();
      hEtaNMixRawPlot->SetName("hEtaNMixRawPlot");
      TH1F *hEtaNMixCorrPlot=(TH1F*)hEtaNMixCorr->Clone();
      hEtaNMixCorrPlot->SetName("hEtaNMixCorrPlot");
      TH1F *hEtaNMixCorr2Plot=(TH1F*)hEtaNMixCorr2->Clone();
      hEtaNMixCorr2Plot->SetName("hEtaNMixCorr2Plot");
      TH1F *hEtaNMixEffPlot=(TH1F*)hEtaNMixEff->Clone();
      hEtaNMixEffPlot->SetName("hEtaNMixEffPlot");

      TH1F *hEtaARawPlot=(TH1F*)hEtaARaw->Clone();
      hEtaARawPlot->SetName("hEtaARawPlot");
      TH1F *hEtaACorrPlot=(TH1F*)hEtaACorr->Clone();
      hEtaACorrPlot->SetName("hEtaACorrPlot");
      TH1F *hEtaACorr2Plot=(TH1F*)hEtaACorr2->Clone();
      hEtaACorr2Plot->SetName("hEtaACorr2Plot");
      TH1F *hEtaAEffPlot=(TH1F*)hEtaAEff->Clone();
      hEtaAEffPlot->SetName("hEtaAEffPlot");
      TH1F *hEtaAMixRawPlot=(TH1F*)hEtaAMixRaw->Clone();
      hEtaAMixRawPlot->SetName("hEtaAMixRawPlot");
      TH1F *hEtaAMixCorrPlot=(TH1F*)hEtaAMixCorr->Clone();
      hEtaAMixCorrPlot->SetName("hEtaAMixCorrPlot");
      TH1F *hEtaAMixCorr2Plot=(TH1F*)hEtaAMixCorr2->Clone();
      hEtaAMixCorr2Plot->SetName("hEtaAMixCorr2Plot");
      TH1F *hEtaAMixEffPlot=(TH1F*)hEtaAMixEff->Clone();
      hEtaAMixEffPlot->SetName("hEtaAMixEffPlot");

      TH2F *hPhiEtaRawPlot=(TH2F*)hPhiEtaRaw->Clone();
      hPhiEtaRawPlot->SetName("hPhiEtaRawPlot");
      TH2F *hPhiEtaCorrPlot=(TH2F*)hPhiEtaCorr->Clone();
      hPhiEtaCorrPlot->SetName("hPhiEtaCorrPlot");
      TH2F *hPhiEtaCorr2Plot=(TH2F*)hPhiEtaCorr2->Clone();
      hPhiEtaCorr2Plot->SetName("hPhiEtaCorr2Plot");
      TH2F *hPhiEtaEffPlot=(TH2F*)hPhiEtaEff->Clone();
      hPhiEtaEffPlot->SetName("hPhiEtaEffPlot");
      TH2F *hPhiEtaMixRawPlot=(TH2F*)hPhiEtaMixRaw->Clone();
      hPhiEtaMixRawPlot->SetName("hPhiEtaMixRawPlot");
      TH2F *hPhiEtaMixCorrPlot=(TH2F*)hPhiEtaMixCorr->Clone();
      hPhiEtaMixCorrPlot->SetName("hPhiEtaMixRawPlot");
      TH2F *hPhiEtaMixCorr2Plot=(TH2F*)hPhiEtaMixCorr2->Clone();
      hPhiEtaMixCorr2Plot->SetName("hPhiEtaMixCorr2Plot");
      TH2F *hPhiEtaMixEffPlot=(TH2F*)hPhiEtaMixEff->Clone();
      hPhiEtaMixEffPlot->SetName("hPhiEtaMixEffPlot");
      
      TrigSum=MPt[2];
      TrigSum2=MPt2[2];

      hPhiRawPlot->Scale(MPt[2]);
      hPhiCorrPlot->Scale(MPt[2]);
      hPhiCorr2Plot->Scale(MPt[2]);
      hPhiEffPlot->Scale(MPt[2]);
      hPhiMCPlot->Scale(MPt2[2]);
      hPhiMixRawPlot->Scale(MPt[2]);
      hPhiMixCorrPlot->Scale(MPt[2]);
      hPhiMixCorr2Plot->Scale(MPt[2]);
      hPhiMixEffPlot->Scale(MPt[2]);
      hPhiMixMCPlot->Scale(MPt2[2]);

      hEtaNRawPlot->Scale(MPt[2]);
      hEtaNCorrPlot->Scale(MPt[2]);
      hEtaNCorr2Plot->Scale(MPt[2]);
      hEtaNEffPlot->Scale(MPt[2]);
      hEtaNMixRawPlot->Scale(MPt[2]);
      hEtaNMixCorrPlot->Scale(MPt[2]);
      hEtaNMixCorr2Plot->Scale(MPt[2]);
      hEtaNMixEffPlot->Scale(MPt[2]);

      hEtaARawPlot->Scale(MPt[2]);
      hEtaACorrPlot->Scale(MPt[2]);
      hEtaACorr2Plot->Scale(MPt[2]);
      hEtaAEffPlot->Scale(MPt[2]);
      hEtaAMixRawPlot->Scale(MPt[2]);
      hEtaAMixCorrPlot->Scale(MPt[2]);
      hEtaAMixCorr2Plot->Scale(MPt[2]);
      hEtaAMixEffPlot->Scale(MPt[2]);
      
      hPhiEtaRawPlot->Scale(MPt[2]);
      hPhiEtaCorrPlot->Scale(MPt[2]);
      hPhiEtaCorr2Plot->Scale(MPt[2]);
      hPhiEtaEffPlot->Scale(MPt[2]);
      hPhiEtaMixRawPlot->Scale(MPt[2]);
      hPhiEtaMixCorrPlot->Scale(MPt[2]);
      hPhiEtaMixCorr2Plot->Scale(MPt[2]);
      hPhiEtaMixEffPlot->Scale(MPt[2]);
    }
    else{//I need to do trigger weighted averaging here so I need the trigger counts
      hPhiRawPlot->Add(hPhiRaw,MPt[2]);
      hPhiCorrPlot->Add(hPhiCorr,MPt[2]);
      hPhiCorr2Plot->Add(hPhiCorr2,MPt[2]);
      hPhiEffPlot->Add(hPhiEff,MPt[2]);
      hPhiMCPlot->Add(hPhiMC,MPt2[2]);
      hPhiMixRawPlot->Add(hPhiMixRaw,MPt[2]);
      hPhiMixCorrPlot->Add(hPhiMixCorr,MPt[2]);
      hPhiMixCorr2Plot->Add(hPhiMixCorr2,MPt[2]);
      hPhiMixEffPlot->Add(hPhiMixEff,MPt[2]);
      hPhiMCPlot->Add(hPhiMixMC,MPt2[2]);

      hEtaNRawPlot->Add(hEtaNRaw,MPt[2]);
      hEtaNCorrPlot->Add(hEtaNCorr,MPt[2]);
      hEtaNCorr2Plot->Add(hEtaNCorr2,MPt[2]);
      hEtaNEffPlot->Add(hEtaNEff,MPt[2]);
      hEtaNMixRawPlot->Add(hEtaNMixRaw,MPt[2]);
      hEtaNMixCorrPlot->Add(hEtaNMixCorr,MPt[2]);
      hEtaNMixCorr2Plot->Add(hEtaNMixCorr2,MPt[2]);
      hEtaNMixEffPlot->Add(hEtaNMixEff,MPt[2]);

      hEtaARawPlot->Add(hEtaARaw,MPt[2]);
      hEtaACorrPlot->Add(hEtaACorr,MPt[2]);
      hEtaACorr2Plot->Add(hEtaACorr2,MPt[2]);
      hEtaAEffPlot->Add(hEtaAEff,MPt[2]);
      hEtaAMixRawPlot->Add(hEtaAMixRaw,MPt[2]);
      hEtaAMixCorrPlot->Add(hEtaAMixCorr,MPt[2]);
      hEtaAMixCorr2Plot->Add(hEtaAMixCorr2,MPt[2]);
      hEtaAMixEffPlot->Add(hEtaAMixEff,MPt[2]);

      hPhiEtaRawPlot->Add(hPhiEtaRaw,MPt[2]);
      hPhiEtaCorrPlot->Add(hPhiEtaCorr,MPt[2]);
      hPhiEtaCorr2Plot->Add(hPhiEtaCorr2,MPt[2]);
      hPhiEtaEffPlot->Add(hPhiEtaEff,MPt[2]);
      hPhiEtaMixRawPlot->Add(hPhiEtaMixRaw,MPt[2]);
      hPhiEtaMixCorrPlot->Add(hPhiEtaCorr,MPt[2]);
      hPhiEtaMixCorr2Plot->Add(hPhiEtaMixCorr2,MPt[2]);
      hPhiEtaMixEffPlot->Add(hPhiEtaMixEff,MPt[2]);
      TrigSum+=MPt[2];
      TrigSum2+=MPt2[2];
    }
  }
  if(TrigSum==0)cout << "TrigSum: " << TrigSum << endl;
  hPhiRawPlot->Scale(1./TrigSum);
  hPhiCorrPlot->Scale(1./TrigSum);
  hPhiCorr2Plot->Scale(1./TrigSum);
  hPhiEffPlot->Scale(1./TrigSum);
  hPhiMCPlot->Scale(1./TrigSum2);
  hPhiMixRawPlot->Scale(1./TrigSum);
  hPhiMixCorrPlot->Scale(1./TrigSum);
  hPhiMixCorr2Plot->Scale(1./TrigSum);
  hPhiMixEffPlot->Scale(1./TrigSum);
  hPhiMixMCPlot->Scale(1./TrigSum);

  hEtaNRawPlot->Scale(1./TrigSum);
  hEtaNCorrPlot->Scale(1./TrigSum);
  hEtaNCorr2Plot->Scale(1./TrigSum);
  hEtaNEffPlot->Scale(1./TrigSum);
  hEtaNMixRawPlot->Scale(1./TrigSum);
  hEtaNMixCorrPlot->Scale(1./TrigSum);
  hEtaNMixCorr2Plot->Scale(1./TrigSum);
  hEtaNMixEffPlot->Scale(1./TrigSum);

  hEtaARawPlot->Scale(1./TrigSum);
  hEtaACorrPlot->Scale(1./TrigSum);
  hEtaACorr2Plot->Scale(1./TrigSum);
  hEtaAEffPlot->Scale(1./TrigSum);
  hEtaAMixRawPlot->Scale(1./TrigSum);
  hEtaAMixCorrPlot->Scale(1./TrigSum);
  hEtaAMixCorr2Plot->Scale(1./TrigSum);
  hEtaAMixEffPlot->Scale(1./TrigSum);

  hPhiEtaRawPlot->Scale(1./TrigSum);
  hPhiEtaCorrPlot->Scale(1./TrigSum);
  hPhiEtaCorr2Plot->Scale(1./TrigSum);
  hPhiEtaEffPlot->Scale(1./TrigSum);
  hPhiEtaMixRawPlot->Scale(1./TrigSum);
  hPhiEtaMixCorrPlot->Scale(1./TrigSum);
  hPhiEtaMixCorr2Plot->Scale(1./TrigSum);
  hPhiEtaMixEffPlot->Scale(1./TrigSum); 
  
  if(noTitle){
    hPhiRawPlot->SetTitle("");
    hPhiCorrPlot->SetTitle("");
    hPhiCorr2Plot->SetTitle("");
    hPhiEffPlot->SetTitle("");
    hPhiMCPlot->SetTitle("");
    hPhiMixRawPlot->SetTitle("");
    hPhiMixCorrPlot->SetTitle("");
    hPhiMixCorr2Plot->SetTitle("");
    hPhiMixEffPlot->SetTitle("");
    hPhiMixMCPlot->SetTitle("");

    hEtaNRawPlot->SetTitle("");
    hEtaNCorrPlot->SetTitle("");
    hEtaNCorr2Plot->SetTitle("");
    hEtaNEffPlot->SetTitle("");
    hEtaNMixRawPlot->SetTitle("");
    hEtaNMixCorrPlot->SetTitle("");
    hEtaNMixCorr2Plot->SetTitle("");
    hEtaNMixEffPlot->SetTitle("");

    hEtaARawPlot->SetTitle("");
    hEtaACorrPlot->SetTitle("");
    hEtaACorr2Plot->SetTitle("");
    hEtaAEffPlot->SetTitle("");
    hEtaAMixRawPlot->SetTitle("");
    hEtaAMixCorrPlot->SetTitle("");
    hEtaAMixCorr2Plot->SetTitle("");
    hEtaAMixEffPlot->SetTitle("");
    
    hPhiEtaRawPlot->SetTitle("");
    hPhiEtaCorrPlot->SetTitle("");
    hPhiEtaCorr2Plot->SetTitle("");
    hPhiEtaEffPlot->SetTitle("");
    hPhiEtaMixRawPlot->SetTitle("");
    hPhiEtaMixCorrPlot->SetTitle("");
    hPhiEtaMixCorr2Plot->SetTitle("");
    hPhiEtaMixEffPlot->SetTitle(""); 
  }
  
  if(ReBinPhi>1){
    rebin1D(hPhiRawPlot,ReBinPhi);
    rebin1D(hPhiCorrPlot,ReBinPhi);
    hPhiCorr2Plot=rebin1D(hPhiCorr2Plot,ReBinPhi);
    hPhiEffPlot=rebin1D(hPhiEffPlot,ReBinPhi);
    hPhiMCPlot=rebins1D(hPhiMCPlot,ReBinPhi);
    hPhiMixRawPlot=rebin1D(hPhiMixRawPlot,ReBinPhi);
    hPhiMixCorrPlot=rebin1D(hPhiMixCorrPlot,ReBinPhi);
    hPhiMixCorr2Plot=rebin1D(hPhiMixCorr2Plot,ReBinPhi);
    hPhiMixEffPlot=rebin1D(hPhiMixEffPlot,ReBinPhi);
    hPhiMixMCPlot=rebin1D(hPhiMixMCPlot,ReBinPhi);
  }
 if(ReBinEta>1){
    rebin1D(hEtaNRawPlot,ReBinEtaN);
    rebin1D(hEtaNCorrPlot,ReBinEtaN);
    hEtaNCorr2Plot=rebin1D(hEtaNCorr2Plot,ReBinEtaN);
    hEtaNEffPlot=rebin1D(hEtaNEffPlot,ReBinEtaN);
    hEtaNMixRawPlot=rebin1D(hEtaNMixRawPlot,ReBinEtaN);
    hEtaNMixCorrPlot=rebin1D(hEtaNMixCorrPlot,ReBinEtaN);
    hEtaNMixCorr2Plot=rebin1D(hEtaNMixCorr2Plot,ReBinEtaN);
    hEtaNMixEffPlot=rebin1D(hEtaNMixEffPlot,ReBinEtaN);

    rebin1D(hEtaARawPlot,ReBinEtaA);
    rebin1D(hEtaACorrPlot,ReBinEtaA);
    hEtaACorr2Plot=rebin1D(hEtaACorr2Plot,ReBinEtaA);
    hEtaAEffPlot=rebin1D(hEtaAEffPlot,ReBinEtaA);
    hEtaAMixRawPlot=rebin1D(hEtaAMixRawPlot,ReBinEtaA);
    hEtaAMixCorrPlot=rebin1D(hEtaAMixCorrPlot,ReBinEtaA);
    hEtaAMixCorr2Plot=rebin1D(hEtaAMixCorr2Plot,ReBinEtaA);
    hEtaAMixEffPlot=rebin1D(hEtaAMixEffPlot,ReBinEtaA);
  }


 ///Phi Plots////////////////////////////
  c1=new TCanvas("c1","PhiRaw");
  SetMargins1D(c1);
  // if(DrawAsHisto)gStyle->SetHistFillColor(2);
 
  // cout << hPhiRawPlot->GetXaxis()->GetTitle() << endl;
  hPhiRawPlot->Draw();
  if(DrawAsHisto)DrawLego1D(c1,hPhiRawPlot);
   
  sprintf(name,"%s/%s%sPhiRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  c1->SaveAs(name);
  sprintf(name,"%s/%s%sPhiRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS)c1->SaveAs(name);
  sprintf(name,"%s/%s%sPhiRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
 if(SavePDF)c1->SaveAs(name,"landscape");
 sprintf(name,"%s/%s%sPhiRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.C",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SaveMacro)c1->SaveAs(name);
  hPhiRawPlot->GetXaxis()->SetLabelOffset(0.005);
  sprintf(name,"%s/%s%sPhiRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.txt",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SaveText)SaveAsText(hPhiRawPlot,name);
  
  c2=new TCanvas("c2","PhiCorr");
  SetMargins1D(c2);
  hPhiCorrPlot->Draw();
  if(DrawAsHisto)DrawLego1D(c2,hPhiCorrPlot);
  sprintf(name,"%s/%s%sPhiCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  c2->SaveAs(name);
  sprintf(name,"%s/%s%sPhiCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS) c2->SaveAs(name);
  hPhiCorrPlot->GetXaxis()->SetLabelOffset(0.005);

   
 if(!notDelta){
 c2bg=new TCanvas("c2bg","PhiCorrBgSub");
 SetMargins1D(c2bg);
 if(ptWeighted==2){
   hPhiCorr2Plot->SetMaximum(APt2);
   hPhiCorr2Plot->SetMinimum(APt1);
 }
  hPhiCorr2Plot->Draw();
  ZeroLine->Draw("same");
  if(DrawAsHisto)DrawLego1D(c2bg,hPhiCorr2Plot);
  sprintf(name,"%s/%s%sPhiCorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  c2bg->SaveAs(name);
  sprintf(name,"%s/%s%sPhiCorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS)c2bg->SaveAs(name);
  hPhiCorr2Plot->GetXaxis()->SetLabelOffset(0.005);
  }
 float maxeff=1.1;
 float mineff=0.5;
  c3=new TCanvas("c3","PhiEff");
  SetMargins1D(c3);
  // hPhiEffPlot->SetMaximum(maxeff);
  // hPhiEffPlot->SetMinimum(mineff);
  hPhiEffPlot->Draw();
  sprintf(name,"%s/%s%sPhiEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  c3->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS)c3->SaveAs(name);

  c3MC=new TCanvas("c3MC","PhiMC");
  SetMargins1D(c3MC);
  hPhiMCPlot->Draw();
  sprintf(name,"%s/%s%sPhiMC_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(DrawMC)c3MC->SaveAs(name);
  sprintf(name,"%s/%s%sPhiMC_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS&&DrawMC)c3MC->SaveAs(name);

  c4=new TCanvas("c4","PhiMixRaw");
  SetMargins1D(c4);
  hPhiMixRawPlot->Draw();
  sprintf(name,"%s/%s%sPhiMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  c4->SaveAs(name);
  sprintf(name,"%s/%s%sPhiMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)c4->SaveAs(name);
  sprintf(name,"%s/%s%sPhiMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.C",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(SaveMacro)c4->SaveAs(name);
  sprintf(name,"%s/%s%sPhiMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.txt",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(SaveText)SaveAsText(hPhiMixRawPlot,name);

  c5=new TCanvas("c5","PhiMixCorr");
  SetMargins1D(c5);
  hPhiMixCorrPlot->Draw();
  sprintf(name,"%s/%s%sPhiMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  c5->SaveAs(name);
  sprintf(name,"%s/%s%sPhiMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)c5->SaveAs(name);
  
  if(!notDelta){
  c5bg=new TCanvas("c5bg","PhiMixCorrNorm");
  SetMargins1D(c5bg);
  hPhiMixCorr2Plot->Draw();
  sprintf(name,"%s/%s%sPhiMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  c5bg->SaveAs(name);
  sprintf(name,"%s/%s%sPhiMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
    if(SavePS)c5bg->SaveAs(name);
  }

  c6=new TCanvas("c6","PhiMixEff");
  SetMargins1D(c6);
  //hPhiMixEffPlot->SetMaximum(maxeff);
  //hPhiMixEffPlot->SetMinimum(mineff);
  hPhiMixEffPlot->Draw();
  sprintf(name,"%s/%s%sPhiMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  c6->SaveAs(name);
  sprintf(name,"%s/%s%sPhiMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)c6->SaveAs(name);

  c6MC=new TCanvas("c6MC","PhiMixMC");
  SetMargins1D(c6MC);
  hPhiMixMCPlot->Draw();
  sprintf(name,"%s/%s%sPhiMixMC_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(DrawMC)c6MC->SaveAs(name);
  sprintf(name,"%s/%s%sPhiMixMC_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS&&DrawMC)c6MC->SaveAs(name);

  c2Cyl=new TCanvas("c2Cyl","PhiCorrCyl");
  SetMargins1DCyl(c2Cyl);
  hPhiCorrPlot->Draw("legocyl");
  sprintf(name,"%s/%s%sPhiCorr%dCyl_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  c2Cyl->SaveAs(name);
  sprintf(name,"%s/%s%sPhiCorr%dCyl_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)c2Cyl->SaveAs(name);

  if(!notDelta){
  c2All=new TCanvas("c2All","PhiCorrAll");
  SetMargins1D(c2All);
  TH1F *hPhiCorrAll=(TH1F*)hPhiCorrPlot->Clone();
  hPhiCorrAll->SetName("hPhiCorrAll");
  hPhiCorrAll->SetMinimum(hPhiMixCorr->GetMinimum()-fabs(hPhiMixCorr2->GetMinimum()-hPhiMixCorr->GetMinimum()));
  if(hPhiCorrAll->GetMinimum()<0)hPhiCorrAll->SetMinimum(0);
  hPhiCorrAll->Draw();
  hPhiMixCorrPlot->Draw("same");
  hPhiMixCorr2Plot->Draw("same");
  hPhiCorrAll->Draw("same");
  sprintf(name,"%s/%s%sPhiCorr%dAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  c2All->SaveAs(name);
  sprintf(name,"%s/%s%sPhiCorr%dAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)c2All->SaveAs(name);
 sprintf(name,"%s/%s%sPhiCorr%dAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePDF)c2All->SaveAs(name);
}

  c2Fit=new TCanvas("c2Fit","PhiCorrFit");
  SetMargins1D(c2Fit);
  TH1F* hPhiCorrFitPlot=(TH1F*)hPhiCorrPlot->Clone();
  hPhiCorrFitPlot->SetName("hPhiCorrFitPlot");
  hPhiCorrFitPlot->Draw();
  hPhiCorrFitPlot->Fit("fit1");
  fit2->SetParameter(2,fit1->GetParameter(4));
  fit2->SetParLimits(2,0,0);
  if(DrawAsHisto)DrawLego1D(c2,hPhiCorrFitPlot);
  sprintf(name,"%s/%s%sPhiCorr%dFit_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  c2Fit->SaveAs(name);
  sprintf(name,"%s/%s%sPhiCorr%dFit_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS) c2Fit->SaveAs(name);
  sprintf(name,"%s/%s%sPhiCorr%dFit_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePDF) c2Fit->SaveAs(name);
  hPhiCorrPlot->GetXaxis()->SetLabelOffset(0.005);

  ////Near-Side Eta Plots////////////////////////////////
  cN1=new TCanvas("cN1","EtaNRaw");
  SetMargins1D(cN1); 
  hEtaNRawPlot->Draw();
  if(DrawAsHisto)DrawLego1D(cN1,hEtaNRawPlot);
  sprintf(name,"%s/%s%sEtaNRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  cN1->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS)cN1->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePDF)cN1->SaveAs(name);
   sprintf(name,"%s/%s%sEtaNRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.C",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SaveMacro)cN1->SaveAs(name);
  hEtaNRawPlot->GetXaxis()->SetLabelOffset(0.005);
  sprintf(name,"%s/%s%sEtaNRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.txt",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SaveText)SaveAsText(hEtaNRawPlot,name);
  
  cN2=new TCanvas("cN2","EtaNCorr");
  SetMargins1D(cN2);
  hEtaNCorrPlot->Draw();
  if(DrawAsHisto)DrawLego1D(cN2,hEtaNCorrPlot);
  sprintf(name,"%s/%s%sEtaNCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  cN2->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS) cN2->SaveAs(name);
 sprintf(name,"%s/%s%sEtaNCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePDF) cN2->SaveAs(name);
  hEtaNCorrPlot->GetXaxis()->SetLabelOffset(0.005);
   
 if(!notDelta){
   cN2bg=new TCanvas("cN2bg","EtaNCorrBgSub");
   SetMargins1D(cN2bg);
   if(ptWeighted==2){
     hEtaNCorr2Plot->SetMaximum(APt2);
     hEtaNCorr2Plot->SetMinimum(APt1);
   }
   hEtaNCorr2Plot->Draw();
   ZeroLine->Draw("same");
   if(DrawAsHisto)DrawLego1D(cN2bg,hEtaNCorr2Plot);
   sprintf(name,"%s/%s%sEtaNCorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
   cN2bg->SaveAs(name);
   sprintf(name,"%s/%s%sEtaNCorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
   if(SavePS)cN2bg->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNCorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
   if(SavePDF)cN2bg->SaveAs(name)
   hEtaNCorr2Plot->GetXaxis()->SetLabelOffset(0.005);
 }
 float maxeff=1.1;
 float mineff=0.5;
  cN3=new TCanvas("cN3","EtaNEff");
  SetMargins1D(cN3);
  //hEtaNEffPlot->SetMaximum(maxeff);
  //hEtaNEffPlot->SetMinimum(mineff);
  hEtaNEffPlot->Draw();
  sprintf(name,"%s/%s%sEtaNEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  cN3->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS)cN3->SaveAs(name);
 sprintf(name,"%s/%s%sEtaNEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePDF)cN3->SaveAs(name);

  cN4=new TCanvas("cN4","EtaNMixRaw");
  SetMargins1D(cN4);
  hEtaNMixRawPlot->Draw();
  sprintf(name,"%s/%s%sEtaNMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cN4->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)cN4->SaveAs(name);
sprintf(name,"%s/%s%sEtaNMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePDF)cN4->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.C",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SaveMacro)cN4->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.txt",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SaveText)SaveAsText(hEtaNMixRawPlot,name);

  cN5=new TCanvas("cN5","EtaNMixCorr");
  SetMargins1D(cN5);
  hEtaNMixCorrPlot->Draw();
  sprintf(name,"%s/%s%sEtaNMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cN5->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)cN5->SaveAs(name);
  
  if(!notDelta){
  cN5bg=new TCanvas("cN5bg","EtaNMixCorrNorm");
  SetMargins1D(cN5bg);
  hEtaNMixCorr2Plot->Draw();
  sprintf(name,"%s/%s%sEtaNMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cN5bg->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
    if(SavePS)cN5bg->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
    if(SavePDF)cN5bg->SaveAs(name);
  }

  cN6=new TCanvas("cN6","EtaNMixEff");
  SetMargins1D(cN6);
  //hEtaNMixEffPlot->SetMaximum(maxeff);
  // hEtaNMixEffPlot->SetMinimum(mineff);
  hEtaNMixEffPlot->Draw();
  sprintf(name,"%s/%s%sEtaNMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cN6->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)cN6->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePDF)cN6->SaveAs(name);

  if(!notDelta){
  cN2All=new TCanvas("cN2All","EtaNCorrAll");
  SetMargins1D(cN2All);
  TH1F *hEtaNCorrAll=(TH1F*)hEtaNCorrPlot->Clone();
  hEtaNCorrAll->SetName("hEtaNCorrAll");
  hEtaNCorrAll->SetMinimum(hEtaNMixCorr->GetMinimum()-fabs(hEtaNMixCorr2->GetMinimum()-hEtaNMixCorr->GetMinimum()));
  if(hEtaNCorrAll->GetMinimum()<0)hEtaNCorrAll->SetMinimum(0);
  hEtaNCorrAll->Draw();
  hEtaNMixCorrPlot->Draw("same");
  hEtaNMixCorr2Plot->Draw("same");
  hEtaNCorrAll->Draw("same");
  sprintf(name,"%s/%s%sEtaNCorr%dAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cN2All->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNCorr%dAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)cN2All->SaveAs(name);
 sprintf(name,"%s/%s%sEtaNCorr%dAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePDF)cN2All->SaveAs(name);
}

 cN2Fit=new TCanvas("cN2Fit","EtaNCorrFit");
  SetMargins1D(cN2Fit);
  TH1F *hEtaNCorrFitPlot=(TH1F*)hEtaNCorrPlot->Clone();
  hEtaNCorrFitPlot->SetName("hEtaNCorrFitPlot");
  hEtaNCorrFitPlot->Draw();
  hEtaNCorrFitPlot->Fit("fit2");
  if(DrawAsHisto)DrawLego1D(cN2,hEtaNCorrPlot);
  sprintf(name,"%s/%s%sEtaNCorr%dFit_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  cN2Fit->SaveAs(name);
  sprintf(name,"%s/%s%sEtaNCorr%dFit_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS) cN2Fit->SaveAs(name);
 sprintf(name,"%s/%s%sEtaNCorr%dFit_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePDF) cN2Fit->SaveAs(name);
  hEtaNCorrPlot->GetXaxis()->SetLabelOffset(0.005);

  /////Away-Side Eta Plots////////////////////////

   cA1=new TCanvas("cA1","EtaARaw");
  SetMargins1D(cA1);
  // if(DrawAsHisto)gStyle->SetHistFillColor(2);
  //cout << hEtaARawPlot->GetXaxis()->GetTitle() << endl;
  hEtaARawPlot->Draw();
  if(DrawAsHisto)DrawLego1D(cA1,hEtaARawPlot);
  sprintf(name,"%s/%s%sEtaARaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  cA1->SaveAs(name);
  sprintf(name,"%s/%s%sEtaARaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS)cA1->SaveAs(name);
 sprintf(name,"%s/%s%sEtaARaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePDF)cA1->SaveAs(name);
   sprintf(name,"%s/%s%sEtaARaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.C",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SaveMacro)cA1->SaveAs(name);
  hEtaARawPlot->GetXaxis()->SetLabelOffset(0.005);
  sprintf(name,"%s/%s%sEtaARaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.txt",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SaveText)SaveAsText(hEtaARawPlot,name);
  
  cA2=new TCanvas("cA2","EtaACorr");
  SetMargins1D(cA2);
  hEtaACorrPlot->Draw();
  if(DrawAsHisto)DrawLego1D(ca2,hEtaACorrPlot);
  sprintf(name,"%s/%s%sEtaACorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  cA2->SaveAs(name);
  sprintf(name,"%s/%s%sEtaACorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS) cA2->SaveAs(name);
  sprintf(name,"%s/%s%sEtaACorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePDF) cA2->SaveAs(name);
  hEtaACorrPlot->GetXaxis()->SetLabelOffset(0.005);
   
 if(!notDelta){
 cA2bg=new TCanvas("cA2bg","EtaACorrBgSub");
 SetMargins1D(cA2bg);
 if(ptWeighted==2){
   hEtaACorr2Plot->SetMaximum(APt2);
   hEtaACorr2Plot->SetMinimum(APt1);
 }
  hEtaACorr2Plot->Draw();
  ZeroLine->Draw("same");
  if(DrawAsHisto)DrawLego1D(c2bg,hEtaACorr2Plot);
  sprintf(name,"%s/%s%sEtaACorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  cA2bg->SaveAs(name);
  sprintf(name,"%s/%s%sEtaACorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS)cA2bg->SaveAs(name);
 sprintf(name,"%s/%s%sEtaACorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePDF)cA2bg->SaveAs(name);
  hEtaACorr2Plot->GetXaxis()->SetLabelOffset(0.005);
  }
 float maxeff=1.1;
 float mineff=0.5;
  cA3=new TCanvas("cA3","EtaAEff");
  SetMargins1D(cA3);
  //hEtaAEffPlot->SetMaximum(maxeff);
  // hEtaAEffPlot->SetMinimum(mineff);
  hEtaAEffPlot->Draw();
  sprintf(name,"%s/%s%sEtaAEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  cA3->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePS)cA3->SaveAs(name);
 sprintf(name,"%s/%s%sEtaAEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],histArray[DrawAsHisto],titArray[noTitle],ZYACent);
  if(SavePDF)cA3->SaveAs(name);

  cA4=new TCanvas("cA4","EtaAMixRaw");
  SetMargins1D(cA4);
  hEtaAMixRawPlot->Draw();
  sprintf(name,"%s/%s%sEtaAMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cA4->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)cA4->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePDF)cA4->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.C",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SaveMacro)cA4->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.txt",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SaveText)SaveAsText(hEtaAMixRawPlot,name);

  cA5=new TCanvas("cA5","EtaAMixCorr");
  SetMargins1D(cA5);
  hEtaAMixCorrPlot->Draw();
  sprintf(name,"%s/%s%sEtaAMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cA5->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)cA5->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePDF)cA5->SaveAs(name);
  
  if(!notDelta){
  cA5bg=new TCanvas("cA5bg","EtaAMixCorrNorm");
  SetMargins1D(cA5bg);
  hEtaAMixCorr2Plot->Draw();
  sprintf(name,"%s/%s%sEtaAMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cA5bg->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
    if(SavePS)cA5bg->SaveAs(name);
sprintf(name,"%s/%s%sEtaAMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
    if(SavePDF)cA5bg->SaveAs(name);
  }

  cA6=new TCanvas("cA6","EtaAMixEff");
  SetMargins1D(cA6);
  //hEtaAMixEffPlot->SetMaximum(maxeff);
  //hEtaAMixEffPlot->SetMinimum(mineff);
  hEtaAMixEffPlot->Draw();
  sprintf(name,"%s/%s%sEtaAMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cA6->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)cA6->SaveAs(name);
  sprintf(name,"%s/%s%sEtaAMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePDF)cA6->SaveAs(name);

  if(!notDelta){
  cA2All=new TCanvas("cA2All","EtaACorrAll");
  SetMargins1D(cA2All);
  TH1F *hEtaACorrAll=(TH1F*)hEtaACorrPlot->Clone();
  hEtaACorrAll->SetName("hEtaACorrAll");
  hEtaACorrAll->SetMinimum(hEtaAMixCorr->GetMinimum()-fabs(hEtaAMixCorr2->GetMinimum()-hEtaAMixCorr->GetMinimum()));
  if(hEtaACorrAll->GetMinimum()<0)hEtaACorrAll->SetMinimum(0);
  hEtaACorrAll->Draw();
  hEtaAMixCorrPlot->Draw("same");
  hEtaAMixCorr2Plot->Draw("same");
  hEtaACorrAll->Draw("same");
  sprintf(name,"%s/%s%sEtaACorr%dAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  cA2All->SaveAs(name);
  sprintf(name,"%s/%s%sEtaACorr%dAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePS)cA2All->SaveAs(name);
  sprintf(name,"%s/%s%sEtaACorr%dAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d%s%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
  if(SavePDF)cA2All->SaveAs(name);
}


 
  ////PhiEtaPlots/////////////////////////////
  c10=new TCanvas("c10","PhiEtaRaw");
  SetMargins2DSurf(c10);
  hPhiEtaRawPlot->Draw("lego2");
  sprintf(name,"%s/%s%sPhiEtaRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c10->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c10->SaveAs(name);
 sprintf(name,"%s/%s%sPhiEtaRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c10->SaveAs(name);


  c20=new TCanvas("c20","PhiEtaCorr");
  SetMargins2DSurf(c20);
  hPhiEtaCorrPlot->Draw("lego2");
  sprintf(name,"%s/%s%sPhiEtaCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c20->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c20->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c20->SaveAs(name);
  
  if(!notDelta){
  c20bg=new TCanvas("c20bg","PhiEtaCorrBgSub");
  SetMargins2DSurf(c20bg);
  if(ptWeighted==2){
    hPhiEtaCorr2Plot->SetMaximum(APt2);
    hPhiEtaCorr2Plot->SetMinimum(APt1);
  }
  hPhiEtaCorr2Plot->Draw("lego2");
  sprintf(name,"%s/%s%sPhiEtaCorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c20bg->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaCorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c20bg->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaCorr%dBgSub_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c20bg->SaveAs(name);
  }

  c30=new TCanvas("c30","PhiEtaEff");
  //hPhiEtaEffPlot->SetMaximum(maxeff);
  //hPhiEtaEffPlot->SetMinimum(mineff);
  SetMargins2DSurf(c30);
  hPhiEtaEffPlot->Draw("lego2");
  sprintf(name,"%s/%s%sPhiEtaEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c30->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c30->SaveAs(name);
 sprintf(name,"%s/%s%sPhiEtaEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c30->SaveAs(name);

  c40=new TCanvas("c40","PhiEtaMixRaw");
  SetMargins2DSurf(c40);
  hPhiEtaMixRawPlot->Draw("lego2");
  sprintf(name,"%s/%s%sPhiEtaMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c40->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c40->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixRaw_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c40->SaveAs(name);

  c50=new TCanvas("c50","PhiEtaMixCorr");
  SetMargins2DSurf(c50);
  hPhiEtaMixCorrPlot->Draw("lego2");
  sprintf(name,"%s/%s%sPhiEtaMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c50->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c50->SaveAs(name);
 sprintf(name,"%s/%s%sPhiEtaMixCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c50->SaveAs(name);
  
  if(!notDelta){
  c50bg=new TCanvas("c50bg","PhiEtaMixCorr");
  SetMargins2DSurf(c50bg);
  hPhiEtaMixCorr2Plot->Draw("lego2");
  sprintf(name,"%s/%s%sPhiEtaMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c50bg->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c50bg->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixCorr%dNorm_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c50bg->SaveAs(name);
  }

  c60=new TCanvas("c60","PhiEtaMixEff");
  // hPhiEtaMixEffPlot->SetMaximum(maxeff);
  //hPhiEtaMixEffPlot->SetMinimum(mineff);
  SetMargins2DSurf(c60);
  hPhiEtaMixEffPlot->Draw("lego2");
  sprintf(name,"%s/%s%sPhiEtaMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c60->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c60->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixEff%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c60->SaveAs(name);

   c100=new TCanvas("c100","PhiEtaRaw");
  SetMargins2D(c100);
  hPhiEtaRawPlot->Draw("colz");
  sprintf(name,"%s/%s%sPhiEtaRaw_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c100->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaRaw_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c100->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaRaw_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c100->SaveAs(name);

  c200=new TCanvas("c200","PhiEtaCorr");
  SetMargins2D(c200);
  hPhiEtaCorrPlot->Draw("colz");
  sprintf(name,"%s/%s%sPhiEtaCorr%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c200->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaCorr%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c200->SaveAs(name);
 sprintf(name,"%s/%s%sPhiEtaCorr%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c200->SaveAs(name);
  
   if(!notDelta){
  c200bg=new TCanvas("c200bg","PhiEtaCorrBgSub");
  SetMargins2D(c200bg);
  hPhiEtaCorr2Plot->Draw("colz");
  sprintf(name,"%s/%s%sPhiEtaCorr%dBgSub_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c200bg->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaCorr%dBgSub_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c200bg->SaveAs(name);
 sprintf(name,"%s/%s%sPhiEtaCorr%dBgSub_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c200bg->SaveAs(name);
  }

  c300=new TCanvas("c300","PhiEtaEff");
  SetMargins2D(c300);
  hPhiEtaEffPlot->Draw("colz");
  sprintf(name,"%s/%s%sPhiEtaEff%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c300->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaEff%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c300->SaveAs(name);

  c400=new TCanvas("c400","PhiEtaMixRaw");
  SetMargins2D(c400);
  hPhiEtaMixRawPlot->Draw("colz");
  sprintf(name,"%s/%s%sPhiEtaMixRaw_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c400->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixRaw_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c400->SaveAs(name);
 sprintf(name,"%s/%s%sPhiEtaMixRaw_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c400->SaveAs(name);

  c500=new TCanvas("c500","PhiEtaMixCorr");
  SetMargins2D(c500);
  hPhiEtaMixCorrPlot->Draw("colz");
  sprintf(name,"%s/%s%sPhiEtaMixCorr%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c500->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixCorr%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c500->SaveAs(name);
 sprintf(name,"%s/%s%sPhiEtaMixCorr%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c500->SaveAs(name);
  
   if(!notDelta){
  c500bg=new TCanvas("c500bg","PhiEtaMixCorrNorm");
  SetMargins2D(c500bg);
  hPhiEtaMixCorr2Plot->Draw("colz");
  sprintf(name,"%s/%s%sPhiEtaMixCorr%dNorm_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c500bg->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixCorr%dNorm_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c500bg->SaveAs(name);
 sprintf(name,"%s/%s%sPhiEtaMixCorr%dNorm_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c500bg->SaveAs(name);
  }
   
  c600=new TCanvas("c600","PhiEtaMixEff");
  SetMargins2D(c600);
  hPhiEtaMixEffPlot->Draw("colz");
  sprintf(name,"%s/%s%sPhiEtaMixEff%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c600->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaMixEff%d_2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c600->SaveAs(name);
  

  
  c20Cyl=new TCanvas("c20Cyl","PhiEtaCorrCyl");
  SetMargins2DSurf(c20Cyl);
  c20Cyl->SetTheta(45);
  c20Cyl->SetPhi(200);
  hPhiEtaCorrPlot->Draw("surf1cyl");
  sprintf(name,"%s/%s%sPhiEtaCorr%dCyl_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.gif",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  c20Cyl->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaCorr%dCyl_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.eps",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePS)c20Cyl->SaveAs(name);
  sprintf(name,"%s/%s%sPhiEtaCorr%dCyl_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d%s_z%1.1f.pdf",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,titArray[noTitle],ZYACent);
  if(SavePDF)c20Cyl->SaveAs(name);

 if(SaveRoot){
   sprintf(name,"%s/DrawDPhi%s%sCorr%d_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_R%d_%d%s%s_z%1.1f.root",Folder,cDelta[notDelta],cPt[ptWeighted],EffMethod,TPt1,TPt2,APt1,APt2,Cent1,Cent2,ReBinPhi,ReBinEta,csign[LSign],titArray[noTitle],ZYACent);
 	TFile *fout=new TFile(name,"recreate");
 	cout << "Writing Histograms to File: " << name << endl;
	hNEvents->Write();
	hNTrig->Write();
 	hPhiRawPlot->Write();
 	hPhiEffPlot->Write();
 	hPhiCorrPlot->Write();
 	hPhiMixRawPlot->Write();
 	hPhiMixEffPlot->Write();
 	hPhiMixCorrPlot->Write();
	hEtaARawPlot->Write();
 	hEtaAEffPlot->Write();
 	hEtaACorrPlot->Write();
 	hEtaAMixRawPlot->Write();
 	hEtaAMixEffPlot->Write();
 	hEtaAMixCorrPlot->Write();
	hEtaNRawPlot->Write();
 	hEtaNEffPlot->Write();
 	hEtaNCorrPlot->Write();
 	hEtaNMixRawPlot->Write();
 	hEtaNMixEffPlot->Write();
 	hEtaNMixCorrPlot->Write();
 	hPhiEtaRawPlot->Write();
 	hPhiEtaEffPlot->Write();
 	hPhiEtaCorrPlot->Write();
 	hPhiEtaMixRawPlot->Write();
 	hPhiEtaMixEffPlot->Write();
 	hPhiEtaMixCorrPlot->Write();
 	if(!notDelta){
 		hPhiCorr2Plot->Write();
 		hPhiMixCorr2Plot->Write();
		hEtaNCorr2Plot->Write();
 		hEtaNMixCorr2Plot->Write();
		hEtaACorr2Plot->Write();
 		hEtaAMixCorr2Plot->Write();
 		hPhiEtaCorr2Plot->Write();
 		hPhiEtaMixCorr2Plot->Write();
 	}
 	fout->Write();
 	fout->Close();
}
}
