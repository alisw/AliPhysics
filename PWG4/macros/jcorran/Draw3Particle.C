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
//Draw 3-Particle correlations plots (using output of AliAnalysisTaskDiHadron)
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
void Draw3Particle(float TPt1, float TPt2, float APt1, float APt2, float Cent1, float Cent2, int LSign=0){//ptWeighted 0 no weight, 1 weighted, 2 <pt>; Lsign 0-all, 1-like-sign 2-unlike-sign
  gROOT->Reset();
  gROOT->LoadMacro("Util9.C");
  Style(1);
  int SaveRoot=1;
  int SavePS=1;
  int SavePDF=0;
  float ZYACent=1.5;//-0to-1.99 2-particle stripes -2 two-particle 3-D min, -3 3-particle PhiPhi Min
  float ZYAWidth=0.2; 
  int ReBinPhi=1;
  int ReBinEta=1;
  int DrawAsHisto=0;
  int noTitle=1;
  int SaveMacro=0;
  int SaveText=1;
  int DrawMC=1;
  int EffMethod=0;//0 no correction, 1 angle dependent, 2 use mix for triggered, 3 <1>, 4 fits, 5  fits with VariablePtLimit
  //int SaveEffFits=1;
  char *fileType=".gif";
  
  //char *Folder="2010-08-17/LHC10b_7pass2";
  //char *Folder="2010-08-17/LHC10c_900pass2";
  //char *Folder="2010-08-17/LHC10c6_900Pythia";
  char *Folder="2010-08-17/7Pythia_LHC10b5"; 

  //char *EffFolder="2010-08-17/LHC10c6_900Pythia";
  char *EffFolder="2010-08-17/7Pythia_LHC10b5"; 
  
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
  
  TH2F *hPhiPhiRaw=new TH2F("hPhiPhiRaw","",1,0,1,1,0,1);
  TH2F *hPhiPhiCorr=new TH2F("hPhiPhiCorr","",1,0,1,1,0,1);
  TH2F *hPhiPhiEff=new TH2F("hPhiPhiEff","",1,0,1,1,0,1);
  TH2F *hPhiPhiMC=new TH2F("hPhiPhiMC","",1,0,1,1,0,1);
  TH2F *hPhiPhiSSRaw=new TH2F("hPhiPhiSSRaw","",1,0,1,1,0,1);
  TH2F *hPhiPhiSSCorr=new TH2F("hPhiPhiSSCorr","",1,0,1,1,0,1);
  TH2F *hPhiPhiSSEff=new TH2F("hPhiPhiSSEff","",1,0,1,1,0,1);
  TH2F *hPhiPhiSSMC=new TH2F("hPhiPhiSSMC","",1,0,1,1,0,1);
  TH2F *hPhiPhiMixRaw=new TH2F("hPhiPhiMixRaw","",1,0,1,1,0,1);
  TH2F *hPhiPhiMixCorr=new TH2F("hPhiPhiMixCorr","",1,0,1,1,0,1);
  TH2F *hPhiPhiMixEff=new TH2F("hPhiPhiMixEff","",1,0,1,1,0,1);
  TH2F *hPhiPhiMixMC=new TH2F("hPhiPhiMixMC","",1,0,1,1,0,1);
  
  TH2F *hEtaEtaRaw=new TH2F("hEtaEtaRaw","",1,0,1,1,0,1);
  TH2F *hEtaEtaCorr=new TH2F("hEtaEtaCorr","",1,0,1,1,0,1);
  TH2F *hEtaEtaEff=new TH2F("hEtaEtaEff","",1,0,1,1,0,1);
  TH2F *hEtaEtaMC=new TH2F("hEtaEtaMC","",1,0,1,1,0,1);
  TH2F *hEtaEtaSSRaw=new TH2F("hEtaEtaSSRaw","",1,0,1,1,0,1);
  TH2F *hEtaEtaSSCorr=new TH2F("hEtaEtaSSCorr","",1,0,1,1,0,1);
  TH2F *hEtaEtaSSEff=new TH2F("hEtaEtaSSEff","",1,0,1,1,0,1);
  TH2F *hEtaEtaSSMC=new TH2F("hEtaEtaSSMC","",1,0,1,1,0,1);
  TH2F *hEtaEtaMixRaw=new TH2F("hEtaEtaMixRaw","",1,0,1,1,0,1);
  TH2F *hEtaEtaMixCorr=new TH2F("hEtaEtaMixCorr","",1,0,1,1,0,1);
  TH2F *hEtaEtaMixEff=new TH2F("hEtaEtaMixEff","",1,0,1,1,0,1);
  TH2F *hEtaEtaMixMC=new TH2F("hEtaEtaMixMC","",1,0,1,1,0,1);

  TH1F *hPhiPhiNOn=new TH1F("hPhiPhiNOn","",1,0,1);
  TH1F *hPhiPhiAOn=new TH1F("hPhiPhiAOn","",1,0,1);
  TH1F *hPhiPhiNOff=new TH1F("hPhiPhiNOff","",1,0,1);
  TH1F *hPhiPhiAOff=new TH1F("hPhiPhiAOff","",1,0,1);

  TH1F *hEtaEtaNOn=new TH1F("hEtaEtaNOn","",1,0,1);
    TH1F *hEtaEtaNOff=new TH1F("hEtaEtaNOff","",1,0,1);

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
    MakeProjections(TPt1,TPt2,APt1,APt2,Cent,inList,hPhiRaw,hPhiMixRaw,hEtaNRaw,hEtaNMixRaw,hEtaARaw,hEtaAMixRaw,hPhiEtaRaw,hPhiEtaMixRaw,MPt,0,0,0,LSign);
    MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiEff,hPhiMixEff,hEtaNEff,hEtaNMixEff,hEtaAEff,hEtaAMixEff,hPhiEtaEff,hPhiEtaMixEff,MPt2,0,0,0,LSign);
    MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiMC,hPhiMixMC,hEtaNMC,hEtaNMixMC,hEtaAMC,hEtaAMixMC,hPhiEtaMC,hPhiEtaMixMC,MPt2,0,1,0,LSign);
    
    if(EffMethod<4){
      EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hPhiEff,hPhiMC,hPhiMixEff,hPhiMixMC,hEtaNEff,hEtaNMC,hEtaNMixEff,hEtaNMixMC,hEtaAEff,hEtaAMC,hEtaAMixEff,hEtaAMixMC,hPhiEtaEff,hPhiEtaMC,hPhiEtaMixEff,hPhiEtaMixMC,EffMethod);
    }
    else{
      EffFit(APt1,APt2,Cent,effList,hPhiEff,hPhiMixEff,hEtaNEff,hEtaNMixEff,hEtaAEff,hEtaAMixEff,hPhiEtaEff,hPhiEtaMixEff,LSign,VariablePtLimit);
    }
    //inFile->Close();
    // TFile *inFile2=new TFile(inName);
    // TList *inList2=inFile2->Get("julery_DiHadron");
    Load3Particle(TPt1,TPt2,APt1,APt2,Cent,inList,hPhiPhiRaw,hPhiPhiSSRaw,hPhiPhiMixRaw,hEtaEtaRaw,hEtaEtaSSRaw,hEtaEtaMixRaw,0,LSign);
    Load3Particle(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiPhiEff,hPhiPhiSSEff,hPhiPhiMixEff,hEtaEtaEff,hEtaEtaSSEff,hEtaEtaMixEff,0,LSign);
    Load3Particle(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiPhiMC,hPhiPhiSSMC,hPhiPhiMixMC,hEtaEtaMC,hEtaEtaSSMC,hEtaEtaMixMC,1,LSign);
    
    EffCorr3Part(TPt1,TPt2,APt1,APt2,Cent,hPhiPhiEff,hPhiPhiMC,hPhiPhiSSEff,hPhiPhiSSMC,hPhiPhiMixEff,hPhiPhiMixMC,hEtaEtaEff,hEtaEtaMC,hEtaEtaSSEff,hEtaEtaSSMC,hEtaEtaMixEff,hEtaEtaMixMC,0);
    
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
    
    hPhiPhiCorr=(TH2F*)hPhiPhiRaw->Clone();
    sprintf(name,"hPhiPhiCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hPhiPhiCorr->SetName(name);
    hPhiPhiCorr->Divide(hPhiPhiEff);
    
    hPhiPhiSSCorr=(TH2F*)hPhiPhiSSRaw->Clone();
    sprintf(name,"hPhiPhiSSCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hPhiPhiSSCorr->SetName(name);
    hPhiPhiSSCorr->Divide(hPhiPhiSSEff);
    
    hPhiPhiMixCorr=(TH2F*)hPhiPhiMixRaw->Clone();
    sprintf(name,"hPhiPhiMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hPhiPhiMixCorr->SetName(name);
    cout << "hPhiPhiMixCorr: " << hPhiPhiMixCorr->GetBinContent(1,1) << " " << hPhiPhiMixCorr->GetBinError(1,1) << endl;
    hPhiPhiMixCorr->Divide(hPhiPhiMixEff);
    hPhiPhiMixCorr->Scale(hPhiPhiSSCorr->GetSum()/hPhiPhiMixCorr->GetSum());
     cout << "hPhiPhiMixCorr: " << hPhiPhiMixCorr->GetBinContent(1,1) << " " << hPhiPhiMixCorr->GetBinError(1,1) << endl;
     cout << "hPhiPhiMixEff: " << hPhiPhiMixEff->GetBinContent(1,1) << " " << hPhiPhiMixEff->GetBinError(1,1) << endl;

    hEtaEtaCorr=(TH2F*)hEtaEtaRaw->Clone();
    sprintf(name,"hEtaEtaCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hEtaEtaCorr->SetName(name);
    hEtaEtaCorr->Divide(hEtaEtaEff);
    
    hEtaEtaSSCorr=(TH2F*)hEtaEtaSSRaw->Clone();
    sprintf(name,"hEtaEtaSSCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hEtaEtaSSCorr->SetName(name);
    hEtaEtaSSCorr->Divide(hEtaEtaSSEff);
    
    hEtaEtaMixCorr=(TH2F*)hEtaEtaMixRaw->Clone();
    sprintf(name,"hEtaEtaMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hEtaEtaMixCorr->SetName(name);
    
    hEtaEtaMixCorr->Divide(hEtaEtaMixEff);
    hEtaEtaMixCorr->Scale(hEtaEtaSSCorr->GetSum()/hEtaEtaMixCorr->GetSum());
  
    
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
    for(int x=1;x<=hPhiPhiCorr->GetNbinsX();x++){
      for(int y=1;y<=hPhiPhiCorr->GetNbinsY();y++){
	if(hPhiPhiCorr->GetBinContent(x,y)<1E-10)hPhiPhiCorr->SetBinContent(x,y,0);
	if(hPhiPhiSSCorr->GetBinContent(x,y)<1E-10)hPhiPhiSSCorr->SetBinContent(x,y,0);
	if(hPhiPhiMixCorr->GetBinContent(x,y)<1E-10)hPhiPhiMixCorr->SetBinContent(x,y,0);
      }
    }
    for(int x=1;x<=hEtaEtaCorr->GetNbinsX();x++){
      for(int y=1;y<=hEtaEtaCorr->GetNbinsY();y++){
	if(hEtaEtaCorr->GetBinContent(x,y)<1E-10)hEtaEtaCorr->SetBinContent(x,y,0);
	if(hEtaEtaSSCorr->GetBinContent(x,y)<1E-10)hEtaEtaSSCorr->SetBinContent(x,y,0);
	if(hEtaEtaMixCorr->GetBinContent(x,y)<1E-10)hEtaEtaMixCorr->SetBinContent(x,y,0);
      }
    }

    // GeoCorr2(hEtaNCorr,hEtaNMixCorr,hEtaACorr,hEtaAMixCorr,hPhiEtaCorr,hPhiEtaMixCorr);
    //GeoCorr3Part2(hEtaEtaCorr,hEtaEtaSSCorr,hEtaEtaMixCorr, hEtaEtaHSCorr);
    MixedCorrect(hPhiCorr,hPhiMixCorr,hEtaNCorr,hEtaNMixCorr,hEtaACorr,hEtaAMixCorr,hPhiEtaCorr,hPhiEtaMixCorr);
    MixedCorrect3(hPhiPhiCorr,hPhiPhiSSCorr,hPhiPhiMixCorr,hEtaEtaCorr,hEtaEtaSSCorr,hEtaEtaMixCorr);

    float a[3];
    Float_t b=(hPhiPhiCorr->GetSum()/pow(hPhiCorr->GetSum(),2))/(hPhiPhiSSCorr->GetSum()/pow(hPhiMixCorr->GetSum(),2));
    //  b=1;
    cout << "b:= " << b << endl;
 float bEta=(hEtaEtaCorr->GetSum()/pow(hEtaNCorr->GetSum(),2))/(hEtaEtaSSCorr->GetSum()/pow(hEtaNMixCorr->GetSum(),2));
 float bEta2=bEta;
 //bEta=b;
    cout << "bEta:= " << bEta << endl;
    if(ZYACent>0){
      ZYA1(hPhiCorr,hPhiMixCorr,a,ZYACent,ZYAWidth);
    cout << "a: "  << a[0] << "  Error: " << a[1] << "  Error Code: " << a[2] << "  <pT>: " << MPt[0]  << endl;
    }
    else if(ZYACent>-2){
      ZYAM2D2(hPhiEtaCorr,hPhiEtaMixCorr,a,fabs(ZYACent),ZYAWidth);
    cout << "a: "  << a[0] << "  Error: " << a[1] << "  Error Code: " << a[2] << "  <pT>: " << MPt[0]  << endl;
    }
    else if(fabs(ZYACent+2)<0.1) ZYAM2D(hPhiEtaCorr,hPhiEtaMixCorr,a,72,10);
    else if(fabs(ZYACent+3)<0.1) ZYAM3(hPhiCorr, hPhiMixCorr, hPhiPhiCorr,hPhiPhiSSCorr,hPhiPhiMixCorr, a, b,36,50);

    hPhiMixCorr->Scale(a[0]);
    hEtaNMixCorr->Scale(a[0]);
    float scale3=a[0]*a[0]*b;
    //bEta2=0;
    float scale3Eta=a[0]*a[0]*bEta;
    float scale3Eta2=a[0]*a[0]*bEta2;
    hPhiPhiMixCorr->Scale(scale3);
    hPhiPhiSSCorr->Scale(scale3);
    // hEtaEtaSSCorr->Add(hEtaEtaMixCorr,-1);
    hEtaEtaMixCorr->Scale(scale3Eta2);
    hEtaEtaSSCorr->Scale(scale3Eta);
    // hEtaEtaSSCorr->Add(hEtaEtaMixCorr);

    TH2F *hPhiPhiHSCorr=(TH2F*)hPhiPhiCorr->Clone();
    hPhiPhiHSCorr->SetName("hPhiPhiHSCorr");
    hPhiPhiHSCorr->SetTitle("Hard-Soft Term");
    
    /*
    TH2F *hEtaEtaHSCorr=(TH2F*)hEtaEtaMixRaw->Clone();
    hEtaEtaHSCorr->SetName("hEtaEtaHSCorr");
    hEtaEtaHSCorr->SetTitle("Hard-Soft Term");
    */
    float binsx=hEtaEtaCorr->GetNbinsX();
    float minx=hEtaEtaCorr->GetBinCenter(1)-hEtaEtaCorr->GetBinWidth(1)/2;
    float maxx=hEtaEtaCorr->GetBinCenter(binsx)+hEtaEtaCorr->GetBinWidth(binsx)/2;
    TH2F *hEtaEtaHSCorr=new TH2F("hEtaEtaHCCorr","Hard-Soft Term",binsx,minx,maxx,binsx,minx,maxx);
    hEtaEtaHSCorr->Sumw2();
 

    TH1F *hPhiCorr2=(TH1F*)hPhiCorr->Clone();
    hPhiCorr2->SetName("hPhiCorr2");
    hPhiCorr2->Add(hPhiMixCorr,-1);
    
    TH1F *hEtaNCorr2=(TH1F*)hEtaNCorr->Clone();
    hEtaNCorr2->SetName("hEtaNCorr2");
    hEtaNCorr2->Add(hEtaNMixCorr,-1);
    
    float sigX, sigY, esigX, esigY, sig;
    float mixX, mixY, emixX, emixY, esig;
    int binRatio=hPhiCorr->GetNbinsX()/hPhiPhiCorr->GetNbinsX();
    float CheckRatio=hPhiCorr->GetNbinsX()/hPhiPhiCorr->GetNbinsX();
    if(fabs(binRatio-CheckRatio)>0.01){
      cout << "Warning Use Bins for 3-Particle That are Divisors of 2-Particle!!!!!!    Phi" << endl;
}
    // cout << "binRatio " << binRatio << endl;
    for(int x=1;x<=hPhiPhiCorr->GetNbinsX();x++){
      sigX=0; esigX=0; mixX=0; emixX=0;
      for(int j=(binRatio*(x-1)+1);j<=(x*binRatio);j++){
      sigX+=hPhiCorr2->GetBinContent(j);
      esigX+=pow(hPhiCorr2->GetBinError(j),2);
      mixX+=hPhiMixCorr->GetBinContent(j);
      emixX+=pow(hPhiMixCorr->GetBinError(j),2);
      }
      for(int y=1;y<=hPhiPhiCorr->GetNbinsX();y++){
	sigY=0; esigY=0; mixY=0; emixY=0;
	for(int k=(binRatio*(y-1)+1);k<=(y*binRatio);k++){
	sigY+=hPhiCorr2->GetBinContent(k);
	esigY+=pow(hPhiCorr2->GetBinError(k),2);
	mixY+=hPhiMixCorr->GetBinContent(k);
	emixY+=pow(hPhiMixCorr->GetBinError(k),2);
	}
	//cout << sigX << " " << sigY << " " << esigX << " " << esigY << " " << mixX << " " << emixX << " " << mixY << " " << emixY << endl;
	sig=sigX*mixY+sigY*mixX;
	esig=pow(pow(sigX*mixY,2)*(esigX/sigX/sigX+emixY/mixY/mixY)+pow(sigY*mixX,2)*(esigY/sigY/sigY+emixX/mixX/mixX),0.5);
	sig=sig/pow(binRatio,2);
	esig=esig/pow(binRatio,2);
	hPhiPhiHSCorr->SetBinContent(x,y,sig);
	hPhiPhiHSCorr->SetBinError(x,y,esig);
      }
    }

 float sigX, sigY, esigX, esigY, sig;
    float mixX, mixY, emixX, emixY, esig;
    float centX, centY;
    int binRatio=hEtaNCorr->GetNbinsX()/hEtaEtaCorr->GetNbinsX();
    float CheckRatio=hEtaNCorr->GetNbinsX()/hEtaEtaCorr->GetNbinsX();
    float etaCut=hEtaNCorr->GetBinCenter(hEtaNCorr->GetNbinsX())+hEtaNCorr->GetBinWidth(1);
    if(fabs(binRatio-CheckRatio)>0.01){
   cout << "Warning Use Bins for 3-Particle That are Divisors of 2-Particle!!!!!!   Eta" << endl;
    }
    for(int x=1;x<=hEtaNCorr->GetNbinsX();x++){
      sigX=0; esigX=0; mixX=0; emixX=0;
      centX=hEtaEtaCorr->GetXaxis()->GetBinCenter(x);
      for(int j=(binRatio*(x-1)+1);j<=(x*binRatio);j++){
	sigX+=hEtaNCorr2->GetBinContent(j);
	esigX+=hEtaNCorr2->GetBinError(j);
	mixX+=hEtaNMixCorr->GetBinContent(j);
	emixX+=hEtaNMixCorr->GetBinError(j);
      }
      for(int y=1;y<=hEtaNCorr->GetNbinsX();y++){
	sigY=0; esigY=0; mixY=0; emixY=0;
	centY=hEtaEtaCorr->GetXaxis()->GetBinCenter(y);
	for(int k=(binRatio*(y-1)+1);k<=(y*binRatio);k++){
	  sigY+=hEtaNCorr2->GetBinContent(k);
	  esigY+=hEtaNCorr2->GetBinError(k);
	  mixY+=hEtaNMixCorr->GetBinContent(k);
	  emixY+=hEtaNMixCorr->GetBinError(k);
	}
	sig=sigX*mixY+sigY*mixX;

	//divide by 0 protection
	if(sigX==0)sigX=1E-5;
	if(sigY==0)sigY=1E-5;
	if(mixX==0)mixX=1E-5;
	if(mixY==0)mixY=1E-5;
	esig=pow(pow(sigX*mixY,2)*(pow(esigX/sigX,2)+pow(emixY/mixY,2))+pow(sigY*mixX,2)*(pow(esigY/sigY,2)+pow(emixX/mixX,2)),0.5);
	//if(fabs(centX-centY)>etaCut){sig=0;esig=0;}
	if(hEtaEtaMixCorr->GetBinContent(x,y)){
	hEtaEtaHSCorr->SetBinContent(x,y,sig/binRatio/binRatio);
	hEtaEtaHSCorr->SetBinError(x,y,esig/binRatio/binRatio);
	}
	//else cout << hEtaEtaSSCorr->GetBinContent(x,y)  << " " << hEtaEtaHSCorr->GetBinError(x,y) << endl;
      }
    }
    
   
  
    // hEtaNCorr2=(TH1F*)hEtaNCorr->Clone();
    //  hEtaNCorr2->SetName("hEtaNCorr2");
    //  hEtaNCorr2->Add(hEtaNMixCorr,-1);    


    TH2F *hPhiPhiCorr2=(TH2F*)hPhiPhiCorr->Clone();
    hPhiPhiCorr2->SetName("hPhiPhiCorr2");
    hPhiPhiCorr2->SetTitle("Background Subtracted");
    hPhiPhiCorr2->Add(hPhiPhiHSCorr,-1);
    hPhiPhiCorr2->Add(hPhiPhiSSCorr,-1);

    TH2F *hEtaEtaCorr2=(TH2F*)hEtaEtaCorr->Clone();
    hEtaEtaCorr2->SetName("hEtaEtaCorr2");
    hEtaEtaCorr2->SetTitle("Background Subtracted");
    hEtaEtaCorr2->Add(hEtaEtaHSCorr,-1);
    hEtaEtaCorr2->Add(hEtaEtaSSCorr,-1);

    //here need to sum over the centralities

    ProjPhiPhi(hPhiPhiCorr2,hPhiPhiNOn,hPhiPhiNOff,hPhiPhiAOn,hPhiPhiAOff,1,0.35);//0.35
    ProjEtaEta(hEtaEtaCorr2,hEtaEtaNOn,hEtaEtaNOff,0.35);

    cPhiAll=new TCanvas("cPhiAll","3-Particle Phi",1500,900);
    cPhiAll->Divide(4,3);
    //  SetMargins1D(cPhiAll_1);
    //1-D plots
    float RM1=0.02;
    float LM1=0.16;
    float TM1=0.1;
    float BM1=0.15;
    //colz
    float RM2=0.25;
    float LM2=0.12;
    float TM2=0.1;
    float BM2=0.15;
    //surf
    float RM3=0.05;
    float LM3=0.21;
    float TM3=0.05;
    float BM3=0.15;

    cPhiAll_1->SetRightMargin(RM1);
    cPhiAll_1->SetLeftMargin(LM1);
    cPhiAll_1->SetTopMargin(TM1);
    cPhiAll_1->SetBottomMargin(BM1);

    cPhiAll_2->SetRightMargin(RM1);
    cPhiAll_2->SetLeftMargin(LM1);
    cPhiAll_2->SetTopMargin(TM1);
    cPhiAll_2->SetBottomMargin(BM1);

    cPhiAll_3->SetRightMargin(RM3);
    cPhiAll_3->SetLeftMargin(LM3);
    cPhiAll_3->SetTopMargin(TM3);
    cPhiAll_3->SetBottomMargin(BM3);

    cPhiAll_4->SetRightMargin(RM2);
    cPhiAll_4->SetLeftMargin(LM2);
    cPhiAll_4->SetTopMargin(TM2);
    cPhiAll_4->SetBottomMargin(BM2);
    
    cPhiAll_5->SetRightMargin(RM3);
    cPhiAll_5->SetLeftMargin(LM3);
    cPhiAll_5->SetTopMargin(TM3);
    cPhiAll_5->SetBottomMargin(BM3);

    cPhiAll_6->SetRightMargin(RM2);
    cPhiAll_6->SetLeftMargin(LM2);
    cPhiAll_6->SetTopMargin(TM2);
    cPhiAll_6->SetBottomMargin(BM2);

    cPhiAll_7->SetRightMargin(RM3);
    cPhiAll_7->SetLeftMargin(LM3);
    cPhiAll_7->SetTopMargin(TM3);
    cPhiAll_7->SetBottomMargin(BM3);

    cPhiAll_8->SetRightMargin(RM2);
    cPhiAll_8->SetLeftMargin(LM2);
    cPhiAll_8->SetTopMargin(TM2);
    cPhiAll_8->SetBottomMargin(BM2);

    cPhiAll_9->SetRightMargin(RM3);
    cPhiAll_9->SetLeftMargin(LM3);
    cPhiAll_9->SetTopMargin(TM3);
    cPhiAll_9->SetBottomMargin(BM3);

    cPhiAll_10->SetRightMargin(RM2);
    cPhiAll_10->SetLeftMargin(LM2);
    cPhiAll_10->SetTopMargin(TM2);
    cPhiAll_10->SetBottomMargin(BM2);
    
    cPhiAll_11->SetRightMargin(RM1);
    cPhiAll_11->SetLeftMargin(LM1);
    cPhiAll_11->SetTopMargin(TM1);
    cPhiAll_11->SetBottomMargin(BM1);

    cPhiAll_12->SetRightMargin(RM1);
    cPhiAll_12->SetLeftMargin(LM1);
    cPhiAll_12->SetTopMargin(TM1);
    cPhiAll_12->SetBottomMargin(BM1);

    hPhiPhiCorr->SetTitle("Raw Signal");
    hPhiPhiSSCorr->SetTitle("Soft-Soft Term");
    hPhiPhiNOn->SetTitle("Near");
    hPhiPhiNOff->SetTitle("Near");
    hPhiPhiAOn->SetTitle("Away");
    hPhiPhiAOff->SetTitle("Away");

    cPhiAll_1->cd();
    SetTitles1D(hPhiCorr);
    hPhiCorr->Draw();
    hPhiMixCorr->Draw("same");
    
    cPhiAll_2->cd();
    SetTitles1D(hPhiCorr2);
    hPhiCorr2->Draw();
    ZeroLine->Draw("same");
    
    cPhiAll_3->cd();
    SetTitles2D(hPhiPhiCorr);
    hPhiPhiCorr->Draw("lego2");

    cPhiAll_4->cd();
    hPhiPhiCorr->Draw("colz");

    cPhiAll_5->cd();
    SetTitles2D(hPhiPhiHSCorr);
    hPhiPhiHSCorr->Draw("lego2");
    
    cPhiAll_6->cd();
    hPhiPhiHSCorr->Draw("colz");

    cPhiAll_7->cd();
    SetTitles2D(hPhiPhiMixCorr);
    hPhiPhiSSCorr->Draw("lego2");
    
    cPhiAll_8->cd();
    hPhiPhiSSCorr->Draw("colz");

    cPhiAll_9->cd();
    SetTitles2D(hPhiPhiCorr2);
    hPhiPhiCorr2->Draw("lego2");
    
    cPhiAll_10->cd();
    hPhiPhiCorr2->Draw("colz");

    cPhiAll_11->cd();
    hPhiPhiNOff->Draw();
    hPhiPhiNOn->Draw("same");
    keySymbol(0.2,0.85,"On",4,25,0.05,1.5);
    keySymbol(0.2,0.8,"Off",2,20,0.05,1.5);

    cPhiAll_12->cd();
    hPhiPhiAOff->Draw();
    hPhiPhiAOn->Draw("same");

    sprintf(name,"%s/Draw3Particle_PhiAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhiAll->SaveAs(name);

    cEtaAll=new TCanvas("cEtaAll","3-Particle Eta",1500,900);
    cEtaAll->Divide(4,3);
    //  SetMargins1D(cEtaAll_1);
    /*
    //1-D plots
    float RM1=0.02;
    float LM1=0.16;
    float TM1=0.1;
    float BM1=0.15;
    //colz
    float RM2=0.24;
    float LM2=0.12;
    float TM2=0.1;
    float BM2=0.15;
    //surf
    float RM3=0.05;
    float LM3=0.2;
    float TM3=0.05;
    float BM3=0.15;
    */

    cEtaAll_1->SetRightMargin(RM1);
    cEtaAll_1->SetLeftMargin(LM1);
    cEtaAll_1->SetTopMargin(TM1);
    cEtaAll_1->SetBottomMargin(BM1);

    cEtaAll_2->SetRightMargin(RM1);
    cEtaAll_2->SetLeftMargin(LM1);
    cEtaAll_2->SetTopMargin(TM1);
    cEtaAll_2->SetBottomMargin(BM1);

    cEtaAll_3->SetRightMargin(RM3);
    cEtaAll_3->SetLeftMargin(LM3);
    cEtaAll_3->SetTopMargin(TM3);
    cEtaAll_3->SetBottomMargin(BM3);

    cEtaAll_4->SetRightMargin(RM2);
    cEtaAll_4->SetLeftMargin(LM2);
    cEtaAll_4->SetTopMargin(TM2);
    cEtaAll_4->SetBottomMargin(BM2);
    
    cEtaAll_5->SetRightMargin(RM3);
    cEtaAll_5->SetLeftMargin(LM3);
    cEtaAll_5->SetTopMargin(TM3);
    cEtaAll_5->SetBottomMargin(BM3);

    cEtaAll_6->SetRightMargin(RM2);
    cEtaAll_6->SetLeftMargin(LM2);
    cEtaAll_6->SetTopMargin(TM2);
    cEtaAll_6->SetBottomMargin(BM2);

    cEtaAll_7->SetRightMargin(RM3);
    cEtaAll_7->SetLeftMargin(LM3);
    cEtaAll_7->SetTopMargin(TM3);
    cEtaAll_7->SetBottomMargin(BM3);

    cEtaAll_8->SetRightMargin(RM2);
    cEtaAll_8->SetLeftMargin(LM2);
    cEtaAll_8->SetTopMargin(TM2);
    cEtaAll_8->SetBottomMargin(BM2);

    cEtaAll_9->SetRightMargin(RM3);
    cEtaAll_9->SetLeftMargin(LM3);
    cEtaAll_9->SetTopMargin(TM3);
    cEtaAll_9->SetBottomMargin(BM3);

    cEtaAll_10->SetRightMargin(RM2);
    cEtaAll_10->SetLeftMargin(LM2);
    cEtaAll_10->SetTopMargin(TM2);
    cEtaAll_10->SetBottomMargin(BM2);
    
    cEtaAll_11->SetRightMargin(RM1);
    cEtaAll_11->SetLeftMargin(LM1);
    cEtaAll_11->SetTopMargin(TM1);
    cEtaAll_11->SetBottomMargin(BM1);

    cEtaAll_12->SetRightMargin(RM1);
    cEtaAll_12->SetLeftMargin(LM1);
    cEtaAll_12->SetTopMargin(TM1);
    cEtaAll_12->SetBottomMargin(BM1);

    hEtaEtaCorr->SetTitle("Raw Signal");
    hEtaEtaSSCorr->SetTitle("Soft-Soft Term");
    hEtaEtaNOn->SetTitle("Projections");
    hEtaEtaNOff->SetTitle("Projections");
 

    cEtaAll_1->cd();
    SetTitles1D(hEtaNCorr);
    hEtaNCorr->Draw();
    hEtaNMixCorr->Draw("same");
    
    cEtaAll_2->cd();
    SetTitles1D(hEtaNCorr2);
    hEtaNCorr2->Draw();
    ZeroLine->Draw("same");
    
    cEtaAll_3->cd();
    SetTitles2D(hEtaEtaCorr);
    hEtaEtaCorr->Draw("lego2");

    cEtaAll_4->cd();
    hEtaEtaCorr->Draw("colz");

    cEtaAll_5->cd();
    SetTitles2D(hEtaEtaHSCorr);
    hEtaEtaHSCorr->Draw("lego2");
    
    cEtaAll_6->cd();
    hEtaEtaHSCorr->Draw("colz");

    cEtaAll_7->cd();
    SetTitles2D(hEtaEtaMixCorr);
    hEtaEtaSSCorr->Draw("lego2");
    
    cEtaAll_8->cd();
    hEtaEtaSSCorr->Draw("colz");

    cEtaAll_9->cd();
    SetTitles2D(hEtaEtaCorr2);
    hEtaEtaCorr2->Draw("lego2");
    
    cEtaAll_10->cd();
    hEtaEtaCorr2->Draw("colz");

    cEtaAll_11->cd();
    hEtaEtaNOff->Draw();
    hEtaEtaNOn->Draw("same");
    keySymbol(0.2,0.85,"On",4,25,0.05,1.5);
    keySymbol(0.2,0.8,"Off",2,20,0.05,1.5);


    sprintf(name,"%s/Draw3Particle_EtaAll_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEtaAll->SaveAs(name);

    cPhi1=new TCanvas("cPhi1","2-Particle",400,300);
    SetMargins1D(cPhi1);
    hPhiCorr->SetTitle("");
    hPhiCorr->Draw();
    hPhiMixCorr->Draw("same");
    sprintf(name,"%s/Draw3Particle_PhiCorr_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi1->SaveAs(name);
    
    cPhi2=new TCanvas("cPhi2","2-Particle BgSub",400,350);
    SetMargins1D(cPhi2);
    hPhiCorr2->SetTitle("");
    hPhiCorr2->Draw();
    ZeroLine->Draw("same");
    sprintf(name,"%s/Draw3Particle_PhiCorr2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi2->SaveAs(name);

    cPhi3=new TCanvas("cPhi3","3-Particle lego",400,300);
    SetMargins2DSurf(cPhi3);
    hPhiPhiCorr->SetTitle("");
    hPhiPhiCorr->Draw("lego2");
    sprintf(name,"%s/Draw3Particle_PhiPhiCorrLego_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi3->SaveAs(name);

    cPhi4=new TCanvas("cPhi4","3-Particle",400,350);
    SetMargins2D(cPhi4);
    hPhiPhiCorr->SetTitle("");
    hPhiPhiCorr->Draw("colz");
    sprintf(name,"%s/Draw3Particle_PhiPhiCorr_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi4->SaveAs(name);

    cPhi5=new TCanvas("cPhi5","HS lego",400,300);
    SetMargins2DSurf(cPhi5);
    hPhiPhiHSCorr->SetTitle("");
    hPhiPhiHSCorr->Draw("lego2");
    sprintf(name,"%s/Draw3Particle_PhiPhiHSCorrLego_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi5->SaveAs(name);

    cPhi6=new TCanvas("cPhi6","HS",400,350);
    SetMargins2D(cPhi6);
    hPhiPhiHSCorr->SetTitle("");
    hPhiPhiHSCorr->Draw("colz");
    sprintf(name,"%s/Draw3Particle_PhiPhiHSCorr_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi6->SaveAs(name);    

     cPhi7=new TCanvas("cPhi3","SS lego",400,300);
    SetMargins2DSurf(cPhi7);
    hPhiPhiSSCorr->SetTitle("");
    hPhiPhiSSCorr->Draw("lego2");
    sprintf(name,"%s/Draw3Particle_PhiPhiSSCorrLego_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi7->SaveAs(name);

    cPhi8=new TCanvas("cPhi8","SS",400,350);
    SetMargins2D(cPhi8);
    hPhiPhiSSCorr->SetTitle("");
    hPhiPhiSSCorr->Draw("colz");
    sprintf(name,"%s/Draw3Particle_PhiPhiSSCorr_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi8->SaveAs(name);

    cPhi9=new TCanvas("cPhi9","3-Particle BGSub lego",400,300);
    SetMargins2DSurf(cPhi9);
    hPhiPhiCorr2->SetTitle("");
    hPhiPhiCorr2->Draw("lego2");
    sprintf(name,"%s/Draw3Particle_PhiPhiCorr2Lego_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi9->SaveAs(name);

    cPhi10=new TCanvas("cPhi10","3-Particle BGSub",400,350);
    SetMargins2D(cPhi10);
    hPhiPhiCorr2->SetTitle("");
    hPhiPhiCorr2->Draw("colz");
    sprintf(name,"%s/Draw3Particle_PhiPhiCorr2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi10->SaveAs(name);

     cPhi11=new TCanvas("cPhi11","Near Projection",400,300);
     SetMargins1D(cPhi11);
     hPhiPhiNOff->SetTitle("");
     hPhiPhiNOff->Draw();
    hPhiPhiNOn->Draw("same");
    keySymbol(0.2,0.85,"On",4,25,0.05,1.5);
    keySymbol(0.2,0.8,"Off",2,20,0.05,1.5);
    sprintf(name,"%s/Draw3Particle_PhiProjN_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi11->SaveAs(name);
    
    cPhi12=new TCanvas("cPhi12","Away Projection",400,300);
    SetMargins1D(cPhi12);
    hPhiPhiAOff->SetTitle("");
    hPhiPhiAOff->Draw();
    hPhiPhiAOn->Draw("same");
    keySymbol(0.2,0.85,"On",4,25,0.05,1.5);
    keySymbol(0.2,0.8,"Off",2,20,0.05,1.5);
    sprintf(name,"%s/Draw3Particle_PhiProjA_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cPhi12->SaveAs(name);

   cEta1=new TCanvas("cEta1","2-Particle",400,300);
    SetMargins1D(cEta1);
    hEtaNCorr->SetTitle("");
    hEtaNCorr->Draw();
    hEtaNMixCorr->Draw("same");
    sprintf(name,"%s/Draw3Particle_EtaNCorr_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta1->SaveAs(name);
    
    cEta2=new TCanvas("cEta2","2-Particle BgSub",400,350);
    SetMargins1D(cEta2);
    hEtaNCorr2->SetTitle("");
    hEtaNCorr2->Draw();
    ZeroLine->Draw("same");
    sprintf(name,"%s/Draw3Particle_EtaNCorr2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta2->SaveAs(name);

    cEta3=new TCanvas("cEta3","3-Particle lego",400,300);
    SetMargins2DSurf(cEta3);
    hEtaEtaCorr->SetTitle("");
    hEtaEtaCorr->Draw("lego2");
    sprintf(name,"%s/Draw3Particle_EtaEtaCorrLego_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta3->SaveAs(name);

    cEta4=new TCanvas("cEta4","3-Particle",400,350);
    SetMargins2D(cEta4);
    hEtaEtaCorr->SetTitle("");
    hEtaEtaCorr->Draw("colz");
    sprintf(name,"%s/Draw3Particle_EtaEtaCorr_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta4->SaveAs(name);

    cEta5=new TCanvas("cEta5","HS lego",400,300);
    SetMargins2DSurf(cEta5);
    hEtaEtaHSCorr->SetTitle("");
    hEtaEtaHSCorr->Draw("lego2");
    sprintf(name,"%s/Draw3Particle_EtaEtaHSCorrLego_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta5->SaveAs(name);

    cEta6=new TCanvas("cEta6","HS",400,350);
    SetMargins2D(cEta6);
    hEtaEtaHSCorr->SetTitle("");
    hEtaEtaHSCorr->Draw("colz");
    sprintf(name,"%s/Draw3Particle_EtaEtaHSCorr_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta6->SaveAs(name);    

     cEta7=new TCanvas("cEta3","SS lego",400,300);
    SetMargins2DSurf(cEta7);
    hEtaEtaSSCorr->SetTitle("");
    hEtaEtaSSCorr->Draw("lego2");
    sprintf(name,"%s/Draw3Particle_EtaEtaSSCorrLego_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta7->SaveAs(name);

    cEta8=new TCanvas("cEta8","SS",400,350);
    SetMargins2D(cEta8);
    hEtaEtaSSCorr->SetTitle("");
    hEtaEtaSSCorr->Draw("colz");
    sprintf(name,"%s/Draw3Particle_EtaEtaSSCorr_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta8->SaveAs(name);

    cEta9=new TCanvas("cEta9","3-Particle BGSub lego",400,300);
    SetMargins2DSurf(cEta9);
    hEtaEtaCorr2->SetTitle("");
    hEtaEtaCorr2->Draw("lego2");
    sprintf(name,"%s/Draw3Particle_EtaEtaCorr2Lego_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta9->SaveAs(name);

    cEta10=new TCanvas("cEta10","3-Particle BGSub",400,350);
    SetMargins2D(cEta10);
    hEtaEtaCorr2->SetTitle("");
    hEtaEtaCorr2->Draw("colz");
    sprintf(name,"%s/Draw3Particle_EtaEtaCorr2_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta10->SaveAs(name);

     cEta11=new TCanvas("cEta11","Near Projection",400,300);
     SetMargins1D(cEta11);
     hEtaEtaNOff->SetTitle("");
     hEtaEtaNOff->Draw();
    hEtaEtaNOn->Draw("same");
    keySymbol(0.2,0.85,"On",4,25,0.05,1.5);
    keySymbol(0.2,0.8,"Off",2,20,0.05,1.5);
    sprintf(name,"%s/Draw3Particle_EtaProjN_%2.1fPt%2.1f_%1.2fpt%1.2f_%dC%d_E%d_Z%1.2f%s",Folder,TPt1,TPt2,APt1,APt2,Cent1,Cent2,EffMethod,ZYACent,fileType);
    cEta11->SaveAs(name);
    
   
  }
}
