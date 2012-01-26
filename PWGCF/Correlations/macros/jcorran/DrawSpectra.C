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
//Draws jet-spectra,pt dependences, centrality dependcies, etc for output of AliAnalysisTaskDiHardron
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

void DrawSpectra(){
  gROOT->Reset();
  gROOT->LoadMacro("Util9.C");
  Style(1);
  gStyle->SetOptFit(1);
  //Parameters to be changed
  Int_t APtTPtMult=0;//0 Associated particle spectra, 1 Trigger particle dependence, 2 Centrality dependence
  Float_t TPt1=4;
  Float_t TPt2=6;
  Float_t APt1=0.25;
  Float_t APt2=2;
  Float_t Mult1=0;
  Float_t Mult2=0;
  Int_t PtWeighted=0;
  Int_t HorizontalErrors=0;
  Int_t NoTitle=1;
  Int_t DrawFit=1;
  Int_t LSign=0;
  Int_t SaveFits=1;
  int EffMethod=0;//0 no correction, 1 angle dependent, 2 use mix for triggered, 3 <1> 
  char *filetype=".gif"; //.gif, .eps etc...
  // const Int_t nTPt=12;
  //Float_t TPtBins[(nTPt+1)]={2.5,3,4,6,8,10,15,20,30,40,50,75,100};
  /*
  const Int_t nTPt=6;
  Float_t TPtBins[(nTPt+1)]={2,2.5,3,4,6,10,20};
  */
  
  //  const Int_t nTPt=11;
  //Float_t TPtBins[(nTPt+1)]={2,2.5,3,4,5,6,8,10,15,20,30,50};
   const Int_t nTPt=9;
   Float_t TPtBins[(nTPt+1)]={2,2.5,3,4,5,6,8,10,15,20};
   // const Int_t nAPt=30;
   // Float_t APtBins[(nAPt+1)]={0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,15,20,25,30,35,40,45,50,60,70,80,90,100};
  // const Int_t nAPt=17;//up to 15-20 for 7 TeV
  // Float_t APtBins[(nAPt+1)]={0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,6,8,10,15,20,30};
  // const Int_t nAPt=4;//from  20-30 for 7 TeV
  //Float_t APtBins[(nAPt+1)]={0.25,1,2,4,30};
  //const Int_t nAPt=9;
  // Float_t APtBins[(nAPt+1)]={0.25,0.5,1,1.5,2,2.5,3,3.5,4,6};
  const Int_t nAPt=6;
  Float_t APtBins[(nAPt+1)]={0.25,0.5,0.75,1,2,3,4};
  const Int_t nMult=2;
  Float_t MultBins1[nMult]={1,2};
  Float_t MultBins2[nMult]={1,2};
  Float_t MultArray1[3]={0,0,20};//>=
  Float_t MultArray2[3]={500,20,500};//<
  Float_t ZYAMCent=1.5;
  Float_t ZYAMWidth=0.2;
  float Pi=3.1415926535898;
  Float_t NearWidthPhi=1.5;
  Float_t AwayWidthPhi=Pi-NearWidthPhi;
  Float_t NearWidthEta=0.5;

  // char *Folder="2009-09-30/v4-16-Rev-08/pdc1_resdb_100_0000";
  // char *EffFolder="2009-09-30/v4-16-Rev-08/pdc1_resdb_100_0000";
  //char *Folder="current/7";
  //  char *EffFolder="current/7";
  //char *Folder="2009-12-15/RealData";//LHC09b8
  // char *EffFolder="2009-12-03/Combined900";
  // char *Folder="2010-01-06/900Data";//LHC09b8
  // char *Folder="2010-03-10/pass5";
  // char *EffFolder="2010-03-10/all_5";
  // char *Folder="2010-04-12/pass5";
  //  char *EffFolder="2010-04-12/900Pythia";
  /*
  char *Folder="2010-08-17/LHC10b_7pass2";
  char *EffFolder="2010-08-17/7Pythia_LHC10b5";
  */
  
  char *Folder="2010-08-17/LHC10c_900pass2";
  //char *Folder="2010-07-02/LHC10c6_900Pythia";
   char *EffFolder="2010-08-17/LHC10c6_900Pythia";
  

 //Plotting 
  Float_t MarkerSize=1.5;
  Int_t ColorNearPhi=kBlack;
  Int_t ColorAwayPhi=kRed;
  Int_t ColorMin=kGreen+2;
  Int_t ColorNearPhiFit=kBlue+2;
  Int_t ColorAwayPhiFit=kRed+2;
  Int_t ColorMinFit=kGreen+3;
  Int_t ColorNearEta=kGreen+2;
  Int_t ColorAwayEta=kMagenta+3;
  Int_t ColorNearEtaFit=kGreen+3;
  Int_t MarkerNearPhi=20;
  Int_t MarkerAwayPhi=21;
  Int_t MarkerMin=22;
  Int_t MarkerNearPhiFit=24;
  Int_t MarkerAwayPhiFit=25;
  Int_t MarkerMinFit=26;
  Int_t MarkerNearEta=29;
  Int_t MarkerAwayEta=23;
  Int_t MarkerNearEtaFit=30;
  
  
  TF1 *fit1=new TF1("fit1","1/sqrt(2*3.1415926)*([0]/[1]*(exp(-0.5*pow(x/[1],2))+exp(-0.5*pow((x-6.29185)/[1],2)))+[2]/[3]*(exp(-0.5*pow((x-3.14159)/[3],2))+exp(-0.5*pow((x+3.14159)/[3],2))))+[4]");
  fit1->SetParameters(.1,.2,.1,.35,.1);
  fit1->SetParNames("Near-Yield", "Near-Width","Away-Yield","Away-Width","Bg-Level");
  fit1->SetParLimits(0,1E-10,100);
  fit1->SetParLimits(1,0.01,6);
  fit1->SetParLimits(2,1E-10,100);
  fit1->SetParLimits(3,0.01,6);
  fit1->SetParLimits(4,1E-10,500);
   TF1 *fit2=new TF1("fit2","1/sqrt(2*3.1415926)*([0]/[1]*(exp(-0.5*pow(x/[1],2))))+[2]");
  fit2->SetParameters(.1,.2,.1);
  fit2->SetParNames("Near-Yield", "Near-Width","Bg-Level");
  fit2->SetParLimits(0,1E-10,100);
  fit2->SetParLimits(1,0.01,6);
  fit2->SetParLimits(2,1E-10,500);
  //Bollow this users should not need to change
  char *cPt[2]={"","PT"};
  char name[250];
  char inName[100];
  char effName[100];
  sprintf(inName,"%s/julery_DiHadron.root",Folder);
  sprintf(effName,"%s/julery_DiHadron.root",EffFolder);
  cout << inName << endl;
  cout << effName << endl;
  Float_t MPt[3];//<pt>, error <pt>, # of triggers
  Float_t MPt2[3];
  Float_t TrigSum=0;
  
  TH1F *hPhiRaw=new TH1F("hPhiRaw","",1,0,1);
  TH1F *hPhiCorr=new TH1F("hPhiCorr","",1,0,1);
  TH1F *hPhiEff=new TH1F("hPhiEff","",1,0,1);
  TH1F *hPhiMC=new TH1F("hPhiMC","",1,0,1);
  TH1F *hPhiMixRaw=new TH1F("hPhiMixRaw","",1,0,1);
  TH1F *hPhiMixCorr=new TH1F("hPhiMixCorr","",1,0,1);
  TH1F *hPhiMixEff=new TH1F("hPhiMixEff","",1,0,1);
  TH1F *hPhiMixMC=new TH1F("hPhiMixMC","",1,0,1);
  TH2F *hPhiEtaRaw=new TH2F("hPhiEtaRaw","",1,0,1,1,0,1);
  TH2F *hPhiEtaCorr=new TH2F("hPhiEtaCorr","",1,0,1,1,0,1);
  TH2F *hPhiEtaEff=new TH2F("hPhiEtaEff","",1,0,1,1,0,1);
  TH2F *hPhiEtaMC=new TH2F("hPhiEtaMC","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixRaw=new TH2F("hPhiEtaMixRaw","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixCorr=new TH2F("hPhiEtaMixCorr","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixEff=new TH2F("hPhiEtaMixEff","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixMC=new TH2F("hPhiEtaMixMC","",1,0,1,1,0,1);

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
  
  TH1F *hTPhiRaw=new TH1F("hTPhiRaw","For near-side pt sum",1,0,1);
  TH1F *hTPhiCorr=new TH1F("hTPhiCorr","",1,0,1);
  TH1F *hTPhiEff=new TH1F("hTPhiErr","",1,0,1);
  TH1F *hTPhiMC=new TH1F("hTPhiMC","",1,0,1);
  TH1F *hTPhiMixRaw=new TH1F("hTPhiMixRaw","",1,0,1);
  TH1F *hTPhiMixCorr=new TH1F("hTPhiMixCorr","",1,0,1);
  TH1F *hTPhiMixEff=new TH1F("hTPhiMixEff","",1,0,1);
  TH1F *hTPhiMixMC=new TH1F("hTPhiMixMC","",1,0,1);
  TH2F *hTPhiEtaRaw=new TH2F("hTPhiEtaRaw","",1,0,1,1,0,1);
  TH2F *hTPhiEtaCorr=new TH2F("hTPhiEtaCorr","",1,0,1,1,0,1);
  TH2F *hTPhiEtaEff=new TH2F("hTPhiEtaEff","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMC=new TH2F("hTPhiEtaMC","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMixRaw=new TH2F("hTPhiEtaMixRaw","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMixCorr=new TH2F("hTPhiEtaMixCorr","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMixEff=new TH2F("hTPhiEtaMixEff","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMixMC=new TH2F("hTPhiEtaMixMC","",1,0,1,1,0,1);

  TH1F *hTEtaNRaw=new TH1F("hTEtaNRaw","",1,0,1);
  TH1F *hTEtaNCorr=new TH1F("hTEtaNCorr","",1,0,1);
  TH1F *hTEtaNEff=new TH1F("hTEtaNEff","",1,0,1);
  TH1F *hTEtaNMC=new TH1F("hTEtaNMC","",1,0,1);
  TH1F *hTEtaNMixRaw=new TH1F("hTEtaNMixRaw","",1,0,1);
  TH1F *hTEtaNMixCorr=new TH1F("hTEtaNMixCorr","",1,0,1);
  TH1F *hTEtaNMixEff=new TH1F("hTEtaNMixEff","",1,0,1);
  TH1F *hTEtaNMixMC=new TH1F("hTEtaNMixMC","",1,0,1);

  TH1F *hTEtaARaw=new TH1F("hTEtaARaw","",1,0,1);
  TH1F *hTEtaACorr=new TH1F("hTEtaACorr","",1,0,1);
  TH1F *hTEtaAEff=new TH1F("hTEtaAEff","",1,0,1);
  TH1F *hTEtaAMC=new TH1F("hTEtaAMC","",1,0,1);
  TH1F *hTEtaAMixRaw=new TH1F("hTEtaAMixRaw","",1,0,1);
  TH1F *hTEtaAMixCorr=new TH1F("hTEtaAMixCorr","",1,0,1);
  TH1F *hTEtaAMixEff=new TH1F("hTEtaAMixEff","",1,0,1);
  TH1F *hTEtaAMixMC=new TH1F("hTEtaAMixMC","",1,0,1);
  
  TH1F *hPhiRawPt=new TH1F("hPhiRawPt","",1,0,1);
  TH1F *hPhiCorrPt=new TH1F("hPhiCorrPt","",1,0,1);
  TH1F *hPhiEffPt=new TH1F("hPhiEffPt","",1,0,1);
  TH1F *hPhiMCPt=new TH1F("hPhiMCPt","",1,0,1);
  TH1F *hPhiMixRawPt=new TH1F("hPhiMixRawPt","",1,0,1);
  TH1F *hPhiMixCorrPt=new TH1F("hPhiMixCorrPt","",1,0,1);
  TH1F *hPhiMixEffPt=new TH1F("hPhiMixEffPt","",1,0,1);
  TH1F *hPhiMixMCPt=new TH1F("hPhiMixMCPt","",1,0,1);
  TH2F *hPhiEtaRawPt=new TH2F("hPhiEtaRawPt","",1,0,1,1,0,1);
  TH2F *hPhiEtaCorrPt=new TH2F("hPhiEtaCorrPt","",1,0,1,1,0,1);
  TH2F *hPhiEtaEffPt=new TH2F("hPhiEtaEffPt","",1,0,1,1,0,1);
  TH2F *hPhiEtaMCPt=new TH2F("hPhiEtaMCPt","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixRawPt=new TH2F("hPhiEtaMixRawPt","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixCorrPt=new TH2F("hPhiEtaMixCorrPt","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixEffPt=new TH2F("hPhiEtaMixEffPt","",1,0,1,1,0,1);
  TH2F *hPhiEtaMixMCPt=new TH2F("hPhiEtaMixMCPt","",1,0,1,1,0,1);

  TH1F *hEtaNRawPt=new TH1F("hEtaNRawPt","",1,0,1);
  TH1F *hEtaNCorrPt=new TH1F("hEtaNCorrPt","",1,0,1);
  TH1F *hEtaNEffPt=new TH1F("hEtaNEffPt","",1,0,1);
  TH1F *hEtaNMCPt=new TH1F("hEtaNMCPt","",1,0,1);
  TH1F *hEtaNMixRawPt=new TH1F("hEtaNMixRawPt","",1,0,1);
  TH1F *hEtaNMixCorrPt=new TH1F("hEtaNMixCorrPt","",1,0,1);
  TH1F *hEtaNMixEffPt=new TH1F("hEtaNMixEffPt","",1,0,1);
  TH1F *hEtaNMixMCPt=new TH1F("hEtaNMixMCPt","",1,0,1);

  TH1F *hEtaARawPt=new TH1F("hEtaARawPt","",1,0,1);
  TH1F *hEtaACorrPt=new TH1F("hEtaACorrPt","",1,0,1);
  TH1F *hEtaAEffPt=new TH1F("hEtaAEffPt","",1,0,1);
  TH1F *hEtaAMCPt=new TH1F("hEtaAMCPt","",1,0,1);
  TH1F *hEtaAMixRawPt=new TH1F("hEtaAMixRawPt","",1,0,1);
  TH1F *hEtaAMixCorrPt=new TH1F("hEtaAMixCorrPt","",1,0,1);
  TH1F *hEtaAMixEffPt=new TH1F("hEtaAMixEffPt","",1,0,1);
  TH1F *hEtaAMixMCPt=new TH1F("hEtaAMixMCPt","",1,0,1);
  
  TH1F *hTPhiRawPt=new TH1F("hTPhiRawPt","",1,0,1);
  TH1F *hTPhiCorrPt=new TH1F("hTPhiCorrPt","",1,0,1);
  TH1F *hTPhiEffPt=new TH1F("hTPhiEffPt","",1,0,1);
  TH1F *hTPhiMCPt=new TH1F("hTPhiMCPt","",1,0,1);
  TH1F *hTPhiMixRawPt=new TH1F("hTPhiMixRawPt","",1,0,1);
  TH1F *hTPhiMixCorrPt=new TH1F("hTPhiMixCorrPt","",1,0,1);
  TH1F *hTPhiMixEffPt=new TH1F("hTPhiMixEffPt","",1,0,1);
  TH1F *hTPhiMixMCPt=new TH1F("hTPhiMixMCPt","",1,0,1);
  TH2F *hTPhiEtaRawPt=new TH2F("hTPhiEtaRawPt","",1,0,1,1,0,1);
  TH2F *hTPhiEtaCorrPt=new TH2F("hTPhiEtaCorrPt","",1,0,1,1,0,1);
  TH2F *hTPhiEtaEffPt=new TH2F("hTPhiEtaEffPt","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMCPt=new TH2F("hTPhiEtaMCPt","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMixRawPt=new TH2F("hTPhiEtaMixRawPt","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMixCorrPt=new TH2F("hTPhiEtaMixCorrPt","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMixEffPt=new TH2F("hTPhiEtaMixEffPt","",1,0,1,1,0,1);
  TH2F *hTPhiEtaMixMCPt=new TH2F("hTPhiEtaMixMCPt","",1,0,1,1,0,1);

  TH1F *hTEtaNRawPt=new TH1F("hTEtaNRawPt","",1,0,1);
  TH1F *hTEtaNCorrPt=new TH1F("hTEtaNCorrPt","",1,0,1);
  TH1F *hTEtaNEffPt=new TH1F("hTEtaNEffPt","",1,0,1);
  TH1F *hTEtaNMCPt=new TH1F("hTEtaNMCPt","",1,0,1);
  TH1F *hTEtaNMixRawPt=new TH1F("hTEtaNMixRawPt","",1,0,1);
  TH1F *hTEtaNMixCorrPt=new TH1F("hTEtaNMixCorrPt","",1,0,1);
  TH1F *hTEtaNMixEffPt=new TH1F("hTEtaNMixEffPt","",1,0,1);
  TH1F *hTEtaNMixMCPt=new TH1F("hTEtaNMixMCPt","",1,0,1);

  TH1F *hTEtaARawPt=new TH1F("hTEtaARawPt","",1,0,1);
  TH1F *hTEtaACorrPt=new TH1F("hTEtaACorrPt","",1,0,1);
  TH1F *hTEtaAEffPt=new TH1F("hTEtaAEffPt","",1,0,1);
  TH1F *hTEtaAMCPt=new TH1F("hTEtaAMCPt","",1,0,1);
  TH1F *hTEtaAMixRawPt=new TH1F("hTEtaAMixRawPt","",1,0,1);
  TH1F *hTEtaAMixCorrPt=new TH1F("hTEtaAMixCorrPt","",1,0,1);
  TH1F *hTEtaAMixEffPt=new TH1F("hTEtaAMixEffPt","",1,0,1);
  TH1F *hTEtaAMixMCPt=new TH1F("hTEtaAMixMCPt","",1,0,1);
  
  TH1F *hMult=new TH1F("hMult","",1,0,1);
  TH1F *hMultTrig=new TH1F("hMultTrig","",1,0,1);
  TH1F *hMult2=new TH1F("hMult2","",1,0,1);
  TH1F *hMultTrig2=new TH1F("hMultTrig2","",1,0,1);
  float plotmin=10,plotmax=11;
  float plotmin2=10,plotmax2=11;
  float plotmin3=10,plotmax3=11;
  
  //histograms for ZYAM
  TH1F *hPhiCorr2;
  TH1F *hPhiMixCorr2;
  TH1F *hEtaNCorr2;
  TH1F *hEtaACorr2;
  TH1F *hEtaMixCorr2;
  TH2F *hPhiEtaCorr2;
  TH2F *hPhiEtaMixCorr2;
  
  TH1F *hPhiCorrPt2;
  TH1F *hPhiMixCorrPt2;
  TH1F *hEtaNCorrPt2;
  TH1F *hEtaACorrPt2;
  TH1F *hEtaMixCorr2;
  TH2F *hPhiEtaCorrPt2;
  TH2F *hPhiEtaMixCorrPt2;
  
  //histograms for summing the centraility bins
  TH1F *hPhiCorrSum;
  TH1F *hPhiCorr2Sum;
  TH1F *hPhiMixCorrSum;
  TH1F *hPhiMixCorr2Sum;
  TH1F *hEtaNCorrSum;
  TH1F *hEtaNCorr2Sum;
  TH1F *hEtaNMixCorrSum;
  TH1F *hEtaNMixCorr2Sum;
  TH1F *hEtaACorrSum;
  TH1F *hEtaACorr2Sum; 
  TH1F *hEtaAMixCorrSum;
  TH1F *hEtaAMixCorr2Sum;
  TH2F *hPhiEtaCorrSum;
  TH2F *hPhiEtaCorr2Sum;
  TH2F *hPhiEtaMixCorrSum;
  TH2F *hPhiEtaMixCorr2Sum;

  TH1F *hPhiCorrPtSum;
  TH1F *hPhiCorrPt2Sum;
  TH1F *hPhiMixCorrPtSum;
  TH1F *hPhiMixCorrPt2Sum;
  TH1F *hEtaNCorrPtSum;
  TH1F *hEtaNCorrPt2Sum;
  TH1F *hEtaNMixCorrPtSum;
  TH1F *hEtaNMixCorrPt2Sum;
  TH1F *hEtaACorrPtSum;
  TH1F *hEtaACorrPt2Sum;
  TH1F *hEtaAMixCorrPtSum;
  TH1F *hEtaAMixCorrPt2Sum;
  TH2F *hPhiEtaCorrPtSum;
  TH2F *hPhiEtaCorrPt2Sum;
  TH2F *hPhiEtaMixCorrPtSum;
  TH2F *hPhiEtaMixCorrPt2Sum;
  TH1F *hTPhiCorrPt2Sum;
  TH2F *hTPhiEtaCorrSum;
  TH2F *hTPhiEtaCorrPt2Sum;
  
  int NearBinsPhi;
  int AwayBinsPhi;
  int MinBinsPhi;
  int NearBinsEta;
  int AwayBinsEta;
  int NearBinsPhiEta1;
  int NearBinsPhiEta2;
  int AwayBinsPhiEta;
  int MinBinsPhiEta;
  float tempsum;
  
  TAxis *phiaxis;
  TAxis *etaaxis;
  
  int loopSize=0;
  int nAPt2=0;
  for(int i=0;i<=nAPt;i++){if(fabs(TPt1-APtBins[i])<0.2)nAPt2=i;}
  if(APtTPtMult==0)loopSize=nAPt2;
  if(APtTPtMult==1)loopSize=nTPt;
  if(APtTPtMult==2)loopSize=nMult;
  const int nMainLoop=loopSize;
  cout << "loopSize " << loopSize << endl;
  if(loopSize==0){
    cout << "Error: Check Arrays" << endl;
    break;
  }
  Float_t MainArray[nMainLoop];
  Float_t eMainArray[nMainLoop];
  Float_t a[3];
  float centx, centy;
  float binwx, binwy;
  const int nx=18;
  //0 (x1+x2)/2, 1 e2, 2 <x>,  3 e2, 4 lower error of 1 and 2, 5 e3, 6 zt1_0, 7 e5, 8 zt1_2, 9 e7, 10 zt1_3, 11 e9, 12 zt2_0, 13 e11, 14 zt2_1, 15 e13, 16 zt2_2, 17 e15  (zt_0 is (pt1+pt2)/2 zt_1 is <pt> zt_2 is lower error)
  float NXPhi[nx][nMainLoop];
  float AXPhi[nx][nMainLoop];
  float NXPhiFit[nx][nMainLoop];
  float AXPhiFit[nx][nMainLoop];
  float NXEta[nx][nMainLoop];
  float AXEta[nx][nMainLoop];
  float NXEtaFit[nx][nMainLoop];
  float AXEtaFit[nx][nMainLoop];
  float NXPhiEta1[nx][nMainLoop];
  float NXPhiEta2[nx][nMainLoop];
  float AXPhiEta[nx][nMainLoop];
  //arrays for our tgraphs to store yields and averge bin contends
  const int nyields=10;//0 yield, 1 eyield, 2 pt yield, 3 ept yield, 4 dN/dpt 5 e4 6 dN/dzt (zt=Trig) yield, 7 ezt yield, 8 dN/dzt (zt=near) yield, 9 ezt2 yield
  float NYieldPhi[nyields][nMainLoop];
  float NYieldPhiZYAM[nyields][nMainLoop];
  float NYieldPhiFit[nyields][nMainLoop];
  float NYieldEta[nyields][nMainLoop];
  float NYieldEtaZYAM[nyields][nMainLoop];
  float NYieldEtaFit[nyields][nMainLoop];
  float NYieldPhiEta1[nyields][nMainLoop];
  float NYieldPhiEta2[nyields][nMainLoop];
  float NYieldPhiEta3[nyields][nMainLoop];//jet only
  float NYieldPhiEta1ZYAM[nyields][nMainLoop];//Jet+Ridge
  float NYieldPhiEta2ZYAM[nyields][nMainLoop];//Ridge
  float NRMSPhi[nyields][nMainLoop];
  float NRMSEta[nyields][nMainLoop];
  float NWidthPhi[nyields][nMainLoop];
  float NWidthEta[nyields][nMainLoop];
  float NJtPhi[nyields][nMainLoop];
  float NJtEta[nyields][nMainLoop];
  float NRMSPhiEtaPhi[nyields][nMainLoop];
  float NRMSPhiEtaEta[nyields][nMainLoop];
  float AYieldPhi[nyields][nMainLoop];
  float AYieldPhiZYAM[nyields][nMainLoop];
  float AYieldPhiFit[nyields][nMainLoop];
  float AYieldEta[nyields][nMainLoop];
  float AYieldEtaZYAM[nyields][nMainLoop];
  float AYieldEtaFit[nyields][nMainLoop];
  float AYieldPhiEta[nyields][nMainLoop];
  float AYieldPhiEtaZYAM[nyields][nMainLoop];
  float ARMSPhi[nyields][nMainLoop];
  float ARMSEta[nyields][nMainLoop];
  float AWidthPhi[nyields][nMainLoop];
  float AWidthEta[nyields][nMainLoop];
  float AJtPhi[nyields][nMainLoop];
  float AJtEta[nyields][nMainLoop];
  float MYieldPhi[nyields][nMainLoop];
  float MYieldPhiFit[nyields][nMainLoop];
  float MYieldPhiEta[nyields][nMainLoop];
  float MAvePhiEta[nyields][nMainLoop];

  float TNYieldPhi[nyields][nMainLoop];
  float TNYieldPhiZYAM[nyields][nMainLoop];
  float TNYieldPhiFit[nyields][nMainLoop];
  float TNYieldPhiEta1[nyields][nMainLoop];
  float TNYieldPhiEta2[nyields][nMainLoop];
  float TNYieldPhiEta3[nyields][nMainLoop];
  float TAYieldPhi[nyields][nMainLoop];
  float TAYieldPhiZYAM[nyields][nMainLoop];
  float TAYieldPhiFit[nyields][nMainLoop];
  float TAYieldPhiEta[nyields][nMainLoop];
  float TAYieldPhiEtaZYAM[nyields][nMainLoop];

  float NAvePhi[nyields][nMainLoop];
  float NAvePhiZYAM[nyields][nMainLoop];
  float NAvePhiFit[nyields][nMainLoop];
  float NAveEta[nyields][nMainLoop];
  float NAveEtaZYAM[nyields][nMainLoop];
  float NAveEtaFit[nyields][nMainLoop];
  float NAvePhiEta1[nyields][nMainLoop];
  float NAvePhiEta2[nyields][nMainLoop];
  float NAvePhiEta1ZYAM[nyields][nMainLoop];
  float NAvePhiEta2ZYAM[nyields][nMainLoop];
  float AAvePhi[nyields][nMainLoop];
  float AAvePhiZYAM[nyields][nMainLoop];
  float AAvePhiFit[nyields][nMainLoop];
  float AAveEta[nyields][nMainLoop];
  float AAveEtaZYAM[nyields][nMainLoop];
  float AAveEtaFit[nyields][nMainLoop];
  float AAvePhiEta[nyields][nMainLoop];
  float AAvePhiEtaZYAM[nyields][nMainLoop];
  float MAvePhi[nyields][nMainLoop];
  float MAvePhiFit[nyields][nMainLoop];
  float MAvePhiEta[nyields][nMainLoop];
  float trms1;
  float trms2;
  float trms3;
  TFile *inFile=new TFile(inName);
  TList *inList=inFile->Get("julery_DiHadron");
  TFile *effFile=new TFile(effName);
  TList *effList=effFile->Get("julery_DiHadron");
  char *FitTit[2]={"","_Fit"};
  //loop for total near-side yield over full pt range
  for(int Cent=Mult1;Cent<=Mult2;Cent++){
 MakeProjections(TPt1,TPt2,0.25,TPt1,Cent,inList,hTPhiRaw,hTPhiMixRaw,hTEtaNRaw,hTEtaNMixRaw,hTEtaARaw,hTEtaAMixRaw,hTPhiEtaRaw,hTPhiEtaMixRaw,MPt,1,0,0,0);
    MakeProjections(TPt1,TPt2,0.25,TPt1,Cent,effList,hTPhiEff,hTPhiMixEff,hTEtaNEff,hTEtaNMixEff,hTEtaAEff,hTEtaAMixEff,hTPhiEtaEff,hTPhiEtaMixEff,MPt2,0,0,0,0);
    MakeProjections(TPt1,TPt2,0.25,TPt1,Cent,effList,hTPhiMC,hTPhiMixMC,hTEtaNMC,hTEtaNMixMC,hTEtaAMC,hTEtaAMixMC,hTPhiEtaMC,hTPhiEtaMixMC,MPt2,0,1,0,0);
    MakeProjections(TPt1,TPt2,0.25,TPt1,Cent,inList,hTPhiRawPt,hTPhiMixRawPt,hTEtaNRawPt,hTEtaNMixRawPt,hTEtaARawPt,hTEtaAMixRawPt,hTPhiEtaRawPt,hTPhiEtaMixRawPt,MPt,1,0,0,0);
    MakeProjections(TPt1,TPt2,0.25,TPt1,Cent,effList,hTPhiEffPt,hTPhiMixEffPt,hTEtaNEffPt,hTEtaNMixEffPt,hTEtaAEffPt,hTEtaAMixEffPt,hTPhiEtaEffPt,hTPhiEtaMixEffPt,MPt2,1,0,0,0);
    MakeProjections(TPt1,TPt2,0.25,TPt1,Cent,effList,hTPhiMCPt,hTPhiMixMCPt,hTEtaNMCPt,hTEtaNMixMCPt,hTEtaAMCPt,hTEtaAMixMCPt,hTPhiEtaMCPt,hTPhiEtaMixMCPt,MPt2,1,1,0,0);
    // cout << "hTPhiRaw " << hTPhiRaw->GetBinContent(1) << endl;
 if(EffMethod<4){
      EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hTPhiEff,hTPhiMC,hTPhiMixEff,hTPhiMixMC,hTEtaNEff,hTEtaNMC,hTEtaNMixEff,hTEtaNMixMC,hTEtaAEff,hTEtaAMC,hTEtaAMixEff,hTEtaAMixMC,hTPhiEtaEff,hTPhiEtaMC,hTPhiEtaMixEff,hTPhiEtaMixMC,EffMethod);
      EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hTPhiEffPt,hTPhiMCPt,hTPhiMixEffPt,hTPhiMixMCPt,hTEtaNEffPt,hTEtaNMCPt,hTEtaNMixEffPt,hTEtaNMixMCPt,hTEtaAEffPt,hTEtaAMCPt,hTEtaAMixEffPt,hTEtaAMixMCPt,hTPhiEtaEffPt,hTPhiEtaMCPt,hTPhiEtaMixEffPt,hTPhiEtaMixMCPt,EffMethod);
    }
    else{
      EffFit(APt1,APt2,Cent,effList,hTPhiEff,hTPhiMixEff,hTEtaNEff,hTEtaNMixEff,hTEtaAEff,hTEtaAMixEff,hTPhiEtaEff,hTPhiEtaMixEff,LSign,VariablePtLimit);
      EffFit(APt1,APt2,Cent,effList,hTPhiEffPt,hTPhiMixEffPt,hTEtaNEffPt,hTEtaNMixEffPt,hTEtaAEffPt,hTEtaAMixEffPt,hTPhiEtaEffPt,hTPhiEtaMixEffPt,LSign,VariablePtLimit);
    }
    /*
    EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hTPhiRaw,hTPhiEff,hTPhiMC,hTPhiMixRaw,hTPhiMixEff,hTPhiMixMC,hTEtaNRaw,hTEtaNEff,hTEtaNMC,hTEtaNMixRaw,hTEtaNMixEff,hTEtaNMixMC,hTEtaARaw,hTEtaAEff,hTEtaAMC,hTEtaAMixRaw,hTEtaAMixEff,hTEtaAMixMC,hTPhiEtaRaw,hTPhiEtaEff,hTPhiEtaMC,hTPhiEtaMixRaw,hTPhiEtaMixEff,hTPhiEtaMixMC,EffMethod);
    EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hTPhiRawPt,hTPhiEffPt,hTPhiMCPt,hTPhiMixRawPt,hTPhiMixEffPt,hTPhiMixMCPt,hTEtaNRawPt,hTEtaNEffPt,hTEtaNMCPt,hTEtaNMixRawPt,hTEtaNMixEffPt,hTEtaNMixMCPt,hTEtaARawPt,hTEtaAEffPt,hTEtaAMCPt,hTEtaAMixRawPt,hTEtaAMixEffPt,hTEtaAMixMCPt,hTPhiEtaRawPt,hTPhiEtaEffPt,hTPhiEtaMCPt,hTPhiEtaMixRawPt,hTPhiEtaMixEffPt,hTPhiEtaMixMCPt,EffMethod);
    */
    
    hTPhiCorr=(TH1F*)hTPhiRaw->Clone();
    sprintf(name,"hTPhiCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTPhiCorr->SetName(name);
    hTPhiCorr->Divide(hTPhiEff);
    
    hTPhiMixCorr=(TH1F*)hTPhiMixRaw->Clone();
    sprintf(name,"hTPhiMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTPhiMixCorr->SetName(name);
    hTPhiMixCorr->Divide(hTPhiMixEff);

    hTEtaNCorr=(TH1F*)hTEtaNRaw->Clone();
    sprintf(name,"hTEtaNCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTEtaNCorr->SetName(name);
    hTEtaNCorr->Divide(hTEtaNEff);
    
    hTEtaNMixCorr=(TH1F*)hTEtaNMixRaw->Clone();
    sprintf(name,"hTEtaNMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTEtaNMixCorr->SetName(name);
    hTEtaNMixCorr->Divide(hTPhiMixEff);

    hTEtaACorr=(TH1F*)hTEtaARaw->Clone();
    sprintf(name,"hTEtaACorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTEtaACorr->SetName(name);
    hTEtaACorr->Divide(hTEtaAEff);
    
    hTEtaAMixCorr=(TH1F*)hTEtaAMixRaw->Clone();
    sprintf(name,"hTEtaAMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTEtaAMixCorr->SetName(name);
    hTEtaAMixCorr->Divide(hTPhiMixEff);
    
    hTPhiEtaCorr=(TH2F*)hTPhiEtaRaw->Clone();
    sprintf(name,"hTPhiEtaCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTPhiEtaCorr->SetName(name);
    hTPhiEtaCorr->Divide(hTPhiEtaEff);
    
    hTPhiEtaMixCorr=(TH2F*)hTPhiEtaMixRaw->Clone();
    sprintf(name,"hTPhiMixCorr_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTPhiEtaMixCorr->SetName(name);
    hTPhiEtaMixCorr->Divide(hTPhiEtaMixEff);
    
    hTPhiCorrPt=(TH1F*)hTPhiRawPt->Clone();
    sprintf(name,"hTPhiCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTPhiCorrPt->SetName(name);
    hTPhiCorrPt->Divide(hTPhiEffPt);
    
    hTPhiMixCorrPt=(TH1F*)hTPhiMixRawPt->Clone();
    sprintf(name,"hTPhiMixCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTPhiMixCorrPt->SetName(name);
    hTPhiMixCorrPt->Divide(hTPhiMixEffPt);

    hTEtaNCorrPt=(TH1F*)hTEtaNRawPt->Clone();
    sprintf(name,"hTEtaNCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTEtaNCorrPt->SetName(name);
    hTEtaNCorrPt->Divide(hTEtaNEffPt);
    
    hTEtaNMixCorrPt=(TH1F*)hTEtaNMixRawPt->Clone();
    sprintf(name,"hTEtaNMixCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTEtaNMixCorrPt->SetName(name);
    hTEtaNMixCorrPt->Divide(hTEtaNMixEffPt);

    hTEtaACorrPt=(TH1F*)hTEtaARawPt->Clone();
    sprintf(name,"hTEtaACorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTEtaACorrPt->SetName(name);
    hTEtaACorrPt->Divide(hTEtaNEffPt);
    
    hTEtaAMixCorrPt=(TH1F*)hTEtaAMixRawPt->Clone();
    sprintf(name,"hTEtaAMixCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTEtaAMixCorrPt->SetName(name);
    hTEtaAMixCorrPt->Divide(hTEtaAMixEffPt);
    
    hTPhiEtaCorrPt=(TH2F*)hTPhiEtaRawPt->Clone();
    sprintf(name,"hTPhiEtaCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTPhiEtaCorrPt->SetName(name);
    hTPhiEtaCorrPt->Divide(hTPhiEtaEffPt);
    
    hTPhiEtaMixCorrPt=(TH2F*)hTPhiEtaMixRawPt->Clone();
    sprintf(name,"hTPhiMixCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
    hTPhiEtaMixCorrPt->SetName(name);
    hTPhiEtaMixCorrPt->Divide(hTPhiEtaMixEffPt);
    //  cout << "hTPhiCorr " << hTPhiCorr->GetBinContent(1) << endl;
    //GeoCorr(hTEtaNCorr,hTEtaNMixCorr,hTEtaNMixMC,hTEtaACorr,hTEtaAMixCorr,hTEtaAMixMC,hTPhiEtaCorr,hTPhiEtaMixCorr,hTPhiEtaMixMC);
    //GeoCorr(hTEtaNCorrPt,hTEtaNMixCorrPt,hTEtaNMixMCPt,hTEtaACorrPt,hTEtaAMixCorrPt,hTEtaAMixMCPt,hTPhiEtaCorrPt,hTPhiEtaMixCorrPt,hTPhiEtaMixMCPt);
    MixedCorrect(hTPhiCorr,hTPhiMixCorr,hTEtaNCorr,hTEtaNMixCorr,hTEtaACorr,hTEtaAMixCorr,hTPhiEtaCorr,hTPhiEtaMixCorr);
  MixedCorrect(hTPhiCorrPt,hTPhiMixCorrPt,hTEtaNCorrPt,hTEtaNMixCorrPt,hTEtaACorrPt,hTEtaAMixCorrPt,hTPhiEtaCorrPt,hTPhiEtaMixCorrPt);
    // cout << "hTPhiCorr " << hTPhiCorr->GetBinContent(1) << endl;
  if(ZYAMCent>0)ZYA1(hTPhiCorr,hTPhiMixCorr,a,ZYAMCent,ZYAMWidth);
  else ZYAM2D(hPhiEtaCorr,hPhiEtaMixCorr,a,72,10);
    cout << "Normalization Factor: " << a[0] << endl;
    hTPhiCorr2=(TH1F*)hTPhiCorr->Clone();
    hTPhiCorr2->SetName("hTPhiCorr2");
    sprintf(name,"Bg Sub. %s",hTPhiCorr2->GetTitle());
    hTPhiCorr2->SetTitle(name);
    
    hTPhiMixCorr2=(TH1F*)hTPhiMixCorr->Clone();
    hTPhiMixCorr2->SetName("hTPhiMixCorr2");
    sprintf(name,"Normalized %1.3f %s",a[0],hTPhiMixCorr2->GetTitle());
    hTPhiMixCorr2->SetTitle(name);
    hTPhiMixCorr2->SetMarkerStyle(25);
    hTPhiMixCorr2->SetMarkerColor(4);
    hTPhiMixCorr2->SetLineColor(4);
    
    hTPhiEtaCorr2=(TH2F*)hTPhiEtaCorr->Clone();
    hTPhiEtaCorr2->SetName("hTPhiEtaCorr2");
    sprintf(name,"BG Sub. %s",hTPhiEtaCorr2->GetTitle());
    hTPhiEtaCorr2->SetTitle(name);
    
    hTPhiEtaMixCorr2=(TH2F*)hTPhiEtaMixCorr->Clone();
    hTPhiEtaMixCorr2->SetName("hTPhiEtaMixCorr2");
    sprintf(name,"Normalized %1.3f %s",a[0],hTPhiEtaMixCorr2->GetTitle());
    hTPhiEtaMixCorr2->SetTitle(name);
    
    hTPhiCorrPt2=(TH1F*)hTPhiCorrPt->Clone();
    hTPhiCorrPt2->SetName("hTPhiCorrPt2");
    sprintf(name,"Bg Sub. %s",hTPhiCorrPt2->GetTitle());
    hTPhiCorrPt2->SetTitle(name);
    
    hTPhiMixCorrPt2=(TH1F*)hTPhiMixCorrPt->Clone();
    hTPhiMixCorrPt2->SetName("hTPhiMixCorrPt2");
    sprintf(name,"Normalized %1.3f %s",a[0],hTPhiMixCorrPt2->GetTitle());
    hTPhiMixCorrPt2->SetTitle(name);
    hTPhiMixCorrPt2->SetMarkerStyle(25);
    hTPhiMixCorrPt2->SetMarkerColor(4);
    hTPhiMixCorrPt2->SetLineColor(4);
    
    hTPhiEtaCorrPt2=(TH2F*)hTPhiEtaCorrPt->Clone();
    hTPhiEtaCorrPt2->SetName("hTPhiEtaCorrPt2");
    sprintf(name,"BG Sub. %s",hTPhiEtaCorrPt2->GetTitle());
    hTPhiEtaCorrPt2->SetTitle(name);
    
    hTPhiEtaMixCorrPt2=(TH2F*)hTPhiEtaMixCorrPt->Clone();
    hTPhiEtaMixCorrPt2->SetName("hTPhiEtaMixCorrPt2");
    sprintf(name,"Normalized %1.3f %s",a[0],hTPhiEtaMixCorrPt2->GetTitle());
    hTPhiEtaMixCorrPt2->SetTitle(name);
    hTPhiMixCorr2->Scale(a[0]);
    hTPhiEtaMixCorr2->Scale(a[0]);
    hTPhiCorr2->Add(hTPhiMixCorr2,-1);
    Add2D(hTPhiEtaCorr2,hTPhiEtaMixCorr2,-1); 	
    hTPhiMixCorrPt2->Scale(a[0]);
    hTPhiEtaMixCorrPt2->Scale(a[0]);
    hTPhiCorrPt2->Add(hTPhiMixCorrPt2,-1);
    Add2D(hTPhiEtaCorrPt2,hTPhiEtaMixCorrPt2,-1); 
    if(Cent==Mult1){
      hTPhiCorrSum=(TH1F*)hTPhiCorr->Clone();
      hTPhiCorrSum->SetName("hTPhiCorrSum");
      hTPhiCorr2Sum=(TH1F*)hTPhiCorr2->Clone();
      hTPhiCorr2Sum->SetName("hTPhiCorr2Sum");
      hTPhiMixCorrSum=(TH1F*)hTPhiMixCorr->Clone();
      hTPhiMixCorrSum->SetName("hTPhiMixCorrSum");
      hTPhiMixCorr2Sum=(TH1F*)hTPhiMixCorr2->Clone();
      hTPhiMixCorr2Sum->SetName("hTPhiMixCorr2Sum");
      hTPhiEtaCorrSum=(TH2F*)hTPhiEtaCorr->Clone();
      hTPhiEtaCorrSum->SetName("hTPhiEtaCorrSum");
      hTPhiEtaCorr2Sum=(TH2F*)hTPhiEtaCorr2->Clone();
      hTPhiEtaCorr2Sum->SetName("hTPhiEtaCorr2Sum");
      hTPhiEtaMixCorrSum=(TH2F*)hTPhiEtaMixCorr->Clone();
      hTPhiEtaMixCorrSum->SetName("hTPhiEtaMixCorrSum");
      hTPhiEtaMixCorr2Sum=(TH2F*)hTPhiEtaMixCorr2->Clone();
      hTPhiEtaMixCorr2Sum->SetName("hTPhiEtaMixCorr2Sum");	
      
      hTPhiCorrPtSum=(TH1F*)hTPhiCorrPt->Clone();
      hTPhiCorrPtSum->SetName("hTPhiCorrPtSum");
      hTPhiCorrPt2Sum=(TH1F*)hTPhiCorrPt2->Clone();
      hTPhiCorrPt2Sum->SetName("hTPhiCorrPt2Sum");
      hTPhiMixCorrPtSum=(TH1F*)hTPhiMixCorrPt->Clone();
      hTPhiMixCorrPtSum->SetName("hTPhiMixCorrPtSum");
      hTPhiMixCorrPt2Sum=(TH1F*)hTPhiMixCorrPt2->Clone();
      hTPhiMixCorrPt2Sum->SetName("hTPhiMixCorrPt2Sum");
      hTPhiEtaCorrPtSum=(TH2F*)hTPhiEtaCorrPt->Clone();
      hTPhiEtaCorrPtSum->SetName("hTPhiEtaCorrPtSum");
      hTPhiEtaCorrPt2Sum=(TH2F*)hTPhiEtaCorrPt2->Clone();
      hTPhiEtaCorrPt2Sum->SetName("hTPhiEtaCorrPt2Sum");
      hTPhiEtaMixCorrPtSum=(TH2F*)hTPhiEtaMixCorrPt->Clone();
      hTPhiEtaMixCorrPtSum->SetName("hTPhiEtaMixCorrPtSum");
      hTPhiEtaMixCorrPt2Sum=(TH2F*)hTPhiEtaMixCorrPt2->Clone();
      hTPhiEtaMixCorrPt2Sum->SetName("hTPhiEtaMixCorrPt2Sum");
      TrigSum=MPt[2];
      hTPhiCorrSum->Scale(TrigSum);
      hTPhiCorr2Sum->Scale(TrigSum);
      hTPhiMixCorrSum->Scale(TrigSum);
      hTPhiMixCorr2Sum->Scale(TrigSum);
      hTPhiEtaCorrSum->Scale(TrigSum);
      hTPhiEtaCorr2Sum->Scale(TrigSum);
      hTPhiEtaMixCorrSum->Scale(TrigSum);
      hTPhiEtaMixCorr2Sum->Scale(TrigSum);
      
      hTPhiCorrPtSum->Scale(TrigSum);
      hTPhiCorrPt2Sum->Scale(TrigSum);
      hTPhiMixCorrPtSum->Scale(TrigSum);
      hTPhiMixCorrPt2Sum->Scale(TrigSum);
      hTPhiEtaCorrPtSum->Scale(TrigSum);
      hTPhiEtaCorrPt2Sum->Scale(TrigSum);
      hTPhiEtaMixCorrPtSum->Scale(TrigSum);
      hTPhiEtaMixCorrPt2Sum->Scale(TrigSum);
    }
    else{     
      hTPhiCorrSum->Add(hTPhiCorr,MPt[2]);
      hTPhiCorr2Sum->Add(hTPhiCorr2,MPt[2]);
      hTPhiMixCorrSum->Add(hTPhiMixCorr,MPt[2]);
      hTPhiMixCorr2Sum->Add(hTPhiMixCorr2,MPt[2]);
      hTPhiEtaCorrSum->Add(hTPhiEtaCorr,MPt[2]);
      hTPhiEtaCorr2Sum->Add(hTPhiEtaCorr2,MPt[2]);
      hTPhiEtaMixCorrSum->Add(hTPhiEtaMixCorr,MPt[2]);
      hTPhiEtaMixCorr2Sum->Add(hTPhiEtaMixCorr2,MPt[2]);
      
      hTPhiCorrPtSum->Add(hTPhiCorr,MPt[2]);
      hTPhiCorrPt2Sum->Add(hTPhiCorr2,MPt[2]);
      hTPhiMixCorrPtSum->Add(hTPhiMixCorr,MPt[2]);
      hTPhiMixCorrPt2Sum->Add(hTPhiMixCorr2,MPt[2]);
      hTPhiEtaCorrPtSum->Add(hTPhiEtaCorr,MPt[2]);
      hTPhiEtaCorrPt2Sum->Add(hTPhiEtaCorr2,MPt[2]);
      hTPhiEtaMixCorrPtSum->Add(hTPhiEtaMixCorr,MPt[2]);
      hTPhiEtaMixCorrPt2Sum->Add(hTPhiEtaMixCorr2,MPt[2]);
      TrigSum+=MPt[2];
    }
  }
  cout << "TrigSum Inclusives: " << TrigSum << endl;
  //rescale to per trigger
  hTPhiCorrSum->Scale(1./TrigSum);
  hTPhiCorr2Sum->Scale(1./TrigSum);
  hTPhiMixCorrSum->Scale(1./TrigSum);
  hTPhiMixCorr2Sum->Scale(1./TrigSum);
  hTPhiEtaCorrSum->Scale(1./TrigSum);
  hTPhiEtaCorr2Sum->Scale(1./TrigSum);
  hTPhiEtaMixCorrSum->Scale(1./TrigSum);
  hTPhiEtaMixCorr2Sum->Scale(1./TrigSum);
  hTPhiCorrPtSum->Scale(1./TrigSum);
  hTPhiCorrPt2Sum->Scale(1./TrigSum);
  hTPhiMixCorrPtSum->Scale(1./TrigSum);
  hTPhiMixCorrPt2Sum->Scale(1./TrigSum);
  hTPhiEtaCorrPtSum->Scale(1./TrigSum);
  hTPhiEtaCorrPt2Sum->Scale(1./TrigSum);
  hTPhiEtaMixCorrPtSum->Scale(1./TrigSum);
  hTPhiEtaMixCorrPt2Sum->Scale(1./TrigSum);

  for(int MainLoop=0;MainLoop<nMainLoop;MainLoop++){
  cout << "nMainLoop " << MainLoop << endl;
    if(APtTPtMult==0){
      APt1=APtBins[MainLoop];
      APt2=APtBins[MainLoop+1];
      MainArray[MainLoop]=(APt1+APt2)/2;
      eMainArray[MainLoop]=(APt2-APt1)/2;
      //cout << "APt1: " << APt1 << " APt2: " << APt2 << endl;
    }
    if(APtTPtMult==1){
      TPt1=TPtBins[MainLoop];
      TPt2=TPtBins[MainLoop+1];
      MainArray[MainLoop]=(TPt1+TPt2)/2;
      eMainArray[MainLoop]=(TPt2-TPt1)/2;
      //cout << "TPt1: " << TPt1 << " TPt2: " << TPt2 << endl;
    }
    if(APtTPtMult==2){
      Mult1=MultBins1[MainLoop];
      Mult2=MultBins2[MainLoop];
      MainArray[MainLoop]=(MultArray1[Mult1]+MultArray2[Mult2]-1)/2;
      eMainArray[MainLoop]=(MultArray2[Mult2]-MultArray1[Mult1]-1)/2;		
    }
    for(int Cent=Mult1;Cent<=Mult2;Cent++){

      MakeProjections(TPt1,TPt2,APt1,APt2,Cent,inList,hPhiRaw,hPhiMixRaw,hEtaNRaw,hEtaNMixRaw,hEtaARaw,hEtaAMixRaw,hPhiEtaRaw,hPhiEtaMixRaw,MPt,0,0,0,LSign);
      MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiEff,hPhiMixEff,hEtaNEff,hEtaNMixEff,hEtaAEff,hEtaAMixEff,hPhiEtaEff,hPhiEtaMixEff,MPt2,0,0,0,LSign);
      MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiMC,hPhiMixMC,hEtaNMC,hEtaNMixMC,hEtaAMC,hEtaAMixMC,hPhiEtaMC,hPhiEtaMixMC,MPt2,0,1,0,LSign);
      
      MakeProjections(TPt1,TPt2,APt1,APt2,Cent,inList,hPhiRawPt,hPhiMixRawPt,hEtaNRawPt,hEtaNMixRawPt,hEtaARawPt,hEtaAMixRawPt,hPhiEtaRawPt,hPhiEtaMixRawPt,MPt,1,0,0,LSign);
      MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiEffPt,hPhiMixEffPt,hEtaNEffPt,hEtaNMixEffPt,hEtaAEffPt,hEtaAMixEffPt,hPhiEtaEffPt,hPhiEtaMixEffPt,MPt2,1,0,0,LSign);
      MakeProjections(TPt1,TPt2,APt1,APt2,Cent,effList,hPhiMCPt,hPhiMixMCPt,hEtaNMCPt,hEtaNMixMCPt,hEtaAMCPt,hEtaAMixMCPt,hPhiEtaMCPt,hPhiEtaMixMCPt,MPt2,1,1,0,LSign);
      /*
      EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hPhiRaw,hPhiEff,hPhiMC,hPhiMixRaw,hPhiMixEff,hPhiMixMC,hEtaNRaw,hEtaNEff,hEtaNMC,hEtaNMixRaw,hEtaNMixEff,hEtaNMixMC,hEtaARaw,hEtaAEff,hEtaAMC,hEtaAMixRaw,hEtaAMixEff,hEtaAMixMC,hPhiEtaRaw,hPhiEtaEff,hPhiEtaMC,hPhiEtaMixRaw,hPhiEtaMixEff,hPhiEtaMixMC,EffMethod);
      EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hPhiRawPt,hPhiEffPt,hPhiMCPt,hPhiMixRawPt,hPhiMixEffPt,hPhiMixMCPt,hEtaNRawPt,hEtaNEffPt,hEtaNMCPt,hEtaNMixRawPt,hEtaNMixEffPt,hEtaNMixMCPt,hEtaARawPt,hEtaAEffPt,hEtaAMCPt,hEtaAMixRawPt,hEtaAMixEffPt,hEtaAMixMCPt,hPhiEtaRawPt,hPhiEtaEffPt,hPhiEtaMCPt,hPhiEtaMixRawPt,hPhiEtaMixEffPt,hPhiEtaMixMCPt,EffMethod);
      */
       if(EffMethod<4){
      EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hPhiEff,hPhiMC,hPhiMixEff,hPhiMixMC,hEtaNEff,hEtaNMC,hEtaNMixEff,hEtaNMixMC,hEtaAEff,hEtaAMC,hEtaAMixEff,hEtaAMixMC,hPhiEtaEff,hPhiEtaMC,hPhiEtaMixEff,hPhiEtaMixMC,EffMethod);
   EffCorr2(TPt1,TPt2,APt1,APt2,Cent,hPhiEffPt,hPhiMCPt,hPhiMixEffPt,hPhiMixMCPt,hEtaNEffPt,hEtaNMCPt,hEtaNMixEffPt,hEtaNMixMCPt,hEtaAEffPt,hEtaAMCPt,hEtaAMixEffPt,hEtaAMixMCPt,hPhiEtaEffPt,hPhiEtaMCPt,hPhiEtaMixEffPt,hPhiEtaMixMCPt,EffMethod);
    }
    else{
      EffFit(APt1,APt2,Cent,effList,hPhiEff,hPhiMixEff,hEtaNEff,hEtaNMixEff,hEtaAEff,hEtaAMixEff,hPhiEtaEff,hPhiEtaMixEff,LSign,VariablePtLimit);
EffFit(APt1,APt2,Cent,effList,hPhiEffPt,hPhiMixEffPt,hEtaNEffPT,hEtaNMixEffPt,hEtaAEffPt,hEtaAMixEffPt,hPhiEtaEffPt,hPhiEtaMixEffPt,LSign,VariablePtLimit);
    }

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
      
      hPhiCorrPt=(TH1F*)hPhiRawPt->Clone();
      sprintf(name,"hPhiCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
      hPhiCorrPt->SetName(name);
      hPhiCorrPt->Divide(hPhiEffPt);
      
      hPhiMixCorrPt=(TH1F*)hPhiMixRawPt->Clone();
      sprintf(name,"hPhiMixCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
      hPhiMixCorrPt->SetName(name);
      hPhiMixCorrPt->Divide(hPhiMixEffPt);

      hEtaNCorrPt=(TH1F*)hEtaNRawPt->Clone();
      sprintf(name,"hEtaNCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
      hEtaNCorrPt->SetName(name);
      hEtaNCorrPt->Divide(hEtaNEffPt);
      
      hEtaNMixCorrPt=(TH1F*)hEtaNMixRawPt->Clone();
      sprintf(name,"hEtaNMixCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
      hEtaNMixCorrPt->SetName(name);
      hEtaNMixCorrPt->Divide(hEtaNMixEffPt);

      hEtaACorrPt=(TH1F*)hEtaARawPt->Clone();
      sprintf(name,"hEtaACorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
      hEtaACorrPt->SetName(name);
      hEtaACorrPt->Divide(hEtaAEffPt);
      
      hEtaAMixCorrPt=(TH1F*)hEtaAMixRawPt->Clone();
      sprintf(name,"hEtaAMixCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
      hEtaAMixCorrPt->SetName(name);
      hEtaAMixCorrPt->Divide(hEtaAMixEffPt);
      
      hPhiEtaCorrPt=(TH2F*)hPhiEtaRawPt->Clone();
      sprintf(name,"hPhiEtaCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
      hPhiEtaCorrPt->SetName(name);
      hPhiEtaCorrPt->Divide(hPhiEtaEffPt);
      
      hPhiEtaMixCorrPt=(TH2F*)hPhiEtaMixRawPt->Clone();
      sprintf(name,"hPhiMixCorrPt_%2.2fPT%2.2f_%2.2fpt%2.2f_%d",TPt1,TPt2,APt1,APt2,Cent);
      hPhiEtaMixCorrPt->SetName(name);
      hPhiEtaMixCorrPt->Divide(hPhiEtaMixEffPt);
   
      //   GeoCorr(hEtaNCorr,hEtaNMixCorr,hEtaNMixMC,hEtaACorr,hEtaAMixCorr,hEtaAMixMC,hPhiEtaCorr,hPhiEtaMixCorr,hPhiEtaMixMC);
      // GeoCorr(hEtaNCorrPt,hEtaNMixCorrPt,hEtaNMixMCPt,hEtaACorrPt,hEtaAMixCorrPt,hEtaAMixMCPt,hPhiEtaCorrPt,hPhiEtaMixCorrPt,hPhiEtaMixMCPt);
      MixedCorrect(hPhiCorr,hPhiMixCorr,hEtaNCorr,hEtaNMixCorr,hEtaACorr,hEtaAMixCorr,hPhiEtaCorr,hPhiEtaMixCorr);
      MixedCorrect(hPhiCorrPt,hPhiMixCorrPt,hEtaNCorrPt,hEtaNMixCorrPt,hEtaACorrPt,hEtaAMixCorrPt,hPhiEtaCorrPt,hPhiEtaMixCorrPt);
      if(ZYAMCent>0)ZYA1(hPhiCorr,hPhiMixCorr,a,ZYAMCent,ZYAMWidth);
      else ZYAM2D(hPhiEtaCorr,hPhiEtaMixCorr,a,72,10);
   

      hPhiCorr2=(TH1F*)hPhiCorr->Clone();
      hPhiCorr2->SetName("hPhiCorr2");
      sprintf(name,"Bg Sub. %s",hPhiCorr2->GetTitle());
      hPhiCorr2->SetTitle(name);
      
      hPhiMixCorr2=(TH1F*)hPhiMixCorr->Clone();
      hPhiMixCorr2->SetName("hPhiMixCorr2");
      sprintf(name,"Normalized %1.3f %s",a[0],hPhiMixCorr2->GetTitle());
      hPhiMixCorr2->SetTitle(name);
      hPhiMixCorr2->SetMarkerStyle(25);
      hPhiMixCorr2->SetMarkerColor(4);
      hPhiMixCorr2->SetLineColor(4);

      hEtaNCorr2=(TH1F*)hEtaNCorr->Clone();
      hEtaNCorr2->SetName("hEtaNCorr2");
      sprintf(name,"Bg Sub. %s",hEtaNCorr2->GetTitle());
      hEtaNCorr2->SetTitle(name);
      
      hEtaNMixCorr2=(TH1F*)hEtaNMixCorr->Clone();
      hEtaNMixCorr2->SetName("hEtaNMixCorr2");
      sprintf(name,"Normalized %1.3f %s",a[0],hEtaNMixCorr2->GetTitle());
      hEtaNMixCorr2->SetTitle(name);
      hEtaNMixCorr2->SetMarkerStyle(25);
      hEtaNMixCorr2->SetMarkerColor(4);
      hEtaNMixCorr2->SetLineColor(4);

      hEtaACorr2=(TH1F*)hEtaACorr->Clone();
      hEtaACorr2->SetName("hEtaACorr2");
      sprintf(name,"Bg Sub. %s",hEtaACorr2->GetTitle());
      hEtaACorr2->SetTitle(name);
      
      hEtaAMixCorr2=(TH1F*)hEtaAMixCorr->Clone();
      hEtaAMixCorr2->SetName("hEtaAMixCorr2");
      sprintf(name,"Normalized %1.3f %s",a[0],hEtaAMixCorr2->GetTitle());
      hEtaAMixCorr2->SetTitle(name);
      hEtaAMixCorr2->SetMarkerStyle(25);
      hEtaAMixCorr2->SetMarkerColor(4);
      hEtaAMixCorr2->SetLineColor(4);
      
      hPhiEtaCorr2=(TH2F*)hPhiEtaCorr->Clone();
      hPhiEtaCorr2->SetName("hPhiEtaCorr2");
      sprintf(name,"BG Sub. %s",hPhiEtaCorr2->GetTitle());
      hPhiEtaCorr2->SetTitle(name);
      
      hPhiEtaMixCorr2=(TH2F*)hPhiEtaMixCorr->Clone();
      hPhiEtaMixCorr2->SetName("hPhiEtaMixCorr2");
      sprintf(name,"Normalized %1.3f %s",a[0],hPhiEtaMixCorr2->GetTitle());
      hPhiEtaMixCorr2->SetTitle(name);
      
      hPhiCorrPt2=(TH1F*)hPhiCorrPt->Clone();
      hPhiCorrPt2->SetName("hPhiCorrPt2");
      sprintf(name,"Bg Sub. %s",hPhiCorrPt2->GetTitle());
      hPhiCorrPt2->SetTitle(name);
      
      hPhiMixCorrPt2=(TH1F*)hPhiMixCorrPt->Clone();
      hPhiMixCorrPt2->SetName("hPhiMixCorrPt2");
      sprintf(name,"Normalized %1.3f %s",a[0],hPhiMixCorrPt2->GetTitle());
      hPhiMixCorrPt2->SetTitle(name);
      hPhiMixCorrPt2->SetMarkerStyle(25);
      hPhiMixCorrPt2->SetMarkerColor(4);
      hPhiMixCorrPt2->SetLineColor(4);

      hEtaNCorrPt2=(TH1F*)hEtaNCorrPt->Clone();
      hEtaNCorrPt2->SetName("hEtaNCorrPt2");
      sprintf(name,"Bg Sub. %s",hEtaNCorrPt2->GetTitle());
      hEtaNCorrPt2->SetTitle(name);
      
      hEtaNMixCorrPt2=(TH1F*)hEtaNMixCorrPt->Clone();
      hEtaNMixCorrPt2->SetName("hEtaNMixCorrPt2");
      sprintf(name,"Normalized %1.3f %s",a[0],hEtaNMixCorrPt2->GetTitle());
      hEtaNMixCorrPt2->SetTitle(name);
      hEtaNMixCorrPt2->SetMarkerStyle(25);
      hEtaNMixCorrPt2->SetMarkerColor(4);
      hEtaNMixCorrPt2->SetLineColor(4);

      hEtaACorrPt2=(TH1F*)hEtaACorrPt->Clone();
      hEtaACorrPt2->SetName("hEtaACorrPt2");
      sprintf(name,"Bg Sub. %s",hEtaACorrPt2->GetTitle());
      hEtaACorrPt2->SetTitle(name);
      
      hEtaAMixCorrPt2=(TH1F*)hEtaAMixCorrPt->Clone();
      hEtaAMixCorrPt2->SetName("hEtaAMixCorrPt2");
      sprintf(name,"Normalized %1.3f %s",a[0],hEtaAMixCorrPt2->GetTitle());
      hEtaAMixCorrPt2->SetTitle(name);
      hEtaAMixCorrPt2->SetMarkerStyle(25);
      hEtaAMixCorrPt2->SetMarkerColor(4);
      hEtaAMixCorrPt2->SetLineColor(4);
      
      hPhiEtaCorrPt2=(TH2F*)hPhiEtaCorrPt->Clone();
      hPhiEtaCorrPt2->SetName("hPhiEtaCorrPt2");
      sprintf(name,"BG Sub. %s",hPhiEtaCorrPt2->GetTitle());
      hPhiEtaCorrPt2->SetTitle(name);
      
      hPhiEtaMixCorrPt2=(TH2F*)hPhiEtaMixCorrPt->Clone();
      hPhiEtaMixCorrPt2->SetName("hPhiEtaMixCorrPt2");
      sprintf(name,"Normalized %1.3f %s",a[0],hPhiEtaMixCorrPt2->GetTitle());
      hPhiEtaMixCorrPt2->SetTitle(name);
      
      hPhiMixCorr2->Scale(a[0]);
      hEtaNMixCorr2->Scale(a[0]);
      hEtaAMixCorr2->Scale(a[0]);
      hPhiEtaMixCorr2->Scale(a[0]);
      hPhiCorr2->Add(hPhiMixCorr2,-1);
      hEtaNCorr2->Add(hEtaNMixCorr2,-1);
      hEtaACorr2->Add(hEtaAMixCorr2,-1);
      Add2D(hPhiEtaCorr2,hPhiEtaMixCorr2,-1); 	
      hPhiMixCorrPt2->Scale(a[0]);
      hPhiEtaMixCorrPt2->Scale(a[0]);
      hEtaNMixCorrPt2->Scale(a[0]);
      hEtaAMixCorrPt2->Scale(a[0]);
      hPhiCorrPt2->Add(hPhiMixCorrPt2,-1);
      hEtaNCorrPt2->Add(hEtaNMixCorrPt2,-1);
      hEtaACorrPt2->Add(hEtaAMixCorrPt2,-1);
      Add2D(hPhiEtaCorrPt2,hPhiEtaMixCorrPt2,-1); 	
      
      //Now we need to add up centrality bins
      if(Cent==Mult1){
	hPhiCorrSum=(TH1F*)hPhiCorr->Clone();
	hPhiCorrSum->SetName("hPhiCorrSum");
	hPhiCorr2Sum=(TH1F*)hPhiCorr2->Clone();
	hPhiCorr2Sum->SetName("hPhiCorr2Sum");
	hPhiMixCorrSum=(TH1F*)hPhiMixCorr->Clone();
	hPhiMixCorrSum->SetName("hPhiMixCorrSum");
	hPhiMixCorr2Sum=(TH1F*)hPhiMixCorr2->Clone();
	hPhiMixCorr2Sum->SetName("hPhiMixCorr2Sum");
	
	hEtaNCorrSum=(TH1F*)hEtaNCorr->Clone();
	hEtaNCorrSum->SetName("hEtaNCorrSum");
	hEtaNCorr2Sum=(TH1F*)hEtaNCorr2->Clone();
	hEtaNCorr2Sum->SetName("hEtaNCorr2Sum");
	hEtaNMixCorrSum=(TH1F*)hEtaNMixCorr->Clone();
	hEtaNMixCorrSum->SetName("hEtaNMixCorrSum");
	hEtaNMixCorr2Sum=(TH1F*)hEtaNMixCorr2->Clone();
	hEtaNMixCorr2Sum->SetName("hEtaNMixCorr2Sum");
	
	hEtaACorrSum=(TH1F*)hEtaACorr->Clone();
	hEtaACorrSum->SetName("hEtaACorrSum");
	hEtaACorr2Sum=(TH1F*)hEtaACorr2->Clone();
	hEtaACorr2Sum->SetName("hEtaACorr2Sum");
	hEtaAMixCorrSum=(TH1F*)hEtaAMixCorr->Clone();
	hEtaAMixCorrSum->SetName("hEtaAMixCorrSum");
	hEtaAMixCorr2Sum=(TH1F*)hEtaAMixCorr2->Clone();
	hEtaAMixCorr2Sum->SetName("hEtaAMixCorr2Sum");

	hPhiEtaCorrSum=(TH2F*)hPhiEtaCorr->Clone();
	hPhiEtaCorrSum->SetName("hPhiEtaCorrSum");
	hPhiEtaCorr2Sum=(TH2F*)hPhiEtaCorr2->Clone();
	hPhiEtaCorr2Sum->SetName("hPhiEtaCorr2Sum");
	hPhiEtaMixCorrSum=(TH2F*)hPhiEtaMixCorr->Clone();
	hPhiEtaMixCorrSum->SetName("hPhiEtaMixCorrSum");
	hPhiEtaMixCorr2Sum=(TH2F*)hPhiEtaMixCorr2->Clone();
	hPhiEtaMixCorr2Sum->SetName("hPhiEtaMixCorr2Sum");	
	
	hPhiCorrPtSum=(TH1F*)hPhiCorrPt->Clone();
	hPhiCorrPtSum->SetName("hPhiCorrPtSum");
	hPhiCorrPt2Sum=(TH1F*)hPhiCorrPt2->Clone();
	hPhiCorrPt2Sum->SetName("hPhiCorrPt2Sum");
	hPhiMixCorrPtSum=(TH1F*)hPhiMixCorrPt->Clone();
	hPhiMixCorrPtSum->SetName("hPhiMixCorrPtSum");
	hPhiMixCorrPt2Sum=(TH1F*)hPhiMixCorrPt2->Clone();
	hPhiMixCorrPt2Sum->SetName("hPhiMixCorrPt2Sum");

	hEtaNCorrPtSum=(TH1F*)hEtaNCorrPt->Clone();
	hEtaNCorrPtSum->SetName("hEtaNCorrPtSum");
	hEtaNCorrPt2Sum=(TH1F*)hEtaNCorrPt2->Clone();
	hEtaNCorrPt2Sum->SetName("hEtaNCorrPt2Sum");
	hEtaNMixCorrPtSum=(TH1F*)hEtaNMixCorrPt->Clone();
	hEtaNMixCorrPtSum->SetName("hEtaNMixCorrPtSum");
	hEtaNMixCorrPt2Sum=(TH1F*)hEtaNMixCorrPt2->Clone();
	hEtaNMixCorrPt2Sum->SetName("hEtaNMixCorrPt2Sum");

	hEtaACorrPtSum=(TH1F*)hEtaACorrPt->Clone();
	hEtaACorrPtSum->SetName("hEtaACorrPtSum");
	hEtaACorrPt2Sum=(TH1F*)hEtaACorrPt2->Clone();
	hEtaACorrPt2Sum->SetName("hEtaACorrPt2Sum");
	hEtaAMixCorrPtSum=(TH1F*)hEtaAMixCorrPt->Clone();
	hEtaAMixCorrPtSum->SetName("hEtaAMixCorrPtSum");
	hEtaAMixCorrPt2Sum=(TH1F*)hEtaAMixCorrPt2->Clone();
	hEtaAMixCorrPt2Sum->SetName("hEtaAMixCorrPt2Sum");

	hPhiEtaCorrPtSum=(TH2F*)hPhiEtaCorrPt->Clone();
	hPhiEtaCorrPtSum->SetName("hPhiEtaCorrPtSum");
	hPhiEtaCorrPt2Sum=(TH2F*)hPhiEtaCorrPt2->Clone();
	hPhiEtaCorrPt2Sum->SetName("hPhiEtaCorrPt2Sum");
	hPhiEtaMixCorrPtSum=(TH2F*)hPhiEtaMixCorrPt->Clone();
	hPhiEtaMixCorrPtSum->SetName("hPhiEtaMixCorrPtSum");
	hPhiEtaMixCorrPt2Sum=(TH2F*)hPhiEtaMixCorrPt2->Clone();
	hPhiEtaMixCorrPt2Sum->SetName("hPhiEtaMixCorrPt2Sum");	
	
	TrigSum=MPt[2];
  	
	hPhiCorrSum->Scale(TrigSum);
	hPhiCorr2Sum->Scale(TrigSum);
	hPhiMixCorrSum->Scale(TrigSum);
	hPhiMixCorr2Sum->Scale(TrigSum);
	hEtaNCorrSum->Scale(TrigSum);
	hEtaNCorr2Sum->Scale(TrigSum);
	hEtaNMixCorrSum->Scale(TrigSum);
	hEtaNMixCorr2Sum->Scale(TrigSum);
	hEtaACorrSum->Scale(TrigSum);
	hEtaACorr2Sum->Scale(TrigSum);
	hEtaAMixCorrSum->Scale(TrigSum);
	hEtaAMixCorr2Sum->Scale(TrigSum);
	hPhiEtaCorrSum->Scale(TrigSum);
	hPhiEtaCorr2Sum->Scale(TrigSum);
	hPhiEtaMixCorrSum->Scale(TrigSum);
	hPhiEtaMixCorr2Sum->Scale(TrigSum);
  	
	hPhiCorrPtSum->Scale(TrigSum);
	hPhiCorrPt2Sum->Scale(TrigSum);
	hPhiMixCorrPtSum->Scale(TrigSum);
	hPhiMixCorrPt2Sum->Scale(TrigSum);
	hEtaNCorrPtSum->Scale(TrigSum);
	hEtaNCorrPt2Sum->Scale(TrigSum);
	hEtaNMixCorrPtSum->Scale(TrigSum);
	hEtaNMixCorrPt2Sum->Scale(TrigSum);
	hEtaACorrPtSum->Scale(TrigSum);
	hEtaACorrPt2Sum->Scale(TrigSum);
	hEtaAMixCorrSum->Scale(TrigSum);
	hEtaAMixCorr2Sum->Scale(TrigSum);
	hPhiEtaCorrPtSum->Scale(TrigSum);
	hPhiEtaCorrPt2Sum->Scale(TrigSum);
	hPhiEtaMixCorrPtSum->Scale(TrigSum);
	hPhiEtaMixCorrPt2Sum->Scale(TrigSum);
      }
      else{//I need to do trigger weighted averaging here so I need the trigger counts
	hPhiCorrSum->Add(hPhiCorr,MPt[2]);
	hPhiCorr2Sum->Add(hPhiCorr2,MPt[2]);
	hPhiMixCorrSum->Add(hPhiMixCorr,MPt[2]);
	hPhiMixCorr2Sum->Add(hPhiMixCorr2,MPt[2]);
	hEtaNCorrSum->Add(hEtaNCorr,MPt[2]);
	hEtaNCorr2Sum->Add(hEtaNCorr2,MPt[2]);
	hEtaNMixCorrSum->Add(hEtaNMixCorr,MPt[2]);
	hEtaNMixCorr2Sum->Add(hEtaNMixCorr2,MPt[2]);
	hEtaACorrSum->Add(hEtaACorr,MPt[2]);
	hEtaACorr2Sum->Add(hEtaNCorr2,MPt[2]);
	hEtaAMixCorrSum->Add(hEtaAMixCorr,MPt[2]);
	hEtaAMixCorr2Sum->Add(hEtaAMixCorr2,MPt[2]);
	hPhiEtaCorrSum->Add(hPhiEtaCorr,MPt[2]);
	hPhiEtaCorr2Sum->Add(hPhiEtaCorr2,MPt[2]);
	hPhiEtaMixCorrSum->Add(hPhiEtaMixCorr,MPt[2]);
	hPhiEtaMixCorr2Sum->Add(hPhiEtaMixCorr2,MPt[2]);
  	
	hPhiCorrPtSum->Add(hPhiCorr,MPt[2]);
	hPhiCorrPt2Sum->Add(hPhiCorr2,MPt[2]);
	hPhiMixCorrPtSum->Add(hPhiMixCorr,MPt[2]);
	hPhiMixCorrPt2Sum->Add(hPhiMixCorr2,MPt[2]);
	hEtaNCorrPtSum->Add(hEtaNCorrPt,MPt[2]);
	hEtaNCorrPt2Sum->Add(hEtaNCorrPt2,MPt[2]);
	hEtaNMixCorrPtSum->Add(hEtaNMixCorrPt,MPt[2]);
	hEtaNMixCorrPt2Sum->Add(hEtaNMixCorrPt2,MPt[2]);
	hEtaACorrPtSum->Add(hEtaACorrPt,MPt[2]);
	hEtaACorrPt2Sum->Add(hEtaNCorrPt2,MPt[2]);
	hEtaAMixCorrPtSum->Add(hEtaAMixCorrPt,MPt[2]);
	hEtaAMixCorrPt2Sum->Add(hEtaAMixCorrPt2,MPt[2]);
	hPhiEtaCorrPtSum->Add(hPhiEtaCorr,MPt[2]);
	hPhiEtaCorrPt2Sum->Add(hPhiEtaCorr2,MPt[2]);
	hPhiEtaMixCorrPtSum->Add(hPhiEtaMixCorr,MPt[2]);
	hPhiEtaMixCorrPt2Sum->Add(hPhiEtaMixCorr2,MPt[2]);
	TrigSum+=MPt[2];
      }
    }
    //rescale to per trigger
    // cout << "TrigSum " << TrigSum << endl;
    hPhiCorrSum->Scale(1./TrigSum);
    hPhiCorr2Sum->Scale(1./TrigSum);
    hPhiMixCorrSum->Scale(1./TrigSum);
    hPhiMixCorr2Sum->Scale(1./TrigSum);
    hEtaNCorrSum->Scale(1./TrigSum);
    hEtaNCorr2Sum->Scale(1./TrigSum);
    hEtaNMixCorrSum->Scale(1./TrigSum);
    hEtaNMixCorr2Sum->Scale(1./TrigSum);
    hEtaACorrSum->Scale(1./TrigSum);
    hEtaACorr2Sum->Scale(1./TrigSum);
    hEtaAMixCorrSum->Scale(1./TrigSum);
    hEtaAMixCorr2Sum->Scale(1./TrigSum);
    hPhiEtaCorrSum->Scale(1./TrigSum);
    hPhiEtaCorr2Sum->Scale(1./TrigSum);
    hPhiEtaMixCorrSum->Scale(1./TrigSum);
    hPhiEtaMixCorr2Sum->Scale(1./TrigSum);
    
    hPhiCorrPtSum->Scale(1./TrigSum);
    hPhiCorrPt2Sum->Scale(1./TrigSum);
    hPhiMixCorrPtSum->Scale(1./TrigSum);
    hPhiMixCorrPt2Sum->Scale(1./TrigSum);
    hEtaNCorrPtSum->Scale(1./TrigSum);
    hEtaNCorrPt2Sum->Scale(1./TrigSum);
    hEtaNMixCorrPtSum->Scale(1./TrigSum);
    hEtaNMixCorrPt2Sum->Scale(1./TrigSum);
    hEtaACorrPtSum->Scale(1./TrigSum);
    hEtaACorrPt2Sum->Scale(1./TrigSum);
    hEtaAMixCorrPtSum->Scale(1./TrigSum);
    hEtaAMixCorrPt2Sum->Scale(1./TrigSum);
    hPhiEtaCorrPtSum->Scale(1./TrigSum);
    hPhiEtaCorrPt2Sum->Scale(1./TrigSum);
    hPhiEtaMixCorrPtSum->Scale(1./TrigSum);
    hPhiEtaMixCorrPt2Sum->Scale(1./TrigSum);
    
    //  for(int xs=0;xs<nx;xs++){xarray[xs][MainLoop]=0;}
    for(int yields=0;yields<nyields;yields++){  
      NYieldPhi[yields][MainLoop]=0;
      NYieldPhiZYAM[yields][MainLoop]=0;
      NYieldPhiFit[yields][MainLoop]=0;
      NYieldEta[yields][MainLoop]=0;
      NYieldEtaZYAM[yields][MainLoop]=0;
      NYieldEtaFit[yields][MainLoop]=0;
      NYieldPhiEta1[yields][MainLoop]=0;
      NYieldPhiEta2[yields][MainLoop]=0;
      NYieldPhiEta1ZYAM[yields][MainLoop]=0;
      NYieldPhiEta2ZYAM[yields][MainLoop]=0;
      NRMSPhi[yields][MainLoop]=0;
      NWidthPhi[yields][MainLoop]=0;
      AYieldPhi[yields][MainLoop]=0;
      AYieldPhiZYAM[yields][MainLoop]=0;
      AYieldPhiFit[yields][MainLoop]=0;
      AYieldEta[yields][MainLoop]=0;
      AYieldEtaZYAM[yields][MainLoop]=0;
      AYieldEtaFit[yields][MainLoop]=0;
      AYieldPhiEta[yields][MainLoop]=0;
      AYieldPhiEtaZYAM[yields][MainLoop]=0;
      MYieldPhi[yields][MainLoop]=0;
      MYieldPhi[yields][MainLoop]=0;
      MAvePhiEta[yields][MainLoop]=0;
      NAvePhi[yields][MainLoop]=0;
      NAvePhiZYAM[yields][MainLoop]=0;
      NAvePhiFit[yields][MainLoop]=0;
      NAveEta[yields][MainLoop]=0;
      NAveEtaZYAM[yields][MainLoop]=0;
      NAveEtaFit[yields][MainLoop]=0;
      NAvePhiEta1[yields][MainLoop]=0;
      NAvePhiEta2[yields][MainLoop]=0;
      NAvePhiEta1ZYAM[yields][MainLoop]=0;
      NAvePhiEta2ZYAM[yields][MainLoop]=0;
      AAvePhi[yields][MainLoop]=0;
      AAvePhiZYAM[yields][MainLoop]=0;
      AAvePhiFit[yields][MainLoop]=0;
      AAveEta[yields][MainLoop]=0;
      AAveEtaZYAM[yields][MainLoop]=0;
      AAveEtaFit[yields][MainLoop]=0;
      AAvePhiEta[yields][MainLoop]=0;
      AAvePhiEtaZYAM[yields][MainLoop]=0;
      MAvePhi[yields][MainLoop]=0;
      MAvePhiFit[yields][MainLoop]=0;
      MAvePhiEta[yields][MainLoop]=0;
      
      TNYieldPhi[yields][MainLoop]=0;
      TNYieldPhiZYAM[yields][MainLoop]=0;
      TNYieldPhiFit[yields][MainLoop]=0;
      TNYieldPhiEta1[yields][MainLoop]=0;
      TNYieldPhiEta2[yields][MainLoop]=0;
      TNYieldPhiEta3[yields][MainLoop]=0;
      TAYieldPhi[yields][MainLoop]=0;
      TAYieldPhiZYAM[yields][MainLoop]=0;
      TAYieldPhiFit[yields][MainLoop]=0;
      TAYieldPhiEta[yields][MainLoop]=0;
      TAYieldPhiEtaZYAM[yields][MainLoop]=0;
    }
    //lets do the dphi histo yeilds
    //Get our yields for the $\Delta\phi$ Histograms (from bin counting)
    float nbinsPhiEtaPhi=hPhiEtaCorr->GetNbinsX();
    float nbinsPhiEtaEta=hPhiEtaCorr->GetNbinsY();
    binwx=hPhiCorr->GetBinWidth(1);
    NearBinsPhi=0;
    AwayBinsPhi=0;
    MinBinsPhi=0;
    NearBinsPhiEta1=0;
    NearBinsPhiEta2=0;
    AwayBinsPhiEta=0;
    MinBinsPhiEta=0;
    
    //Bins with 0 content have messed up error bars and screw up the fits
    for(int xx=1;xx<=hPhiCorrSum->GetNbinsX();xx++){
      if(hPhiCorrSum->GetBinContent(xx)==0)hPhiCorrSum->SetBinError(xx,1./hPhiCorrSum->GetBinWidth(xx)/TrigSum);
      if(hPhiCorrPtSum->GetBinContent(xx)==0)hPhiCorrSum->SetBinError(xx,1./pow(hPhiCorrPtSum->GetBinWidth(xx),2)/TrigSum);
    }
 for(int xx=1;xx<=hEtaNCorrSum->GetNbinsX();xx++){
   if(hEtaNCorrSum->GetBinContent(xx)==0){
     float acp=1-fabs(hEtaNCorr->GetBinCenter(xx))/1.8;
     hEtaNCorrSum->SetBinError(xx,1./hEtaNCorrSum->GetBinWidth(xx)/TrigSum/acp);
     // cout << "Error Eta" << hEtaNCorrSum->GetBinContent(xx) << " " << hEtaNCorrSum->GetBinError(xx) << endl;
   }
   if(hEtaNCorrPtSum->GetBinContent(xx)==0)hEtaNCorrPtSum->SetBinError(xx,1./hEtaNCorrPtSum->GetBinWidth(xx)/TrigSum/acp);
    }
    
    //0 yield, 1 eyield, 2 pt yield, 3 ept yield, 4 pt/pt_trig yield, 5 ezt yield, 6 pt/pt_near yield, 7 ezt2 yield
    trms1=0;trms2=0;trms3=0;
    char outName[200];
    cfit1=new TCanvas("cfit1","",800,600);
    //Fit The Histograms
    fit1->SetParameters(1,0.2,1,0.4,0.5);
    hPhiCorrSum->Fit("fit1");
    NYieldPhiFit[0][MainLoop]=fit1->GetParameter(0);
    NYieldPhiFit[1][MainLoop]=fit1->GetParError(0);
    NWidthPhi[0][MainLoop]=fit1->GetParameter(1);
    NWidthPhi[1][MainLoop]=fit1->GetParError(1);
    AYieldPhiFit[0][MainLoop]=fit1->GetParameter(2);
    AYieldPhiFit[1][MainLoop]=fit1->GetParError(2);
    AWidthPhi[0][MainLoop]=fit1->GetParameter(3);
    AWidthPhi[1][MainLoop]=fit1->GetParError(3);
    MYieldPhiFit[0][MainLoop]=fit1->GetParameter(4)*2*Pi;
    MYieldPhiFit[1][MainLoop]=fit1->GetParError(4)*2*Pi;
    sprintf(outName,"%s/DrawSpectra_FitPhiCorr%d_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_%d_%d%s",Folder,EffMethod,TPt1,TPt2,APt1,APt2,Mult1,Mult2,APtTPtMult,MainLoop,filetype);
    if(SaveFits)cfit1->SaveAs(outName);

    cfit2=new TCanvas("cfit2","",800,600);
    //float par4=fit1->GetParameter(4)*1.5*4/(hEtaNCorrSum->GetNbinsX()*hEtaNCorrSum->GetBinWidth(1));//proper scaling (selfcorrects for different bin etc)
    //  float par4=fit1->GetParameter(4)*hPhiCorrSum->GetBinWidth(1)/hEtaNCorrSum->GetBinWidth(1)*0.8/3.14159;
    //float par4=fit1->GetParameter(4)*3/3.2;
 float par4=fit1->GetParameter(4)*6.4/3;
    fit2->SetParameter(2,par4);
    // fit2->SetParLimits(2,par4,par4);
    fit2->SetParameter(0,fit1->GetParameter(0));
    fit2->SetParameter(1,fit1->GetParameter(1));
    hEtaNCorrSum->Fit("fit2");
    NYieldEtaFit[0][MainLoop]=fit2->GetParameter(0);
    NYieldEtaFit[1][MainLoop]=fit2->GetParError(0);
    NWidthEta[0][MainLoop]=fit2->GetParameter(1);
    NWidthEta[1][MainLoop]=fit2->GetParError(1);
    sprintf(outName,"%s/DrawSpectra_FitEtaCorr%d_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_%d_%d%s",Folder,EffMethod,TPt1,TPt2,APt1,APt2,Mult1,Mult2,APtTPtMult,MainLoop,filetype);
    if(SaveFits)cfit2->SaveAs(outName);

    hPhiCorrPtSum->Fit("fit1");
    NYieldPhiFit[2][MainLoop]=fit1->GetParameter(0);
    NYieldPhiFit[3][MainLoop]=fit1->GetParError(0);
    NWidthPhi[2][MainLoop]=fit1->GetParameter(1);
    NWidthPhi[3][MainLoop]=fit1->GetParError(1);
    AYieldPhiFit[2][MainLoop]=fit1->GetParameter(2);
    AYieldPhiFit[3][MainLoop]=fit1->GetParError(2);
    AWidthPhi[2][MainLoop]=fit1->GetParameter(3);
    AWidthPhi[3][MainLoop]=fit1->GetParError(3);
    MYieldPhiFit[2][MainLoop]=fit1->GetParameter(4)*2*Pi;
    MYieldPhiFit[3][MainLoop]=fit1->GetParError(4)*2*Pi;

  
    // par4=fit1->GetParameter(4)*hPhiCorrSum->GetBinWidth(1)/hEtaNCorrSum->GetBinWidth(1)*1.5/3.1415962;
    //float par4=fit1->GetParameter(4)*1.5*2/1.8;
    //float par4=fit1->GetParameter(4)*1.5*4/(hEtaNCorrPtSum->GetNbinsX()*hEtaNCorrPtSum->GetBinWidth(1));
    //float par4=fit1->GetParameter(4)*hPhiCorrSum->GetBinWidth(1)/hEtaNCorrSum->GetBinWidth(1)*1.5/3.14159;
    float par4=fit1->GetParameter(4)*6.4/3;
    fit2->SetParameter(2,par4);
    //fit2->SetParLimits(2,par4,par4);
    fit2->SetParameter(0,fit1->GetParameter(0));
    fit2->SetParameter(1,fit1->GetParameter(1));
    //hEtaNCorrSum->Fit("fit2");
    // fit2->SetParLimits(2,fit1->GetParameter(4),fit1->GetParameter(4));
    hEtaNCorrPtSum->Fit("fit2");
    NYieldEtaFit[2][MainLoop]=fit2->GetParameter(0);
    NYieldEtaFit[3][MainLoop]=fit2->GetParError(0);
    NWidthEta[2][MainLoop]=fit2->GetParameter(1);
    NWidthEta[3][MainLoop]=fit2->GetParError(1);

    //Not doing any away-side fits now but setting up non-zero values so no divide by zero later as code is inplace in case there are fits
    AYieldEtaFit[0][MainLoop]=1;
    AYieldEtaFit[1][MainLoop]=1;
    AWidthEta[0][MainLoop]=1;
    AWidthEta[1][MainLoop]=1;
    AYieldEtaFit[2][MainLoop]=1;
    AYieldEtaFit[3][MainLoop]=1;
    AWidthEta[2][MainLoop]=1;
    AWidthEta[3][MainLoop]=1;

    //Average pT from fits
    NAvePhiFit[0][MainLoop]=NYieldPhiFit[2][MainLoop]/NYieldPhiFit[0][MainLoop];
    NAvePhiFit[1][MainLoop]=NAvePhiFit[0][MainLoop]*pow(pow(NYieldPhiFit[3][MainLoop]/NYieldPhiFit[2][MainLoop],2)+pow(NYieldPhiFit[1][MainLoop]/NYieldPhiFit[0][MainLoop],2),0.5);
    AAvePhiFit[0][MainLoop]=AYieldPhiFit[2][MainLoop]/AYieldPhiFit[0][MainLoop];
    AAvePhiFit[1][MainLoop]=AAvePhiFit[0][MainLoop]*pow(pow(AYieldPhiFit[3][MainLoop]/AYieldPhiFit[2][MainLoop],2)+pow(AYieldPhiFit[1][MainLoop]/AYieldPhiFit[0][MainLoop],2),0.5);
    MAvePhiFit[0][MainLoop]=MYieldPhiFit[2][MainLoop]/MYieldPhiFit[0][MainLoop];
    MAvePhiFit[1][MainLoop]=MAvePhiFit[0][MainLoop]*pow(pow(MYieldPhiFit[3][MainLoop]/MYieldPhiFit[2][MainLoop],2)+pow(MYieldPhiFit[1][MainLoop]/MYieldPhiFit[0][MainLoop],2),0.5);
    NAveEtaFit[0][MainLoop]=NYieldEtaFit[2][MainLoop]/NYieldEtaFit[0][MainLoop];
    NAveEtaFit[1][MainLoop]=NAveEtaFit[0][MainLoop]*pow(pow(NYieldEtaFit[3][MainLoop]/NYieldEtaFit[2][MainLoop],2)+pow(NYieldEtaFit[1][MainLoop]/NYieldEtaFit[0][MainLoop],2),0.5);

    //Computer Jt from the fits
    float meanptA=NYieldPhiFit[2][MainLoop]/NYieldPhiFit[0][MainLoop];
    float EmeanptA=meanptA*sqrt(pow(NYieldPhiFit[3][MainLoop]/NYieldPhiFit[2][MainLoop],2)+pow(NYieldPhiFit[1][MainLoop]/NYieldPhiFit[0][MainLoop],2));
    NJtPhi[0][MainLoop]=NWidthPhi[0][MainLoop]*meanptA*MPt[0]*pow(2/(meanptA*meanptA+MPt[0]*MPt[0]),0.5);
    float temp1=pow(MPt[0]*MPt[0]+meanptA*meanptA,0.5);
    float temp2=2*EmeanptA*meanptA+2*MPt[0]*MPt[1];
    float temp3=0.5*temp2/temp1;
    NJtPhi[1][MainLoop]=NJtPhi[0][MainLoop]*sqrt(pow(NWidthPhi[1][MainLoop]/NWidthPhi[0][MainLoop],2)+pow(EmeanptA/meanptA,2)+pow(MPt[1]/MPt[0],2)+temp3*temp3);

    meanptA=NYieldEtaFit[2][MainLoop]/NYieldEtaFit[0][MainLoop];
    EmeanptA=meanptA*sqrt(pow(NYieldEtaFit[3][MainLoop]/NYieldEtaFit[2][MainLoop],2)+pow(NYieldEtaFit[1][MainLoop]/NYieldEtaFit[0][MainLoop],2));
    NJtEta[0][MainLoop]=NWidthEta[0][MainLoop]*meanptA*MPt[0]*pow(2/(meanptA*meanptA+MPt[0]*MPt[0]),0.5);
    temp1=pow(MPt[0]*MPt[0]+meanptA*meanptA,0.5);
    temp2=2*EmeanptA*meanptA+2*MPt[0]*MPt[1];
    temp3=0.5*temp2/temp1;
    NJtEta[1][MainLoop]=NJtEta[0][MainLoop]*sqrt(pow(NWidthEta[1][MainLoop]/NWidthEta[0][MainLoop],2)+pow(EmeanptA/meanptA,2)+pow(MPt[1]/MPt[0],2)+temp3*temp3);

    //PhiLoop
    for(int x1=0;x1<hPhiCorrSum->GetNbinsX();x1++){
      centx=hPhiCorrSum->GetBinCenter(x1);
      //Near-Side 1d
      if(fabs(centx)<NearWidthPhi||fabs(centx-2*Pi)<NearWidthPhi||fabs(centx+2*Pi)<NearWidthPhi){
	NYieldPhi[0][MainLoop]+=binwx*hPhiCorrSum->GetBinContent(x1);
	NYieldPhi[1][MainLoop]+=pow(binwx*hPhiCorrSum->GetBinError(x1),2);
	NYieldPhi[2][MainLoop]+=binwx*hPhiCorrPtSum->GetBinContent(x1);
	NYieldPhi[3][MainLoop]+=pow(binwx*hPhiCorrPtSum->GetBinError(x1),2);
	NYieldPhiZYAM[0][MainLoop]+=binwx*hPhiCorr2Sum->GetBinContent(x1);
	NYieldPhiZYAM[1][MainLoop]+=pow(binwx*hPhiCorr2Sum->GetBinError(x1),2);
	NYieldPhiZYAM[2][MainLoop]+=binwx*hPhiCorrPt2Sum->GetBinContent(x1);
	NYieldPhiZYAM[3][MainLoop]+=pow(binwx*hPhiCorrPt2Sum->GetBinError(x1),2);
	if(APtTPtMult==0){
	  TNYieldPhi[0][MainLoop]+=binwx*hTPhiCorrSum->GetBinContent(x1);
	  TNYieldPhi[1][MainLoop]+=pow(binwx*hTPhiCorrSum->GetBinError(x1),2);
	  TNYieldPhi[2][MainLoop]+=binwx*hTPhiCorrPtSum->GetBinContent(x1);
	  TNYieldPhi[3][MainLoop]+=pow(binwx*hTPhiCorrPtSum->GetBinError(x1),2);
	  TNYieldPhiZYAM[0][MainLoop]+=binwx*hTPhiCorr2Sum->GetBinContent(x1);
	  TNYieldPhiZYAM[1][MainLoop]+=pow(binwx*hTPhiCorr2Sum->GetBinError(x1),2);
	  TNYieldPhiZYAM[2][MainLoop]+=binwx*hTPhiCorrPt2Sum->GetBinContent(x1);
	  TNYieldPhiZYAM[3][MainLoop]+=pow(binwx*hTPhiCorrPt2Sum->GetBinError(x1),2);
	}
	NRMSPhi[0][MainLoop]+=pow(binwx*centx*hPhiCorr2Sum->GetBinContent(x1),2);
	NRMSPhi[1][MainLoop]+=pow(2*binwx*centx*hPhiCorr2Sum->GetBinError(x1),2);
	NRMSPhi[2][MainLoop]+=pow(binwx*centx*hPhiCorrPt2Sum->GetBinContent(x1),2);
	NRMSPhi[3][MainLoop]+=pow(2*binwx*centx*hPhiCorrPt2Sum->GetBinError(x1),2);
	//how do I want to do the errors
	//1/RMS*(2*sum cont*centx/sum cont)(1/sum cont^2)*err
	//later do NRMSPhi[0][MainLoop]=sqrt(NRMSPhi[0][Mainloop]/NYieldPhi[0][MainLoop]) and NRMSPhi[1][MainLoop]=0.5*sqrt(NRMSPhi[1]/NRMSPhi[0]^2+NYieldPhi[1]/NYieldPhi[0]^2) I'm assuming no sqrts on errors yet
	NearBinsPhi++;
      }
      //AwaySide 1d
      if(fabs(centx-Pi)<AwayWidthPhi||fabs(centx-3*Pi)<AwayWidthPhi||fabs(centx+Pi)<AwayWidthPhi){
	trms1=(centx-Pi);
	if(trms1<Pi)trms1+=2*Pi;
	if(trms1>Pi)trms1-=2*Pi;
	AYieldPhi[0][MainLoop]+=binwx*hPhiCorrSum->GetBinContent(x1);
	AYieldPhi[1][MainLoop]+=pow(binwx*hPhiCorrSum->GetBinContent(x1),2);
	AYieldPhi[2][MainLoop]+=binwx*hPhiCorrPtSum->GetBinContent(x1);
	AYieldPhi[3][MainLoop]+=pow(binwx*hPhiCorrPtSum->GetBinError(x1),2);
	AYieldPhiZYAM[0][MainLoop]+=binwx*hPhiCorr2Sum->GetBinContent(x1);
	AYieldPhiZYAM[1][MainLoop]+=pow(binwx*hPhiCorr2Sum->GetBinError(x1),2);
	AYieldPhiZYAM[2][MainLoop]+=binwx*hPhiCorrPt2Sum->GetBinContent(x1);
	AYieldPhiZYAM[3][MainLoop]+=pow(binwx*hPhiCorrPt2Sum->GetBinError(x1),2);
	ARMSPhi[0][MainLoop]+=pow(binwx*trms1*hPhiCorr2Sum->GetBinContent(x1),2);
	ARMSPhi[1][MainLoop]+=pow(2*binwx*trms1*hPhiCorr2Sum->GetBinError(x1),2);
	ARMSPhi[2][MainLoop]+=pow(binwx*trms1*hPhiCorrPt2Sum->GetBinContent(x1),2);
	ARMSPhi[3][MainLoop]+=pow(2*binwx*trms1*hPhiCorrPt2Sum->GetBinError(x1),2);
	AwayBinsPhi++;
      }
      //Minimum
      if(fabs(centx-ZYAMCent)<ZYAMWidth||fabs(centx+ZYAMCent)<ZYAMWidth||fabs(centx-ZYAMCent+2*Pi)<ZYAMWidth||fabs(centx-ZYAMCent-2*Pi)<ZYAMWidth||fabs(centx+ZYAMCent+2*Pi)<ZYAMWidth||fabs(centx+ZYAMCent-2*Pi)<ZYAMWidth){
	MYieldPhi[0][MainLoop]+=hPhiCorrSum->GetBinContent(x1);
	MYieldPhi[1][MainLoop]+=pow(hPhiCorrSum->GetBinError(x1),2);
	MYieldPhi[2][MainLoop]+=hPhiCorrPtSum->GetBinContent(x1);
	MYieldPhi[3][MainLoop]+=pow(hPhiCorrPtSum->GetBinError(x1),2);
	MinBinsPhi++;
      }
    }//Phi Loop

    //EtaLoop
    binwx=hEtaNCorrSum->GetBinWidth(1);
    for(int x1=1;x1<=hEtaNCorrSum->GetNbinsX();x1++){
      NYieldEta[0][MainLoop]+=binwx*hEtaNCorrSum->GetBinContent(x1);
      NYieldEta[1][MainLoop]+=pow(binwx*hEtaNCorrSum->GetBinError(x1),2);
      NYieldEta[2][MainLoop]+=binwx*hEtaNCorrPtSum->GetBinContent(x1);
      NYieldEta[3][MainLoop]+=pow(binwx*hEtaNCorrPtSum->GetBinError(x1),2);
      NYieldEtaZYAM[0][MainLoop]+=binwx*hEtaNCorr2Sum->GetBinContent(x1);
      NYieldEtaZYAM[1][MainLoop]+=pow(binwx*hEtaNCorr2Sum->GetBinError(x1),2);
      NYieldEtaZYAM[2][MainLoop]+=binwx*hEtaNCorrPt2Sum->GetBinContent(x1);
      NYieldEtaZYAM[3][MainLoop]+=pow(binwx*hEtaNCorrPt2Sum->GetBinError(x1),2);
      AYieldEta[0][MainLoop]+=binwx*hEtaACorrSum->GetBinContent(x1);
      AYieldEta[1][MainLoop]+=pow(binwx*hEtaACorrSum->GetBinError(x1),2);
      AYieldEta[2][MainLoop]+=binwx*hEtaACorrPtSum->GetBinContent(x1);
      AYieldEta[3][MainLoop]+=pow(binwx*hEtaACorrPtSum->GetBinError(x1),2);
      AYieldEtaZYAM[0][MainLoop]+=binwx*hEtaACorr2Sum->GetBinContent(x1);
      AYieldEtaZYAM[1][MainLoop]+=pow(binwx*hEtaACorr2Sum->GetBinError(x1),2);
      AYieldEtaZYAM[2][MainLoop]+=binwx*hEtaACorrPt2Sum->GetBinContent(x1);
      AYieldEtaZYAM[3][MainLoop]+=pow(binwx*hEtaACorrPt2Sum->GetBinError(x1),2);
    }//EtaLoop
      

    phiaxis=hPhiEtaCorrSum->GetXaxis();
    etaaxis=hPhiEtaCorrSum->GetYaxis();
    binwx=phiaxis->GetBinWidth(1);
    binwy=etaaxis->GetBinWidth(1);
    for(int x1=0;x1<hPhiEtaCorrSum->GetNbinsX();x1++){
      centx=phiaxis->GetBinCenter(x1);
      for(int y1=0;y1<hPhiEtaCorrSum->GetNbinsY();y1++){
	centy=etaaxis->GetBinCenter(y1);
	//Near-Side 2 d
	if(fabs(centx)<NearWidthPhi||fabs(centx-2*Pi)<NearWidthPhi||fabs(centx+2*Pi)<NearWidthPhi){
	  //Jet Peak Area
	  if(fabs(centy)<NearWidthEta){
	    NYieldPhiEta1[0][MainLoop]+=binwx*binwy*hPhiEtaCorrSum->GetBinContent(x1,y1);
	    NYieldPhiEta1[1][MainLoop]+=pow(binwx*binwy*hPhiEtaCorrSum->GetBinError(x1,y1),2);
	    NYieldPhiEta1[2][MainLoop]+=binwx*binwy*hPhiEtaCorrPtSum->GetBinContent(x1,y1);
	    NYieldPhiEta1[3][MainLoop]+=pow(binwx*binwy*hPhiEtaCorrPtSum->GetBinError(x1,y1),2);
	    NYieldPhiEta1ZYAM[0][MainLoop]+=binwx*binwy*hPhiEtaCorr2Sum->GetBinContent(x1,y1);
	    NYieldPhiEta1ZYAM[1][MainLoop]+=pow(binwx*binwy*hPhiEtaCorr2Sum->GetBinError(x1,y1),2);
	    NYieldPhiEta1ZYAM[2][MainLoop]+=binwx*binwy*hPhiEtaCorrPt2Sum->GetBinContent(x1,y1);
	    NYieldPhiEta1ZYAM[3][MainLoop]+=pow(binwx*binwy*hPhiEtaCorrPt2Sum->GetBinError(x1,y1),2);
	    NearBinsPhiEta1++;
	  }
	  //Ridge Area
	  else{
	    NYieldPhiEta2[0][MainLoop]+=binwx*binwy*hPhiEtaCorrSum->GetBinContent(x1,y1);
	    NYieldPhiEta2[1][MainLoop]+=pow(binwx*binwy*hPhiEtaCorrSum->GetBinError(x1,y1),2);
	    NYieldPhiEta2[2][MainLoop]+=binwx*binwy*hPhiEtaCorrPtSum->GetBinContent(x1,y1);
	    NYieldPhiEta2[3][MainLoop]+=pow(binwx*binwy*hPhiEtaCorrPtSum->GetBinError(x1,y1),2);
	    NYieldPhiEta2ZYAM[0][MainLoop]+=binwx*binwy*hPhiEtaCorr2Sum->GetBinContent(x1,y1);
	    NYieldPhiEta2ZYAM[1][MainLoop]+=pow(binwx*binwy*hPhiEtaCorr2Sum->GetBinError(x1,y1),2);
	    NYieldPhiEta2ZYAM[2][MainLoop]+=binwx*binwy*hPhiEtaCorrPt2Sum->GetBinContent(x1,y1);
	    NYieldPhiEta2ZYAM[3][MainLoop]+=pow(binwx*binwy*hPhiEtaCorrPt2Sum->GetBinError(x1,y1),2);
	    NearBinsPhiEta2++;
	  }
	}
	//Away-side 2D
	if(fabs(centx-Pi)<AwayWidthPhi||fabs(centx-3*Pi)<AwayWidthPhi||fabs(centx+Pi)<AwayWidthPhi){
	  AYieldPhiEta[0][MainLoop]+=binwx*binwy*hPhiEtaCorrSum->GetBinContent(x1,y1);
	  AYieldPhiEta[1][MainLoop]+=pow(binwx*binwy*hPhiEtaCorrSum->GetBinError(x1,y1),2);
	  AYieldPhiEta[2][MainLoop]+=binwx*binwy*hPhiEtaCorrPtSum->GetBinContent(x1,y1);
	  AYieldPhiEta[3][MainLoop]+=pow(binwx*binwy*hPhiEtaCorrPtSum->GetBinError(x1,y1),2);
	  AYieldPhiEtaZYAM[0][MainLoop]+=binwx*binwy*hPhiEtaCorr2Sum->GetBinContent(x1,y1);
	  AYieldPhiEtaZYAM[1][MainLoop]+=pow(binwx*binwy*hPhiEtaCorr2Sum->GetBinError(x1,y1),2);
	  AYieldPhiEtaZYAM[2][MainLoop]+=binwx*binwy*hPhiEtaCorrPt2Sum->GetBinContent(x1,y1);
	  AYieldPhiEtaZYAM[3][MainLoop]+=pow(binwx*binwy*hPhiEtaCorrPt2Sum->GetBinError(x1,y1),2);
	  AwayBinsPhiEta++;
	}
	//Minimum 2D
	if(fabs(centx-ZYAMCent)<ZYAMWidth||fabs(centx+ZYAMCent)<ZYAMWidth||fabs(centx-ZYAMCent+2*Pi)<ZYAMWidth||fabs(centx-ZYAMCent-2*Pi)<ZYAMWidth||fabs(centx+ZYAMCent+2*Pi)<ZYAMWidth||fabs(centx+ZYAMCent-2*Pi)<ZYAMWidth){
	  MYieldPhiEta[0][MainLoop]+=hPhiEtaCorrSum->GetBinContent(x1,y1);
	  MYieldPhiEta[1][MainLoop]+=pow(hPhiEtaCorrSum->GetBinError(x1,y1),2);
	  MYieldPhiEta[2][MainLoop]+=hPhiEtaCorrPtSum->GetBinContent(x1,y1);
	  MYieldPhiEta[3][MainLoop]+=pow(hPhiEtaCorrPtSum->GetBinError(x1,y1),2);
	  MinBinsPhiEta++;
	}//Min 2d
      }//y
    }//x

    //Averages
    
    cout << NYieldPhiZYAM[0][MainLoop] << endl;
    cout << AYieldPhiZYAM[0][MainLoop] << endl;
    cout << MYieldPhi[0][MainLoop] << endl;
    cout << NYieldEtaZYAM[0][MainLoop] << endl;
    cout << AYieldEtaZYAM[0][MainLoop] << endl;
    
    NAvePhi[0][MainLoop]=NYieldPhiZYAM[2][MainLoop]/NYieldPhiZYAM[0][MainLoop];
    AAvePhi[0][MainLoop]=AYieldPhiZYAM[2][MainLoop]/AYieldPhiZYAM[0][MainLoop];
    if(MYieldPhi[0][MainLoop]) MAvePhi[0][MainLoop]=MYieldPhi[2][MainLoop]/MYieldPhi[0][MainLoop];
    else MAvePhi[0][MainLoop]=0;
    NAveEta[0][MainLoop]=NYieldEtaZYAM[2][MainLoop]/NYieldEtaZYAM[0][MainLoop];
    AAveEta[0][MainLoop]=AYieldEtaZYAM[2][MainLoop]/AYieldEtaZYAM[0][MainLoop];
    NAvePhi[1][MainLoop]=NAvePhi[0][MainLoop]*pow(NYieldPhiZYAM[1][MainLoop]/pow(NYieldPhiZYAM[0][MainLoop],2)+NYieldPhiZYAM[3][MainLoop]/pow(NYieldPhiZYAM[2][MainLoop],2),0.5);
    AAvePhi[1][MainLoop]=AAvePhiZYAM[0][MainLoop]*pow(AYieldPhiZYAM[1][MainLoop]/pow(AYieldPhiZYAM[0][MainLoop],2)+AYieldPhiZYAM[3][MainLoop]/pow(AYieldPhiZYAM[2][MainLoop],2),0.5);
    if(MYieldPhi[0][MainLoop])MAvePhi[1][MainLoop]=MAvePhi[0][MainLoop]*pow(MYieldPhi[1][MainLoop]/pow(MYieldPhi[0][MainLoop],2)+MYieldPhi[3][MainLoop]/pow(MYieldPhi[2][MainLoop],2),0.5);
    else (MAvePhi[1][MainLoop])=1;
    NAveEta[1][MainLoop]=NAveEta[0][MainLoop]*pow(NYieldEtaZYAM[1][MainLoop]/pow(NYieldEtaZYAM[0][MainLoop],2)+NYieldEtaZYAM[3][MainLoop]/pow(NYieldEtaZYAM[2][MainLoop],2),0.5);
    AAveEta[1][MainLoop]=AAveEta[0][MainLoop]*pow(AYieldEtaZYAM[1][MainLoop]/pow(AYieldEtaZYAM[0][MainLoop],2)+AYieldEtaZYAM[3][MainLoop]/pow(AYieldEtaZYAM[2][MainLoop],2),0.5);

    //Set up the Error bars correctly
    NYieldPhi[1][MainLoop]=sqrt(NYieldPhi[1][MainLoop]);
    NYieldPhi[3][MainLoop]=sqrt(NYieldPhi[3][MainLoop]);
    AYieldPhi[1][MainLoop]=sqrt(AYieldPhi[1][MainLoop]);
    AYieldPhi[3][MainLoop]=sqrt(AYieldPhi[3][MainLoop]);
    NYieldPhiZYAM[1][MainLoop]=sqrt(NYieldPhiZYAM[1][MainLoop]);
    NYieldPhiZYAM[3][MainLoop]=sqrt(NYieldPhiZYAM[3][MainLoop]);
    AYieldPhiZYAM[1][MainLoop]=sqrt(AYieldPhiZYAM[1][MainLoop]);
    AYieldPhiZYAM[3][MainLoop]=sqrt(AYieldPhiZYAM[3][MainLoop]);
    MYieldPhi[1][MainLoop]=sqrt(MYieldPhi[1][MainLoop]);
    MYieldPhi[3][MainLoop]=sqrt(MYieldPhi[3][MainLoop]);

    NYieldEta[1][MainLoop]=sqrt(NYieldEta[1][MainLoop]);
    NYieldEta[3][MainLoop]=sqrt(NYieldEta[3][MainLoop]);
    AYieldEta[1][MainLoop]=sqrt(AYieldEta[1][MainLoop]);
    AYieldEta[3][MainLoop]=sqrt(AYieldEta[3][MainLoop]);
    NYieldEtaZYAM[1][MainLoop]=sqrt(NYieldEtaZYAM[1][MainLoop]);
    NYieldEtaZYAM[3][MainLoop]=sqrt(NYieldEtaZYAM[3][MainLoop]);
    AYieldEtaZYAM[1][MainLoop]=sqrt(AYieldEtaZYAM[1][MainLoop]);
    AYieldEtaZYAM[3][MainLoop]=sqrt(AYieldEtaZYAM[3][MainLoop]);
    
    //Set Minimums to total number of charged particles
    MYieldPhi[0][MainLoop]=MYieldPhi[0][MainLoop]/MinBinsPhi*2*Pi;
    MYieldPhi[1][MainLoop]=MYieldPhi[1][MainLoop]/MinBinsPhi*2*Pi;
    MYieldPhi[2][MainLoop]=MYieldPhi[2][MainLoop]/MinBinsPhi*2*Pi;
    MYieldPhi[3][MainLoop]=MYieldPhi[3][MainLoop]/MinBinsPhi*2*Pi;
    MYieldPhiEta[0][MainLoop]=MYieldPhiEta[0][MainLoop]/MinBinsPhiEta*2*Pi;
    MYieldPhiEta[1][MainLoop]=MYieldPhiEta[1][MainLoop]/MinBinsPhiEta*2*Pi;
    MYieldPhiEta[2][MainLoop]=MYieldPhiEta[2][MainLoop]/MinBinsPhiEta*2*Pi;
    MYieldPhiEta[3][MainLoop]=MYieldPhiEta[3][MainLoop]/MinBinsPhiEta*2*Pi;
  
    //0 (x1+x2)/2, 1 e2, 2 <x>,  3 e2, 4 lower error of 1 and 2, 5 e3, 6 zt1_0, 7 e5, 8 zt1_2, 9 e7, 10 zt1_3, 11 e9, 12 zt2_0, 13 e11, 14 zt2_1, 15 e13, 16 zt2_2, 17 e15  (zt_0 is (pt1+pt2)/2 zt_1 is <pt> zt_2 is lower error)
    NXPhi[0][MainLoop]=MainArray[MainLoop];
    NXPhi[1][MainLoop]=eMainArray[MainLoop];
    AXPhi[0][MainLoop]=MainArray[MainLoop];
    AXPhi[1][MainLoop]=eMainArray[MainLoop];
    NXPhiFit[0][MainLoop]=MainArray[MainLoop];
    NXPhiFit[1][MainLoop]=eMainArray[MainLoop];
    AXPhiFit[0][MainLoop]=MainArray[MainLoop];
    AXPhiFit[1][MainLoop]=eMainArray[MainLoop];

    NXEta[0][MainLoop]=MainArray[MainLoop];
    NXEta[1][MainLoop]=eMainArray[MainLoop];
    AXEta[0][MainLoop]=MainArray[MainLoop];
    AXEta[1][MainLoop]=eMainArray[MainLoop];
    NXEtaFit[0][MainLoop]=MainArray[MainLoop];
    NXEtaFit[1][MainLoop]=eMainArray[MainLoop];
    AXEtaFit[0][MainLoop]=MainArray[MainLoop];
    AXEtaFit[1][MainLoop]=eMainArray[MainLoop];

    NXPhiEta1[0][MainLoop]=MainArray[MainLoop];
    NXPhiEta1[1][MainLoop]=eMainArray[MainLoop];
    NXPhiEta2[0][MainLoop]=MainArray[MainLoop];
    NXPhiEta2[1][MainLoop]=eMainArray[MainLoop];
    AXPhiEta[0][MainLoop]=MainArray[MainLoop];
    AXPhiEta[1][MainLoop]=eMainArray[MainLoop];
  
    if(APtTPtMult==0){//spectra so <pt_assoc>
      if(NYieldPhiZYAM[0][MainLoop]==0) NYieldPhiZYAM[0][MainLoop]=0.00000001;
      if(NYieldPhiZYAM[2][MainLoop]==0) NYieldPhiZYAM[2][MainLoop]=0.00000001;
      if(AYieldPhiZYAM[0][MainLoop]==0) AYieldPhiZYAM[0][MainLoop]=0.00000001;
      if(AYieldPhiZYAM[2][MainLoop]==0) AYieldPhiZYAM[2][MainLoop]=0.00000001;
      NXPhi[2][MainLoop]=NYieldPhiZYAM[2][MainLoop]/NYieldPhiZYAM[0][MainLoop];
      NXPhi[3][MainLoop]=NXPhi[2][MainLoop]*sqrt(pow(NYieldPhiZYAM[1][MainLoop]/NYieldPhiZYAM[0][MainLoop],2)+pow(NYieldPhiZYAM[3][MainLoop]/NYieldPhiZYAM[2][MainLoop],2));
      AXPhi[2][MainLoop]=AYieldPhiZYAM[2][MainLoop]/AYieldPhiZYAM[0][MainLoop];
      AXPhi[3][MainLoop]=AXPhi[2][MainLoop]*sqrt(pow(AYieldPhiZYAM[1][MainLoop]/AYieldPhiZYAM[0][MainLoop],2)+pow(AYieldPhiZYAM[3][MainLoop]/AYieldPhiZYAM[2][MainLoop],2));
      NXPhiFit[2][MainLoop]=NYieldPhiFit[2][MainLoop]/NYieldPhiFit[0][MainLoop];
      // cout << "Pt " << NXPhiFit[2][MainLoop] << endl;
      NXPhiFit[3][MainLoop]=NXPhiFit[2][MainLoop]*sqrt(pow(NYieldPhiFit[1][MainLoop]/NYieldPhiFit[0][MainLoop],2)+pow(NYieldPhiFit[3][MainLoop]/NYieldPhiFit[2][MainLoop],2));
      AXPhiFit[2][MainLoop]=AYieldPhiFit[2][MainLoop]/AYieldPhiFit[0][MainLoop];
      AXPhiFit[3][MainLoop]=AXPhiFit[2][MainLoop]*sqrt(pow(AYieldPhiFit[1][MainLoop]/AYieldPhiFit[0][MainLoop],2)+pow(AYieldPhiFit[3][MainLoop]/AYieldPhiFit[2][MainLoop],2));

      NXEta[2][MainLoop]=NYieldEtaZYAM[2][MainLoop]/NYieldEtaZYAM[0][MainLoop];
      NXEta[3][MainLoop]=NXEta[2][MainLoop]*sqrt(pow(NYieldEtaZYAM[1][MainLoop]/NYieldEtaZYAM[0][MainLoop],2)+pow(NYieldEtaZYAM[3][MainLoop]/NYieldEtaZYAM[2][MainLoop],2));
      AXEta[2][MainLoop]=AYieldEtaZYAM[2][MainLoop]/AYieldEtaZYAM[0][MainLoop];
      AXEta[3][MainLoop]=AXEta[2][MainLoop]*sqrt(pow(AYieldEtaZYAM[1][MainLoop]/AYieldEtaZYAM[0][MainLoop],2)+pow(AYieldEtaZYAM[3][MainLoop]/AYieldEtaZYAM[2][MainLoop],2));
      NXEtaFit[2][MainLoop]=NYieldEtaFit[2][MainLoop]/NYieldEtaFit[0][MainLoop];
      NXEtaFit[3][MainLoop]=NXEtaFit[2][MainLoop]*sqrt(pow(NYieldEtaFit[1][MainLoop]/NYieldEtaFit[0][MainLoop],2)+pow(NYieldEtaFit[3][MainLoop]/NYieldEtaFit[2][MainLoop],2));
      AXEtaFit[2][MainLoop]=AYieldEtaFit[2][MainLoop]/AYieldEtaFit[0][MainLoop];
      AXEtaFit[3][MainLoop]=AXEtaFit[2][MainLoop]*sqrt(pow(AYieldEtaFit[1][MainLoop]/AYieldEtaFit[0][MainLoop],2)+pow(AYieldEtaFit[3][MainLoop]/AYieldEtaFit[2][MainLoop],2));

      
      /*
	NXPhiEta1[2][MainLoop]=NYieldPhiEta1[2][MainLoop]/NYieldPhiEta1[0][MainLoop];
	NXPhiEta1[3][MainLoop]=NXPhiEta1[2][MainLoop][0]*sqrt(pow(NYieldPhiEta1[1][MainLoop]/NYieldPhiEta1[0][MainLoop],2)+pow(NYieldPhiEta1[3][MainLoop]/NYieldPhiEta1[2][MainLoop],2));
	NXPhiEta2[2][MainLoop]=NYieldPhiEta2[2][MainLoop]/NYieldPhiEta2[0][MainLoop];
	NXPhiEta2[3][MainLoop]=NXPhiEta2[2][MainLoop][0]*sqrt(pow(NYieldPhiEta2[1][MainLoop]/NYieldPhiEta2[0][MainLoop],2)+pow(NYieldPhiEta2[3][MainLoop]/NYieldPhiEta2[2][MainLoop],2));
	AXPhiEta[2][MainLoop]=AYieldPhiEta[2][MainLoop]/AYieldPhiEta[0][MainLoop];
	AXPhiEta[3][MainLoop]=AXPhiEta[2][MainLoop][0]*sqrt(pow(AYieldPhiEta[1][MainLoop]/AYieldPhiEta[0][MainLoop],2)+pow(AYieldPhiEta[3][MainLoop]/AYieldPhiEta[2][MainLoop],2));
      */
   
      if(NXPhi[1][MainLoop]<NXPhi[3][MainLoop]){
	NXPhi[4][MainLoop]=NXPhi[0][MainLoop];
	NXPhi[5][MainLoop]=NXPhi[1][MainLoop];
      }
      else{
	NXPhi[4][MainLoop]=NXPhi[2][MainLoop];
	NXPhi[5][MainLoop]=NXPhi[3][MainLoop];
      } 
      if(AXPhi[1][MainLoop]<AXPhi[3][MainLoop]){
	AXPhi[4][MainLoop]=AXPhi[0][MainLoop];
	AXPhi[5][MainLoop]=AXPhi[1][MainLoop];
      }
      else{
	AXPhi[4][MainLoop]=AXPhi[2][MainLoop];
	AXPhi[5][MainLoop]=AXPhi[3][MainLoop];
      }	
      if(NXPhiFit[1][MainLoop]<NXPhiFit[3][MainLoop]){
	NXPhiFit[4][MainLoop]=NXPhiFit[0][MainLoop];
	NXPhiFit[5][MainLoop]=NXPhiFit[1][MainLoop];
      }
      else{
	NXPhiFit[4][MainLoop]=NXPhiFit[2][MainLoop];
	NXPhiFit[5][MainLoop]=NXPhiFit[3][MainLoop];
      } 
      if(AXPhiFit[1][MainLoop]<AXPhiFit[3][MainLoop]){
	AXPhiFit[4][MainLoop]=AXPhiFit[0][MainLoop];
	AXPhiFit[5][MainLoop]=AXPhiFit[1][MainLoop];
      }
      else{
	AXPhiFit[4][MainLoop]=AXPhiFit[2][MainLoop];
	AXPhiFit[5][MainLoop]=AXPhiFit[3][MainLoop];
      }	

      if(NXEta[1][MainLoop]<NXEta[3][MainLoop]){
	NXEta[4][MainLoop]=NXEta[0][MainLoop];
	NXEta[5][MainLoop]=NXEta[1][MainLoop];
      }
      else{
	NXEta[4][MainLoop]=NXEta[2][MainLoop];
	NXEta[5][MainLoop]=NXEta[3][MainLoop];
      } 
      if(AXEta[1][MainLoop]<AXEta[3][MainLoop]){
	AXEta[4][MainLoop]=AXEta[0][MainLoop];
	AXEta[5][MainLoop]=AXEta[1][MainLoop];
      }
      else{
	AXEta[4][MainLoop]=AXEta[2][MainLoop];
	AXEta[5][MainLoop]=AXEta[3][MainLoop];
      }	
      if(NXEtaFit[1][MainLoop]<NXEtaFit[3][MainLoop]){
	NXEtaFit[4][MainLoop]=NXEtaFit[0][MainLoop];
	NXEtaFit[5][MainLoop]=NXEtaFit[1][MainLoop];
      }
      else{
	NXEtaFit[4][MainLoop]=NXEtaFit[2][MainLoop];
	NXEtaFit[5][MainLoop]=NXEtaFit[3][MainLoop];
      } 
      if(AXEtaFit[1][MainLoop]<AXEtaFit[3][MainLoop]){
	AXEtaFit[4][MainLoop]=AXEtaFit[0][MainLoop];
	AXEtaFit[5][MainLoop]=AXEtaFit[1][MainLoop];
      }
      else{
	AXEtaFit[4][MainLoop]=AXEtaFit[2][MainLoop];
	AXEtaFit[5][MainLoop]=AXEtaFit[3][MainLoop];
      }	
      
      for(int tloop=0;tloop<=2; tloop++){
	NXPhi[6+2*tloop][MainLoop]=NXPhi[0+2*tloop][MainLoop]/MPt[0];
	NXPhi[7+2*tloop][MainLoop]=NXPhi[6+2*tloop][MainLoop]*sqrt(pow(NXPhi[1+2*tloop][MainLoop]/NXPhi[0+2*tloop][MainLoop],2)+pow(MPt[1]/MPt[0],2));
	NXPhi[12+2*tloop][MainLoop]=NXPhi[0+2*tloop][MainLoop]/(MPt[0]+TNYieldPhiZYAM[2][MainLoop]);
	NXPhi[13+2*tloop][MainLoop]=NXPhi[12+2*tloop][MainLoop]*sqrt(pow(NXPhi[1+2*tloop][MainLoop]/NXPhi[0+2*tloop][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldPhiZYAM[3][MainLoop])/pow(MPt[0]+TNYieldPhiZYAM[2][MainLoop],2));
	AXPhi[6+2*tloop][MainLoop]=AXPhi[0+2*tloop][MainLoop]/MPt[0];
	AXPhi[7+2*tloop][MainLoop]=AXPhi[6+2*tloop][MainLoop]*sqrt(pow(AXPhi[1+2*tloop][MainLoop]/AXPhi[0+2*tloop][MainLoop],2)+pow(MPt[1]/MPt[0],2));
	AXPhi[12+2*tloop][MainLoop]=AXPhi[0+2*tloop][MainLoop]/(MPt[0]+TNYieldPhiZYAM[2][MainLoop]);
	AXPhi[13+2*tloop][MainLoop]=AXPhi[12+2*tloop][MainLoop]*sqrt(pow(AXPhi[1+2*tloop][MainLoop]/AXPhi[0+2*tloop][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldPhiZYAM[3][MainLoop])/pow(MPt[0]+TNYieldPhiZYAM[2][MainLoop],2));

	NXPhiFit[6+2*tloop][MainLoop]=NXPhiFit[0+2*tloop][MainLoop]/MPt[0];
	NXPhiFit[7+2*tloop][MainLoop]=NXPhiFit[6+2*tloop][MainLoop]*sqrt(pow(NXPhiFit[1+2*tloop][MainLoop]/NXPhiFit[0+2*tloop][MainLoop],2)+pow(MPt[1]/MPt[0],2));
	//	NXPhi[12+2*tloop][MainLoop]=NXPhi[0+2*tloop][MainLoop]/(MPt[0]+TNYieldPhiZYAM[2][MainLoop]);
	//	NXPhi[13+2*tloop][MainLoop]=NXPhi[12+2*tloop][MainLoop]*sqrt(pow(NXPhi[1+2*tloop][MainLoop]/NXPhi[0+2*tloop][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldPhiZYAM[3][MainLoop])/pow(MPt[0]+TNYieldPhiZYAM[2][MainLoop],2));
	AXPhiFit[6+2*tloop][MainLoop]=AXPhiFit[0+2*tloop][MainLoop]/MPt[0];
	AXPhiFit[7+2*tloop][MainLoop]=AXPhiFit[6+2*tloop][MainLoop]*sqrt(pow(AXPhiFit[1+2*tloop][MainLoop]/AXPhiFit[0+2*tloop][MainLoop],2)+pow(MPt[1]/MPt[0],2));
	//	AXPhi[12+2*tloop][MainLoop]=AXPhi[0+2*tloop][MainLoop]/(MPt[0]+TNYieldPhiZYAM[2][MainLoop]);
	//	AXPhi[13+2*tloop][MainLoop]=AXPhi[12+2*tloop][MainLoop]*sqrt(pow(AXPhi[1+2*tloop][MainLoop]/AXPhi[0+2*tloop][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldPhiZYAM[3][MainLoop])/pow(MPt[0]+TNYieldPhiZYAM[2][MainLoop],2));

	NXEta[6+2*tloop][MainLoop]=NXEta[0+2*tloop][MainLoop]/MPt[0];
	NXEta[7+2*tloop][MainLoop]=NXEta[6+2*tloop][MainLoop]*sqrt(pow(NXEta[1+2*tloop][MainLoop]/NXEta[0+2*tloop][MainLoop],2)+pow(MPt[1]/MPt[0],2));
	//	NXEta[12+2*tloop][MainLoop]=NXEta[0+2*tloop][MainLoop]/(MPt[0]+TNYieldEtaZYAM[2][MainLoop]);
	//	NXEta[13+2*tloop][MainLoop]=NXEta[12+2*tloop][MainLoop]*sqrt(pow(NXEta[1+2*tloop][MainLoop]/NXEta[0+2*tloop][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldEtaZYAM[3][MainLoop])/pow(MPt[0]+TNYieldEtaZYAM[2][MainLoop],2));
	AXEta[6+2*tloop][MainLoop]=AXEta[0+2*tloop][MainLoop]/MPt[0];
	AXEta[7+2*tloop][MainLoop]=AXEta[6+2*tloop][MainLoop]*sqrt(pow(AXEta[1+2*tloop][MainLoop]/AXEta[0+2*tloop][MainLoop],2)+pow(MPt[1]/MPt[0],2));
	//AXEta[12+2*tloop][MainLoop]=AXEta[0+2*tloop][MainLoop]/(MPt[0]+TNYieldEtaZYAM[2][MainLoop]);
	//	AXEta[13+2*tloop][MainLoop]=AXEta[12+2*tloop][MainLoop]*sqrt(pow(AXEta[1+2*tloop][MainLoop]/AXEta[0+2*tloop][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldEtaZYAM[3][MainLoop])/pow(MPt[0]+TNYieldEtaZYAM[2][MainLoop],2));

	NXEtaFit[6+2*tloop][MainLoop]=NXEtaFit[0+2*tloop][MainLoop]/MPt[0];
	NXEtaFit[7+2*tloop][MainLoop]=NXEtaFit[6+2*tloop][MainLoop]*sqrt(pow(NXEtaFit[1+2*tloop][MainLoop]/NXEtaFit[0+2*tloop][MainLoop],2)+pow(MPt[1]/MPt[0],2));
	AXEtaFit[6+2*tloop][MainLoop]=AXEtaFit[0+2*tloop][MainLoop]/MPt[0];
	AXEtaFit[7+2*tloop][MainLoop]=AXEtaFit[6+2*tloop][MainLoop]*sqrt(pow(AXEtaFit[1+2*tloop][MainLoop]/AXEtaFit[0+2*tloop][MainLoop],2)+pow(MPt[1]/MPt[0],2));

      }
    
      NYieldPhi[4][MainLoop]=NYieldPhi[0][MainLoop]/NXPhi[1][MainLoop];
      NYieldPhi[5][MainLoop]=NYieldPhi[1][MainLoop]/NXPhi[1][MainLoop];
      NYieldPhi[6][MainLoop]=NYieldPhi[0][MainLoop]/NXPhi[1][MainLoop]*MPt[0];
      NYieldPhi[7][MainLoop]=NYieldPhi[6][MainLoop]*sqrt(pow(NYieldPhi[1][MainLoop]/NYieldPhi[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));
      NYieldPhi[8][MainLoop]=NYieldPhi[0][MainLoop]/NXPhi[1][MainLoop]*(MPt[0]+TNYieldPhi[2][MainLoop]);
      NYieldPhi[9][MainLoop]=NYieldPhi[8][MainLoop]*sqrt(pow(NYieldPhi[1][MainLoop]/NYieldPhi[0][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldPhi[3][MainLoop]*TNYieldPhi[3][MainLoop])/pow(MPt[0]+TNYieldPhi[2][MainLoop],2));
      NYieldPhiFit[7][MainLoop]=NYieldPhiFit[6][MainLoop]*sqrt(pow(NYieldPhiFit[1][MainLoop]/NYieldPhiFit[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));
      NYieldPhiZYAM[4][MainLoop]=NYieldPhiZYAM[0][MainLoop]/NXPhi[1][MainLoop];
      NYieldPhiZYAM[5][MainLoop]=NYieldPhiZYAM[1][MainLoop]/NXPhi[1][MainLoop];
      NYieldPhiZYAM[6][MainLoop]=NYieldPhiZYAM[0][MainLoop]/NXPhi[1][MainLoop]*MPt[0];
      NYieldPhiZYAM[7][MainLoop]=NYieldPhiZYAM[6][MainLoop]*sqrt(pow(NYieldPhiZYAM[1][MainLoop]/NYieldPhiZYAM[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));
      NYieldPhiZYAM[8][MainLoop]=NYieldPhiZYAM[0][MainLoop]/NXPhi[1][MainLoop]*(MPt[0]+TNYieldPhi[2][MainLoop]);
      NYieldPhiZYAM[9][MainLoop]=NYieldPhiZYAM[8][MainLoop]*sqrt(pow(NYieldPhiZYAM[1][MainLoop]/NYieldPhiZYAM[0][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldPhi[3][MainLoop]*TNYieldPhi[3][MainLoop])/pow(MPt[0]+TNYieldPhi[2][MainLoop],2));
      AYieldPhiZYAM[4][MainLoop]=AYieldPhiZYAM[0][MainLoop]/NXPhi[1][MainLoop];
      AYieldPhiZYAM[5][MainLoop]=AYieldPhiZYAM[1][MainLoop]/NXPhi[1][MainLoop];
      AYieldPhiZYAM[6][MainLoop]=AYieldPhiZYAM[0][MainLoop]/NXPhi[1][MainLoop]*MPt[0];
      AYieldPhiZYAM[7][MainLoop]=AYieldPhiZYAM[6][MainLoop]*sqrt(pow(NYieldPhiZYAM[1][MainLoop]/NYieldPhiZYAM[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));
      AYieldPhiZYAM[8][MainLoop]=AYieldPhiZYAM[0][MainLoop]/NXPhi[1][MainLoop]*(MPt[0]+TNYieldPhi[2][MainLoop]);
      AYieldPhiZYAM[9][MainLoop]=AYieldPhiZYAM[8][MainLoop]*sqrt(pow(NYieldPhiZYAM[1][MainLoop]/NYieldPhiZYAM[0][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldPhi[3][MainLoop]*TNYieldPhi[3][MainLoop])/pow(MPt[0]+TNYieldPhi[2][MainLoop],2));
      
      NYieldPhiFit[4][MainLoop]=NYieldPhiFit[0][MainLoop]/NXPhiFit[1][MainLoop];
      NYieldPhiFit[5][MainLoop]=NYieldPhiFit[1][MainLoop]/NXPhiFit[1][MainLoop];
      NYieldPhiFit[6][MainLoop]=NYieldPhiFit[0][MainLoop]/NXPhiFit[1][MainLoop]*MPt[0];
      AYieldPhiFit[4][MainLoop]=AYieldPhiFit[0][MainLoop]/AXPhiFit[1][MainLoop];
      AYieldPhiFit[5][MainLoop]=AYieldPhiFit[1][MainLoop]/AXPhiFit[1][MainLoop];
      AYieldPhiFit[6][MainLoop]=AYieldPhiFit[0][MainLoop]/AXPhiFit[1][MainLoop]*MPt[0];
      AYieldPhiFit[7][MainLoop]=AYieldPhiFit[6][MainLoop]*sqrt(pow(NYieldPhiFit[1][MainLoop]/NYieldPhiFit[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));

      NYieldEta[4][MainLoop]=NYieldEta[0][MainLoop]/NXEta[1][MainLoop];
      NYieldEta[5][MainLoop]=NYieldEta[1][MainLoop]/NXEta[1][MainLoop];
      NYieldEta[6][MainLoop]=NYieldEta[0][MainLoop]/NXEta[1][MainLoop]*MPt[0];
      NYieldEta[7][MainLoop]=NYieldEta[6][MainLoop]*sqrt(pow(NYieldEta[1][MainLoop]/NYieldEta[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));
      // NYieldEta[8][MainLoop]=NYieldEta[0][MainLoop]/NXEta[1][MainLoop]*(MPt[0]+TNYieldEta[2][MainLoop]);
      //  NYieldEta[9][MainLoop]=NYieldEta[8][MainLoop]*sqrt(pow(NYieldEta[1][MainLoop]/NYieldEta[0][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldEta[3][MainLoop]*TNYieldEta[3][MainLoop])/pow(MPt[0]+TNYieldEta[2][MainLoop],2));
      NYieldEtaFit[7][MainLoop]=NYieldEtaFit[6][MainLoop]*sqrt(pow(NYieldEtaFit[1][MainLoop]/NYieldEtaFit[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));
      NYieldEtaZYAM[4][MainLoop]=NYieldEtaZYAM[0][MainLoop]/NXEta[1][MainLoop];
      NYieldEtaZYAM[5][MainLoop]=NYieldEtaZYAM[1][MainLoop]/NXEta[1][MainLoop];
      NYieldEtaZYAM[6][MainLoop]=NYieldEtaZYAM[0][MainLoop]/NXEta[1][MainLoop]*MPt[0];
      NYieldEtaZYAM[7][MainLoop]=NYieldEtaZYAM[6][MainLoop]*sqrt(pow(NYieldEtaZYAM[1][MainLoop]/NYieldEtaZYAM[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));
      // NYieldEtaZYAM[8][MainLoop]=NYieldEtaZYAM[0][MainLoop]/NXEta[1][MainLoop]*(MPt[0]+TNYieldEta[2][MainLoop]);
      //  NYieldEtaZYAM[9][MainLoop]=NYieldEtaZYAM[8][MainLoop]*sqrt(pow(NYieldEtaZYAM[1][MainLoop]/NYieldEtaZYAM[0][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldEta[3][MainLoop]*TNYieldEta[3][MainLoop])/pow(MPt[0]+TNYieldEta[2][MainLoop],2));
      AYieldEtaZYAM[4][MainLoop]=AYieldEtaZYAM[0][MainLoop]/NXEta[1][MainLoop];
      AYieldEtaZYAM[5][MainLoop]=AYieldEtaZYAM[1][MainLoop]/NXEta[1][MainLoop];
      AYieldEtaZYAM[6][MainLoop]=AYieldEtaZYAM[0][MainLoop]/NXEta[1][MainLoop]*MPt[0];
      AYieldEtaZYAM[7][MainLoop]=AYieldEtaZYAM[6][MainLoop]*sqrt(pow(NYieldEtaZYAM[1][MainLoop]/NYieldEtaZYAM[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));
      //  AYieldEtaZYAM[8][MainLoop]=AYieldEtaZYAM[0][MainLoop]/NXEta[1][MainLoop]*(MPt[0]+TNYieldEta[2][MainLoop]);
      //  AYieldEtaZYAM[9][MainLoop]=AYieldEtaZYAM[8][MainLoop]*sqrt(pow(NYieldEtaZYAM[1][MainLoop]/NYieldEtaZYAM[0][MainLoop],2)+(MPt[1]*MPt[1]+TNYieldEta[3][MainLoop]*TNYieldEta[3][MainLoop])/pow(MPt[0]+TNYieldEta[2][MainLoop],2));
      
      NYieldEtaFit[4][MainLoop]=NYieldEtaFit[0][MainLoop]/NXEtaFit[1][MainLoop];
      NYieldEtaFit[5][MainLoop]=NYieldEtaFit[1][MainLoop]/NXEtaFit[1][MainLoop];
      NYieldEtaFit[6][MainLoop]=NYieldEtaFit[0][MainLoop]/NXEtaFit[1][MainLoop]*MPt[0];
      AYieldEtaFit[4][MainLoop]=AYieldEtaFit[0][MainLoop]/AXEtaFit[1][MainLoop];
      AYieldEtaFit[5][MainLoop]=AYieldEtaFit[1][MainLoop]/AXEtaFit[1][MainLoop];
      AYieldEtaFit[6][MainLoop]=AYieldEtaFit[0][MainLoop]/AXEtaFit[1][MainLoop]*MPt[0];
      AYieldEtaFit[7][MainLoop]=AYieldEtaFit[6][MainLoop]*sqrt(pow(NYieldEtaFit[1][MainLoop]/NYieldEtaFit[0][MainLoop],2)+pow(MPt[1]/MPt[0],2));
    } 
  
    if(APtTPtMult==1){//Trigger pt dependence so <pt_Trig>
      NXPhi[2][MainLoop]=MPt[0];
      NXPhi[3][MainLoop]=MPt[1];
      AXPhi[2][MainLoop]=MPt[0];
      AXPhi[3][MainLoop]=MPt[1];
      NXPhiFit[2][MainLoop]=MPt[0];
      NXPhiFit[3][MainLoop]=MPt[1];
      AXPhiFit[2][MainLoop]=MPt[0];
      AXPhiFit[3][MainLoop]=MPt[1];
      
      NXEta[2][MainLoop]=MPt[0];
      NXEta[3][MainLoop]=MPt[1];
      AXEta[2][MainLoop]=MPt[0];
      AXEta[3][MainLoop]=MPt[1];
      NXEtaFit[2][MainLoop]=MPt[0];
      NXEtaFit[3][MainLoop]=MPt[1];
      AXEtaFit[2][MainLoop]=MPt[0];
      AXEtaFit[3][MainLoop]=MPt[1];

      NXPhiEta1[2][MainLoop]=MPt[0];
      NXPhiEta1[3][MainLoop]=MPt[1];
      NXPhiEta2[2][MainLoop]=MPt[0];
      NXPhiEta2[3][MainLoop]=MPt[1];
      AXPhiEta[2][MainLoop]=MPt[0];
      AXPhiEta[3][MainLoop]=MPt[1];

      if(NXPhi[1][MainLoop]<NXPhi[3][MainLoop]){
	NXPhi[4][MainLoop]=NXPhi[0][MainLoop];
	NXPhi[5][MainLoop]=NXPhi[1][MainLoop];
      }
      else{
	NXPhi[4][MainLoop]=NXPhi[2][MainLoop];
	NXPhi[5][MainLoop]=NXPhi[3][MainLoop];
      }
      if(AXPhi[1][MainLoop]<AXPhi[3][MainLoop]){
	AXPhi[4][MainLoop]=AXPhi[0][MainLoop];
	AXPhi[5][MainLoop]=AXPhi[1][MainLoop];
      }
      else{
	AXPhi[4][MainLoop]=AXPhi[2][MainLoop];
	AXPhi[5][MainLoop]=AXPhi[3][MainLoop];
      }	
      if(NXPhiFit[1][MainLoop]<NXPhiFit[3][MainLoop]){
	NXPhiFit[4][MainLoop]=NXPhiFit[0][MainLoop];
	NXPhiFit[5][MainLoop]=NXPhiFit[1][MainLoop];
      }
      else{
	NXPhiFit[4][MainLoop]=NXPhiFit[2][MainLoop];
	NXPhiFit[5][MainLoop]=NXPhiFit[3][MainLoop];
      }
      if(AXPhiFit[1][MainLoop]<AXPhiFit[3][MainLoop]){
	AXPhiFit[4][MainLoop]=AXPhiFit[0][MainLoop];
	AXPhiFit[5][MainLoop]=AXPhiFit[1][MainLoop];
      }
      else{
	AXPhiFit[4][MainLoop]=AXPhiFit[2][MainLoop];
	AXPhiFit[5][MainLoop]=AXPhiFit[3][MainLoop];
      }	

      if(NXEta[1][MainLoop]<NXEta[3][MainLoop]){
	NXEta[4][MainLoop]=NXEta[0][MainLoop];
	NXEta[5][MainLoop]=NXEta[1][MainLoop];
      }
      else{
	NXEta[4][MainLoop]=NXEta[2][MainLoop];
	NXEta[5][MainLoop]=NXEta[3][MainLoop];
      }
      if(AXEta[1][MainLoop]<AXEta[3][MainLoop]){
	AXEta[4][MainLoop]=AXEta[0][MainLoop];
	AXEta[5][MainLoop]=AXEta[1][MainLoop];
      }
      else{
	AXEta[4][MainLoop]=AXEta[2][MainLoop];
	AXEta[5][MainLoop]=AXEta[3][MainLoop];
      }	
      if(NXEtaFit[1][MainLoop]<NXEtaFit[3][MainLoop]){
	NXEtaFit[4][MainLoop]=NXEtaFit[0][MainLoop];
	NXEtaFit[5][MainLoop]=NXEtaFit[1][MainLoop];
      }
      else{
	NXEtaFit[4][MainLoop]=NXEtaFit[2][MainLoop];
	NXEtaFit[5][MainLoop]=NXEtaFit[3][MainLoop];
      }
      if(AXEtaFit[1][MainLoop]<AXEtaFit[3][MainLoop]){
	AXEtaFit[4][MainLoop]=AXEtaFit[0][MainLoop];
	AXEtaFit[5][MainLoop]=AXEtaFit[1][MainLoop];
      }
      else{
	AXEtaFit[4][MainLoop]=AXEtaFit[2][MainLoop];
	AXEtaFit[5][MainLoop]=AXEtaFit[3][MainLoop];
      }	

      
    }
    if(APtTPtMult==2){//Centrality so Mult
      NXPhi[2][MainLoop]=0;
      NXPhi[3][MainLoop]=0;
      tempsum=0;
      for(int multx=(MultArray1[Mult1]+1);multx<(MultArray2[Mult2]+1);multx++){
	NXPhi[2][MainLoop]+=hMult->GetBinCenter(multx)*(multx-1);
	NXPhi[3][MainLoop]+=pow(hMult->GetBinError(multx),2);
	tempsum+=hMult->GetBinCenter(multx);
      }
      NXPhi[2][MainLoop]/=tempsum;
      NXPhi[3][MainLoop]=sqrt(NXPhi[3][MainLoop])/tempsum;
      AXPhi[2][MainLoop]=NXPhi[2][MainLoop];
      AXPhi[3][MainLoop]=NXPhi[3][MainLoop];
      NXPhiFit[2][MainLoop]=NXPhi[2][MainLoop];
      NXPhiFit[3][MainLoop]=NXPhi[3][MainLoop];
      AXPhiFit[2][MainLoop]=AXPhi[2][MainLoop];
      AXPhiFit[3][MainLoop]=AXPhi[3][MainLoop];
      
      NXPhiEta1[2][MainLoop]=NXPhi[2][MainLoop];
      NXPhiEta1[3][MainLoop]=NXPhi[3][MainLoop];
      NXPhiEta2[2][MainLoop]=NXPhi[2][MainLoop];
      NXPhiEta2[3][MainLoop]=NXPhi[3][MainLoop];
      AXPhiEta[2][MainLoop]=NXPhi[2][MainLoop];
      AXPhiEta[3][MainLoop]=NXPhi[2][MainLoop];
    }
    
  }//Sould end mainloop
  //right now only finish up normal spectra and zT with assoc/trigger pt maybe do the assoc/near-side pt auch
  //0 (x1+x2)/2, 1 e2, 2 <x>,  3 e2, 4 lower error of 1 and 2, 5 e3, 6 zt1_0, 7 e5, 8 zt1_2, 9 e7, 10 zt1_3, 11 e9, 12 zt2_0, 13 e11, 14 zt2_1, 15 e13, 16 zt2_2, 17 e15  (zt_0 is (pt1+pt2)/2 zt_1 is <pt> zt_2 is lower error)
  //0 yield, 1 eyield, 2 pt yield, 3 ept yield, 4 dN/dpt 5 e4 6 dN/dzt (zt=Trig) yield, 7 ezt yield, 8 dN/dzt (zt=near) yield, 9 ezt2 yield
  if(APtTPtMult==0){
    if(HorizontalErrors==0){
      for(int mloop=0;mloop<nMainLoop;mloop++){
	NXPhi[1][mloop]=0;
	NXPhi[3][mloop]=0;
	NXPhi[5][mloop]=0;
	NXPhi[7][mloop]=0;
	NXPhi[9][mloop]=0;
	NXPhi[11][mloop]=0;
	NXPhi[13][mloop]=0;
	NXPhi[15][mloop]=0;
	NXPhi[17][mloop]=0;
	AXPhi[1][mloop]=0;
	AXPhi[3][mloop]=0;
	AXPhi[5][mloop]=0;
	AXPhi[7][mloop]=0;
	AXPhi[9][mloop]=0;
	AXPhi[11][mloop]=0;
	AXPhi[13][mloop]=0;
	AXPhi[15][mloop]=0;
	AXPhi[17][mloop]=0;
	NXPhiFit[1][mloop]=0;
	NXPhiFit[3][mloop]=0;
	NXPhiFit[5][mloop]=0;
	NXPhiFit[7][mloop]=0;
	NXPhiFit[9][mloop]=0;
	NXPhiFit[11][mloop]=0;
	NXPhiFit[13][mloop]=0;
	NXPhiFit[15][mloop]=0;
	NXPhiFit[17][mloop]=0;
	AXPhiFit[1][mloop]=0;
	AXPhiFit[3][mloop]=0;
	AXPhiFit[5][mloop]=0;
	AXPhiFit[7][mloop]=0;
	AXPhiFit[9][mloop]=0;
	AXPhiFit[11][mloop]=0;
	AXPhiFit[13][mloop]=0;
	AXPhiFit[15][mloop]=0;
	AXPhiFit[17][mloop]=0;
	NXEta[1][mloop]=0;
	NXEta[3][mloop]=0;
	NXEta[5][mloop]=0;
	NXEta[7][mloop]=0;
	NXEta[9][mloop]=0;
	NXEta[11][mloop]=0;
	NXEta[13][mloop]=0;
	NXEta[15][mloop]=0;
	NXEta[17][mloop]=0;
	AXEta[1][mloop]=0;
	AXEta[3][mloop]=0;
	AXEta[5][mloop]=0;
	AXEta[7][mloop]=0;
	AXEta[9][mloop]=0;
	AXEta[11][mloop]=0;
	AXEta[13][mloop]=0;
	AXEta[15][mloop]=0;
	AXEta[17][mloop]=0;
	NXEtaFit[1][mloop]=0;
	NXEtaFit[3][mloop]=0;
	NXEtaFit[5][mloop]=0;
	NXEtaFit[7][mloop]=0;
	NXEtaFit[9][mloop]=0;
	NXEtaFit[11][mloop]=0;
	NXEtaFit[13][mloop]=0;
	NXEtaFit[15][mloop]=0;
	NXEtaFit[17][mloop]=0;
	AXEtaFit[1][mloop]=0;
	AXEtaFit[3][mloop]=0;
	AXEtaFit[5][mloop]=0;
	AXEtaFit[7][mloop]=0;
	AXEtaFit[9][mloop]=0;
	AXEtaFit[11][mloop]=0;
	AXEtaFit[13][mloop]=0;
	AXEtaFit[15][mloop]=0;
	AXEtaFit[17][mloop]=0;
      }
    }
    gNPhiSpectraZYAM1=new TGraphErrors(nMainLoop,NXPhi[0],NYieldPhiZYAM[4],NXPhi[1],NYieldPhiZYAM[5]); 
    gNPhiSpectraZYAM2=new TGraphErrors(nMainLoop,NXPhi[2],NYieldPhiZYAM[4],NXPhi[3],NYieldPhiZYAM[5]); 
    gNPhiSpectraZYAM3=new TGraphErrors(nMainLoop,NXPhi[4],NYieldPhiZYAM[4],NXPhi[5],NYieldPhiZYAM[5]); 
    gAPhiSpectraZYAM1=new TGraphErrors(nMainLoop,AXPhi[0],AYieldPhiZYAM[4],AXPhi[1],AYieldPhiZYAM[5]); 
    gAPhiSpectraZYAM2=new TGraphErrors(nMainLoop,AXPhi[2],AYieldPhiZYAM[4],AXPhi[3],AYieldPhiZYAM[5]); 
    gAPhiSpectraZYAM3=new TGraphErrors(nMainLoop,AXPhi[4],AYieldPhiZYAM[4],AXPhi[5],AYieldPhiZYAM[5]); 

    gNPhiZt1ZYAM1=new TGraphErrors(nMainLoop,NXPhi[6],NYieldPhiZYAM[6],NXPhi[7],NYieldPhiZYAM[7]); 
    gNPhiZt1ZYAM2=new TGraphErrors(nMainLoop,NXPhi[8],NYieldPhiZYAM[6],NXPhi[9],NYieldPhiZYAM[7]); 
    gNPhiZt1ZYAM3=new TGraphErrors(nMainLoop,NXPhi[10],NYieldPhiZYAM[6],NXPhi[11],NYieldPhiZYAM[7]); 
    gAPhiZt1ZYAM1=new TGraphErrors(nMainLoop,AXPhi[6],AYieldPhiZYAM[6],AXPhi[7],AYieldPhiZYAM[7]); 
    gAPhiZt1ZYAM2=new TGraphErrors(nMainLoop,AXPhi[8],AYieldPhiZYAM[6],AXPhi[9],AYieldPhiZYAM[7]); 
    gAPhiZt1ZYAM3=new TGraphErrors(nMainLoop,AXPhi[10],AYieldPhiZYAM[6],AXPhi[11],AYieldPhiZYAM[7]); 

    gNPhiZt2ZYAM1=new TGraphErrors(nMainLoop,NXPhi[12],NYieldPhiZYAM[8],NXPhi[13],NYieldPhiZYAM[9]); 
    gNPhiZt2ZYAM2=new TGraphErrors(nMainLoop,NXPhi[14],NYieldPhiZYAM[8],NXPhi[15],NYieldPhiZYAM[9]); 
    gNPhiZt2ZYAM3=new TGraphErrors(nMainLoop,NXPhi[16],NYieldPhiZYAM[8],NXPhi[17],NYieldPhiZYAM[9]); 
    gAPhiZt2ZYAM1=new TGraphErrors(nMainLoop,AXPhi[12],AYieldPhiZYAM[8],AXPhi[13],AYieldPhiZYAM[9]); 
    gAPhiZt2ZYAM2=new TGraphErrors(nMainLoop,AXPhi[14],AYieldPhiZYAM[8],AXPhi[15],AYieldPhiZYAM[9]); 
    gAPhiZt2ZYAM3=new TGraphErrors(nMainLoop,AXPhi[16],AYieldPhiZYAM[8],AXPhi[17],AYieldPhiZYAM[9]); 
    // cout << "NXPhiFit[2] " << NXPhiFit[2][0] << " " << NXPhiFit[2][1] << endl;
    gNPhiSpectraFit1=new TGraphErrors(nMainLoop,NXPhiFit[0],NYieldPhiFit[4],NXPhiFit[1],NYieldPhiFit[5]); 
    gNPhiSpectraFit2=new TGraphErrors(nMainLoop,NXPhiFit[2],NYieldPhiFit[4],NXPhiFit[3],NYieldPhiFit[5]); 
    gNPhiSpectraFit3=new TGraphErrors(nMainLoop,NXPhiFit[4],NYieldPhiFit[4],NXPhiFit[5],NYieldPhiFit[5]); 
    gAPhiSpectraFit1=new TGraphErrors(nMainLoop,AXPhiFit[0],AYieldPhiFit[4],AXPhiFit[1],AYieldPhiFit[5]); 
    gAPhiSpectraFit2=new TGraphErrors(nMainLoop,AXPhiFit[2],AYieldPhiFit[4],AXPhiFit[3],AYieldPhiFit[5]); 
    gAPhiSpectraFit3=new TGraphErrors(nMainLoop,AXPhiFit[4],AYieldPhiFit[4],AXPhiFit[5],AYieldPhiFit[5]); 

    gNPhiZt1Fit1=new TGraphErrors(nMainLoop,NXPhiFit[6],NYieldPhiFit[6],NXPhiFit[7],NYieldPhiFit[7]); 
    gNPhiZt1Fit2=new TGraphErrors(nMainLoop,NXPhiFit[8],NYieldPhiFit[6],NXPhiFit[9],NYieldPhiFit[7]); 
    gNPhiZt1Fit3=new TGraphErrors(nMainLoop,NXPhiFit[10],NYieldPhiFit[6],NXPhiFit[11],NYieldPhiFit[7]); 
    gAPhiZt1Fit1=new TGraphErrors(nMainLoop,AXPhiFit[6],AYieldPhiFit[6],AXPhiFit[7],AYieldPhiFit[7]); 
    gAPhiZt1Fit2=new TGraphErrors(nMainLoop,AXPhiFit[8],AYieldPhiFit[6],AXPhiFit[9],AYieldPhiFit[7]); 
    gAPhiZt1Fit3=new TGraphErrors(nMainLoop,AXPhiFit[10],AYieldPhiFit[6],AXPhiFit[11],AYieldPhiFit[7]); 

    gNPhiWidth1=new TGraphErrors(nMainLoop,NXPhiFit[0],NWidthPhi[0],NXPhiFit[1],NWidthPhi[1]);
    gNPhiWidth2=new TGraphErrors(nMainLoop,NXPhiFit[2],NWidthPhi[0],NXPhiFit[3],NWidthPhi[1]);
    gNPhiWidth3=new TGraphErrors(nMainLoop,NXPhiFit[4],NWidthPhi[0],NXPhiFit[5],NWidthPhi[1]);
    gAPhiWidth1=new TGraphErrors(nMainLoop,AXPhiFit[0],AWidthPhi[0],AXPhiFit[1],AWidthPhi[1]);
    gAPhiWidth2=new TGraphErrors(nMainLoop,AXPhiFit[2],AWidthPhi[0],AXPhiFit[3],AWidthPhi[1]);
    gAPhiWidth3=new TGraphErrors(nMainLoop,AXPhiFit[4],AWidthPhi[0],AXPhiFit[5],AWidthPhi[1]);


    gNEtaSpectraZYAM1=new TGraphErrors(nMainLoop,NXEta[0],NYieldEtaZYAM[4],NXEta[1],NYieldEtaZYAM[5]); 
    gNEtaSpectraZYAM2=new TGraphErrors(nMainLoop,NXEta[2],NYieldEtaZYAM[4],NXEta[3],NYieldEtaZYAM[5]); 
    gNEtaSpectraZYAM3=new TGraphErrors(nMainLoop,NXEta[4],NYieldEtaZYAM[4],NXEta[5],NYieldEtaZYAM[5]); 
    gAEtaSpectraZYAM1=new TGraphErrors(nMainLoop,AXEta[0],AYieldEtaZYAM[4],AXEta[1],AYieldEtaZYAM[5]); 
    gAEtaSpectraZYAM2=new TGraphErrors(nMainLoop,AXEta[2],AYieldEtaZYAM[4],AXEta[3],AYieldEtaZYAM[5]); 
    gAEtaSpectraZYAM3=new TGraphErrors(nMainLoop,AXEta[4],AYieldEtaZYAM[4],AXEta[5],AYieldEtaZYAM[5]); 

    gNEtaZt1ZYAM1=new TGraphErrors(nMainLoop,NXEta[6],NYieldEtaZYAM[6],NXEta[7],NYieldEtaZYAM[7]); 
    gNEtaZt1ZYAM2=new TGraphErrors(nMainLoop,NXEta[8],NYieldEtaZYAM[6],NXEta[9],NYieldEtaZYAM[7]); 
    gNEtaZt1ZYAM3=new TGraphErrors(nMainLoop,NXEta[10],NYieldEtaZYAM[6],NXEta[11],NYieldEtaZYAM[7]); 
    gAEtaZt1ZYAM1=new TGraphErrors(nMainLoop,AXEta[6],AYieldEtaZYAM[6],AXEta[7],AYieldEtaZYAM[7]); 
    gAEtaZt1ZYAM2=new TGraphErrors(nMainLoop,AXEta[8],AYieldEtaZYAM[6],AXEta[9],AYieldEtaZYAM[7]); 
    gAEtaZt1ZYAM3=new TGraphErrors(nMainLoop,AXEta[10],AYieldEtaZYAM[6],AXEta[11],AYieldEtaZYAM[7]); 

    gNEtaZt2ZYAM1=new TGraphErrors(nMainLoop,NXEta[12],NYieldEtaZYAM[8],NXEta[13],NYieldEtaZYAM[9]); 
    gNEtaZt2ZYAM2=new TGraphErrors(nMainLoop,NXEta[14],NYieldEtaZYAM[8],NXEta[15],NYieldEtaZYAM[9]); 
    gNEtaZt2ZYAM3=new TGraphErrors(nMainLoop,NXEta[16],NYieldEtaZYAM[8],NXEta[17],NYieldEtaZYAM[9]); 
    gAEtaZt2ZYAM1=new TGraphErrors(nMainLoop,AXEta[12],AYieldEtaZYAM[8],AXEta[13],AYieldEtaZYAM[9]); 
    gAEtaZt2ZYAM2=new TGraphErrors(nMainLoop,AXEta[14],AYieldEtaZYAM[8],AXEta[15],AYieldEtaZYAM[9]); 
    gAEtaZt2ZYAM3=new TGraphErrors(nMainLoop,AXEta[16],AYieldEtaZYAM[8],AXEta[17],AYieldEtaZYAM[9]); 

    gNEtaSpectraFit1=new TGraphErrors(nMainLoop,NXEtaFit[0],NYieldEtaFit[4],NXEtaFit[1],NYieldEtaFit[5]); 
    //cout << "NXEtaFit[2] " << NXEtaFit[2][0] << " " << NXEtaFit[2][1] << endl;
    gNEtaSpectraFit2=new TGraphErrors(nMainLoop,NXEtaFit[2],NYieldEtaFit[4],NXEtaFit[3],NYieldEtaFit[5]); 
    gNEtaSpectraFit3=new TGraphErrors(nMainLoop,NXEtaFit[4],NYieldEtaFit[4],NXEtaFit[5],NYieldEtaFit[5]); 
    gAEtaSpectraFit1=new TGraphErrors(nMainLoop,AXEtaFit[0],AYieldEtaFit[4],AXEtaFit[1],AYieldEtaFit[5]); 
    gAEtaSpectraFit2=new TGraphErrors(nMainLoop,AXEtaFit[2],AYieldEtaFit[4],AXEtaFit[3],AYieldEtaFit[5]); 
    gAEtaSpectraFit3=new TGraphErrors(nMainLoop,AXEtaFit[4],AYieldEtaFit[4],AXEtaFit[5],AYieldEtaFit[5]); 

    gNEtaZt1Fit1=new TGraphErrors(nMainLoop,NXEtaFit[6],NYieldEtaFit[6],NXEtaFit[7],NYieldEtaFit[7]); 
    gNEtaZt1Fit2=new TGraphErrors(nMainLoop,NXEtaFit[8],NYieldEtaFit[6],NXEtaFit[9],NYieldEtaFit[7]); 
    gNEtaZt1Fit3=new TGraphErrors(nMainLoop,NXEtaFit[10],NYieldEtaFit[6],NXEtaFit[11],NYieldEtaFit[7]); 
    gAEtaZt1Fit1=new TGraphErrors(nMainLoop,AXEtaFit[6],AYieldEtaFit[6],AXEtaFit[7],AYieldEtaFit[7]); 
    gAEtaZt1Fit2=new TGraphErrors(nMainLoop,AXEtaFit[8],AYieldEtaFit[6],AXEtaFit[9],AYieldEtaFit[7]); 
    gAEtaZt1Fit3=new TGraphErrors(nMainLoop,AXEtaFit[10],AYieldEtaFit[6],AXEtaFit[11],AYieldEtaFit[7]); 

    gNEtaWidth1=new TGraphErrors(nMainLoop,NXEtaFit[0],NWidthEta[0],NXEtaFit[1],NWidthEta[1]);
    gNEtaWidth2=new TGraphErrors(nMainLoop,NXEtaFit[2],NWidthEta[0],NXEtaFit[3],NWidthEta[1]);
    gNEtaWidth3=new TGraphErrors(nMainLoop,NXEtaFit[4],NWidthEta[0],NXEtaFit[5],NWidthEta[1]);
  
    float plotmin=100, plotmax=-100;
    float plotminW=100,plotmaxW=-100;
    float plotminZt=100, plotmaxZt=-100;
    // cout << plotmin << endl;
    for(int i=4;i<=8;i+=2){
      for(int j=0;j<nMainLoop;j++){
	if(NYieldPhiZYAM[i][j]<plotmin&&NYieldPhiZYAM[i][j]>0)plotmin=NYieldPhiZYAM[i][j];
	if(NYieldPhiZYAM[i][j]>plotmax&&NYieldPhiZYAM[i][j]>0)plotmax=NYieldPhiZYAM[i][j];
	if(AYieldPhiZYAM[i][j]<plotmin&&AYieldPhiZYAM[i][j]>0)plotmin=AYieldPhiZYAM[i][j];
	if(AYieldPhiZYAM[i][j]>plotmax&&AYieldPhiZYAM[i][j]>0)plotmax=AYieldPhiZYAM[i][j];
	if(NYieldEtaZYAM[i][j]<plotmin&&NYieldEtaZYAM[i][j]>0)plotmin=NYieldEtaZYAM[i][j];
	if(NYieldEtaZYAM[i][j]>plotmax&&NYieldEtaZYAM[i][j]>0)plotmax=NYieldEtaZYAM[i][j];
	if(AYieldEtaZYAM[i][j]<plotmin&&AYieldEtaZYAM[i][j]>0)plotmin=AYieldEtaZYAM[i][j];
	if(AYieldEtaZYAM[i][j]>plotmax&&AYieldEtaZYAM[i][j]>0)plotmax=AYieldEtaZYAM[i][j];
	
	if(i==6){
	  if(NYieldPhiZYAM[i][j]<plotminZt&&NYieldPhiZYAM[i][j]>0)plotminZt=NYieldPhiZYAM[i][j];
	  if(NYieldPhiZYAM[i][j]>plotmaxZt&&NYieldPhiZYAM[i][j]>0)plotmaxZt=NYieldPhiZYAM[i][j];
	  if(AYieldPhiZYAM[i][j]<plotminZt&&AYieldPhiZYAM[i][j]>0)plotminZt=AYieldPhiZYAM[i][j];
	  if(AYieldPhiZYAM[i][j]>plotmaxZt&&AYieldPhiZYAM[i][j]>0)plotmaxZt=AYieldPhiZYAM[i][j];
	  /*
	  if(NYieldEtaZYAM[i][j]<plotminZt&&NYieldEtaZYAM[i][j]>0)plotminZt=NYieldEtaZYAM[i][j];
	  if(NYieldEtaZYAM[i][j]>plotmaxZt&&NYieldEtaZYAM[i][j]>0)plotmaxZt=NYieldEtaZYAM[i][j];
	  if(AYieldEtaZYAM[i][j]<plotminZt&&AYieldEtaZYAM[i][j]>0)plotminZt=AYieldEtaZYAM[i][j];
	  if(AYieldEtaZYAM[i][j]>plotmaxZt&&AYieldEtaZYAM[i][j]>0)plotmaxZt=AYieldEtaZYAM[i][j];
	  */
	}

	if(DrawFit){
	  if(NYieldPhiFit[i][j]<plotmin&&NYieldPhiFit[i][j]>0)plotmin=NYieldPhiFit[i][j];
	  if(NYieldPhiFit[i][j]>plotmax&&NYieldPhiFit[i][j]>0)plotmax=NYieldPhiFit[i][j];
	  if(AYieldPhiFit[i][j]<plotmin&&AYieldPhiFit[i][j]>0)plotmin=AYieldPhiFit[i][j];
	  if(AYieldPhiFit[i][j]>plotmax&&AYieldPhiFit[i][j]>0)plotmax=AYieldPhiFit[i][j];
	}
	if(i==4){
	  if(NWidthPhi[0][j]<plotminW)plotminW=NWidthPhi[0][j];
	  if(NWidthPhi[0][j]>plotmaxW&&NWidthPhi[0][j]<5)plotmaxW=NWidthPhi[0][j];
	  if(AWidthPhi[0][j]<plotminW)plotminW=AWidthPhi[0][j];
	  if(AWidthPhi[0][j]>plotmaxW&&AWidthPhi[0][j]<5)plotmaxW=AWidthPhi[0][j];
	  if(NWidthEta[0][j]<plotminW)plotminW=NWidthEta[0][j];
	  if(NWidthEta[0][j]>plotmaxW&&NWidthEta[0][j]<5)plotmaxW=NWidthEta[0][j];
	}
      }
    }
    plotmin/=2;
    plotmax*=100+1;
    plotminZt/=2;
    plotmaxZt*=2.2;
    //cout << plotmin << endl;
    //if(plotmin<=0)plotmin=0.0001;
    plotminW=0;
    plotmaxW*=1.1;
    plotmaxW=1;

    gNPhiSpectraZYAM1->SetMarkerStyle(MarkerNearPhi);
    gNPhiSpectraZYAM1->SetMarkerSize(MarkerSize);
    gNPhiSpectraZYAM1->SetMarkerColor(ColorNearPhi);
    gNPhiSpectraZYAM1->SetLineColor(ColorNearPhi);
    gNPhiSpectraZYAM2->SetMarkerStyle(MarkerNearPhi);
    gNPhiSpectraZYAM2->SetMarkerSize(MarkerSize);
    gNPhiSpectraZYAM2->SetMarkerColor(ColorNearPhi);
    gNPhiSpectraZYAM2->SetLineColor(ColorNearPhi);
    gNPhiSpectraZYAM3->SetMarkerStyle(MarkerNearPhi);
    gNPhiSpectraZYAM3->SetMarkerSize(MarkerSize);
    gNPhiSpectraZYAM3->SetMarkerColor(ColorNearPhi);
    gNPhiSpectraZYAM3->SetLineColor(ColorNearPhi);
    gNPhiZt1ZYAM1->SetMarkerStyle(MarkerNearPhi);
    gNPhiZt1ZYAM1->SetMarkerSize(MarkerSize);
    gNPhiZt1ZYAM1->SetLineColor(ColorNearPhi);
    gNPhiZt1ZYAM1->SetMarkerColor(ColorNearPhi);
    gNPhiZt1ZYAM2->SetMarkerStyle(MarkerNearPhi);
    gNPhiZt1ZYAM2->SetMarkerSize(MarkerSize);
    gNPhiZt1ZYAM2->SetLineColor(ColorNearPhi);
    gNPhiZt1ZYAM2->SetMarkerColor(ColorNearPhi);
    gNPhiZt1ZYAM3->SetMarkerStyle(MarkerNearPhi);
    gNPhiZt1ZYAM3->SetMarkerSize(MarkerSize);
    gNPhiZt1ZYAM3->SetLineColor(ColorNearPhi);
    gNPhiZt1ZYAM3->SetMarkerColor(ColorNearPhi);
    gNPhiZt2ZYAM3->SetMarkerStyle(MarkerNearPhi);
    gNPhiZt2ZYAM3->SetMarkerSize(MarkerSize);
    gNPhiZt2ZYAM1->SetLineColor(ColorNearPhi);
    gNPhiZt2ZYAM1->SetMarkerColor(ColorNearPhi);
    gNPhiZt2ZYAM2->SetMarkerStyle(MarkerNearPhi);
    gNPhiZt2ZYAM2->SetMarkerSize(MarkerSize);
    gNPhiZt2ZYAM2->SetLineColor(ColorNearPhi);
    gNPhiZt2ZYAM2->SetMarkerColor(ColorNearPhi);
    gNPhiZt2ZYAM3->SetMarkerStyle(MarkerNearPhi);
    gNPhiZt2ZYAM3->SetMarkerSize(MarkerSize);
    gNPhiZt2ZYAM3->SetLineColor(ColorNearPhi);
    gNPhiZt2ZYAM3->SetMarkerColor(ColorNearPhi);

    gAPhiSpectraZYAM1->SetMarkerStyle(MarkerAwayPhi);
    gAPhiSpectraZYAM1->SetMarkerSize(MarkerSize);
    gAPhiSpectraZYAM1->SetMarkerColor(ColorAwayPhi);
    gAPhiSpectraZYAM1->SetLineColor(ColorAwayPhi);
    gAPhiSpectraZYAM2->SetMarkerStyle(MarkerAwayPhi);
    gAPhiSpectraZYAM2->SetMarkerSize(MarkerSize);
    gAPhiSpectraZYAM2->SetMarkerColor(ColorAwayPhi);
    gAPhiSpectraZYAM2->SetLineColor(ColorAwayPhi);
    gAPhiSpectraZYAM3->SetMarkerStyle(MarkerAwayPhi);
    gAPhiSpectraZYAM3->SetMarkerSize(MarkerSize);
    gAPhiSpectraZYAM3->SetMarkerColor(ColorAwayPhi);
    gAPhiSpectraZYAM3->SetLineColor(ColorAwayPhi);
    gAPhiZt1ZYAM1->SetMarkerStyle(MarkerAwayPhi);
    gAPhiZt1ZYAM1->SetMarkerSize(MarkerSize);
    gAPhiZt1ZYAM1->SetLineColor(ColorAwayPhi);
    gAPhiZt1ZYAM1->SetMarkerColor(ColorAwayPhi);
    gAPhiZt1ZYAM2->SetMarkerStyle(MarkerAwayPhi);
    gAPhiZt1ZYAM2->SetMarkerSize(MarkerSize);
    gAPhiZt1ZYAM2->SetLineColor(ColorAwayPhi);
    gAPhiZt1ZYAM2->SetMarkerColor(ColorAwayPhi);
    gAPhiZt1ZYAM3->SetMarkerStyle(MarkerAwayPhi);
    gAPhiZt1ZYAM3->SetMarkerSize(MarkerSize);
    gAPhiZt1ZYAM3->SetLineColor(ColorAwayPhi);
    gAPhiZt1ZYAM3->SetMarkerColor(ColorAwayPhi);
    gAPhiZt2ZYAM3->SetMarkerStyle(MarkerAwayPhi);
    gAPhiZt2ZYAM3->SetMarkerSize(MarkerSize);
    gAPhiZt2ZYAM1->SetLineColor(ColorAwayPhi);
    gAPhiZt2ZYAM1->SetMarkerColor(ColorAwayPhi);
    gAPhiZt2ZYAM2->SetMarkerStyle(MarkerAwayPhi);
    gAPhiZt2ZYAM2->SetMarkerSize(MarkerSize);
    gAPhiZt2ZYAM2->SetLineColor(ColorAwayPhi);
    gAPhiZt2ZYAM2->SetMarkerColor(ColorAwayPhi);
    gAPhiZt2ZYAM3->SetMarkerStyle(MarkerAwayPhi);
    gAPhiZt2ZYAM3->SetMarkerSize(MarkerSize);
    gAPhiZt2ZYAM3->SetLineColor(ColorAwayPhi);
    gAPhiZt2ZYAM3->SetMarkerColor(ColorAwayPhi);

    gNPhiSpectraFit1->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiSpectraFit1->SetMarkerSize(MarkerSize);
    gNPhiSpectraFit1->SetMarkerColor(ColorNearPhiFit);
    gNPhiSpectraFit1->SetLineColor(ColorNearPhiFit);
    gNPhiSpectraFit2->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiSpectraFit2->SetMarkerSize(MarkerSize);
    gNPhiSpectraFit2->SetMarkerColor(ColorNearPhiFit);
    gNPhiSpectraFit2->SetLineColor(ColorNearPhiFit);
    gNPhiSpectraFit3->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiSpectraFit3->SetMarkerSize(MarkerSize);
    gNPhiSpectraFit3->SetMarkerColor(ColorNearPhiFit);
    gNPhiSpectraFit3->SetLineColor(ColorNearPhiFit);
    gNPhiZt1Fit1->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiZt1Fit1->SetMarkerSize(MarkerSize);
    gNPhiZt1Fit1->SetLineColor(ColorNearPhiFit);
    gNPhiZt1Fit1->SetMarkerColor(ColorNearPhiFit);
    gNPhiZt1Fit2->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiZt1Fit2->SetMarkerSize(MarkerSize);
    gNPhiZt1Fit2->SetLineColor(ColorNearPhiFit);
    gNPhiZt1Fit2->SetMarkerColor(ColorNearPhiFit);
    gNPhiZt1Fit3->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiZt1Fit3->SetMarkerSize(MarkerSize);
    gNPhiZt1Fit3->SetLineColor(ColorNearPhiFit);
    gNPhiZt1Fit3->SetMarkerColor(ColorNearPhiFit);
    /*
      gNPhiZt2Fit3->SetMarkerStyle(24);
      gNPhiZt2Fit1->SetLineColor(kBlue+3);
      gNPhiZt2Fit1->SetMarkerColor(kBlue+3);
      gNPhiZt2Fit2->SetMarkerStyle(24);
      gNPhiZt2Fit2->SetLineColor(kBlue+3);
      gNPhiZt2Fit2->SetMarkerColor(kBlue+3);
      gNPhiZt2Fit3->SetMarkerStyle(24);
      gNPhiZt2Fit3->SetLineColor(kBlue+3);
      gNPhiZt2Fit3->SetMarkerColor(kBlue+3);
    */
    gAPhiSpectraFit1->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiSpectraFit1->SetMarkerSize(MarkerSize);
    gAPhiSpectraFit1->SetMarkerColor(ColorAwayPhiFit);
    gAPhiSpectraFit1->SetLineColor(ColorAwayPhiFit);
    gAPhiSpectraFit2->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiSpectraFit2->SetMarkerSize(MarkerSize);
    gAPhiSpectraFit2->SetMarkerColor(ColorAwayPhiFit);
    gAPhiSpectraFit2->SetLineColor(ColorAwayPhiFit);
    gAPhiSpectraFit3->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiSpectraFit3->SetMarkerSize(MarkerSize);
    gAPhiSpectraFit3->SetMarkerColor(ColorAwayPhiFit);
    gAPhiSpectraFit3->SetLineColor(ColorAwayPhiFit);
    gAPhiZt1Fit1->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiZt1Fit1->SetMarkerSize(MarkerSize);
    gAPhiZt1Fit1->SetLineColor(ColorAwayPhiFit);
    gAPhiZt1Fit1->SetMarkerColor(ColorAwayPhiFit);
    gAPhiZt1Fit2->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiZt1Fit2->SetMarkerSize(MarkerSize);
    gAPhiZt1Fit2->SetLineColor(ColorAwayPhiFit);
    gAPhiZt1Fit2->SetMarkerColor(ColorAwayPhiFit);
    gAPhiZt1Fit3->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiZt1Fit3->SetMarkerSize(MarkerSize);
    gAPhiZt1Fit3->SetLineColor(ColorAwayPhiFit);
    gAPhiZt1Fit3->SetMarkerColor(ColorAwayPhiFit);
    /*
      gAPhiZt2Fit3->SetMarkerStyle(25);
      gAPhiZt2Fit1->SetLineColor(kRed+2);
      gAPhiZt2Fit1->SetMarkerColor(kRed+2);
      gAPhiZt2Fit2->SetMarkerStyle(25);
      gAPhiZt2Fit2->SetLineColor(kRed+2);
      gAPhiZt2Fit2->SetMarkerColor(kRed+2);
      gAPhiZt2Fit3->SetMarkerStyle(25);
      gAPhiZt2Fit3->SetLineColor(kRed+2);
      gAPhiZt2Fit3->SetMarkerColor(kRed+2);
    */
  
    gNPhiWidth1->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiWidth1->SetMarkerSize(MarkerSize);
    gNPhiWidth1->SetMarkerColor(ColorNearPhiFit);
    gNPhiWidth1->SetLineColor(ColorNearPhiFit);
    gNPhiWidth2->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiWidth2->SetMarkerSize(MarkerSize);
    gNPhiWidth2->SetMarkerColor(ColorNearPhiFit);
    gNPhiWidth2->SetLineColor(ColorNearPhiFit);
    gNPhiWidth3->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiWidth3->SetMarkerSize(MarkerSize);
    gNPhiWidth3->SetMarkerColor(ColorNearPhiFit);
    gNPhiWidth3->SetLineColor(ColorNearPhiFit);

    gAPhiWidth1->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiWidth1->SetMarkerSize(MarkerSize);
    gAPhiWidth1->SetMarkerColor(ColorAwayPhiFit);
    gAPhiWidth1->SetLineColor(ColorAwayPhiFit);
    gAPhiWidth2->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiWidth2->SetMarkerSize(MarkerSize);
    gAPhiWidth2->SetMarkerColor(ColorAwayPhiFit);
    gAPhiWidth2->SetLineColor(ColorAwayPhiFit);
    gAPhiWidth3->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiWidth3->SetMarkerSize(MarkerSize);
    gAPhiWidth3->SetMarkerColor(ColorAwayPhiFit);
    gAPhiWidth3->SetLineColor(ColorAwayPhiFit);


    gNEtaSpectraZYAM1->SetMarkerStyle(MarkerNearEta);
    gNEtaSpectraZYAM1->SetMarkerSize(MarkerSize);
    gNEtaSpectraZYAM1->SetMarkerColor(ColorNearEta);
    gNEtaSpectraZYAM1->SetLineColor(ColorNearEta);
    gNEtaSpectraZYAM2->SetMarkerStyle(MarkerNearEta);
    gNEtaSpectraZYAM2->SetMarkerSize(MarkerSize);
    gNEtaSpectraZYAM2->SetMarkerColor(ColorNearEta);
    gNEtaSpectraZYAM2->SetLineColor(ColorNearEta);
    gNEtaSpectraZYAM3->SetMarkerStyle(MarkerNearEta);
    gNEtaSpectraZYAM3->SetMarkerSize(MarkerSize);
    gNEtaSpectraZYAM3->SetMarkerColor(ColorNearEta);
    gNEtaSpectraZYAM3->SetLineColor(ColorNearEta);
    gNEtaZt1ZYAM1->SetMarkerStyle(MarkerNearEta);
    gNEtaZt1ZYAM1->SetMarkerSize(MarkerSize);
    gNEtaZt1ZYAM1->SetLineColor(ColorNearEta);
    gNEtaZt1ZYAM1->SetMarkerColor(ColorNearEta);
    gNEtaZt1ZYAM2->SetMarkerStyle(MarkerNearEta);
    gNEtaZt1ZYAM2->SetMarkerSize(MarkerSize);
    gNEtaZt1ZYAM2->SetLineColor(ColorNearEta);
    gNEtaZt1ZYAM2->SetMarkerColor(ColorNearEta);
    gNEtaZt1ZYAM3->SetMarkerStyle(MarkerNearEta);
    gNEtaZt1ZYAM3->SetMarkerSize(MarkerSize);
    gNEtaZt1ZYAM3->SetLineColor(ColorNearEta);
    gNEtaZt1ZYAM3->SetMarkerColor(ColorNearEta);
    gNEtaZt2ZYAM3->SetMarkerStyle(MarkerNearEta);
    gNEtaZt2ZYAM3->SetMarkerSize(MarkerSize);
    gNEtaZt2ZYAM1->SetLineColor(ColorNearEta);
    gNEtaZt2ZYAM1->SetMarkerColor(ColorNearEta);
    gNEtaZt2ZYAM2->SetMarkerStyle(MarkerNearEta);
    gNEtaZt2ZYAM2->SetMarkerSize(MarkerSize);
    gNEtaZt2ZYAM2->SetLineColor(ColorNearEta);
    gNEtaZt2ZYAM2->SetMarkerColor(ColorNearEta);
    gNEtaZt2ZYAM3->SetMarkerStyle(MarkerNearEta);
    gNEtaZt2ZYAM3->SetMarkerSize(MarkerSize);
    gNEtaZt2ZYAM3->SetLineColor(ColorNearEta);
    gNEtaZt2ZYAM3->SetMarkerColor(ColorNearEta);

    gAEtaSpectraZYAM1->SetMarkerStyle(MarkerAwayEta);
    gAEtaSpectraZYAM1->SetMarkerSize(MarkerSize);
    gAEtaSpectraZYAM1->SetMarkerColor(ColorAwayEta);
    gAEtaSpectraZYAM1->SetLineColor(ColorAwayEta);
    gAEtaSpectraZYAM2->SetMarkerStyle(MarkerAwayEta);
    gAEtaSpectraZYAM2->SetMarkerSize(MarkerSize);
    gAEtaSpectraZYAM2->SetMarkerColor(ColorAwayEta);
    gAEtaSpectraZYAM2->SetLineColor(ColorAwayEta);
    gAEtaSpectraZYAM3->SetMarkerStyle(MarkerAwayEta);
    gAEtaSpectraZYAM3->SetMarkerSize(MarkerSize);
    gAEtaSpectraZYAM3->SetMarkerColor(ColorAwayEta);
    gAEtaSpectraZYAM3->SetLineColor(ColorAwayEta);
    gAEtaZt1ZYAM1->SetMarkerStyle(MarkerAwayEta);
    gAEtaZt1ZYAM1->SetMarkerSize(MarkerSize);
    gAEtaZt1ZYAM1->SetLineColor(ColorAwayEta);
    gAEtaZt1ZYAM1->SetMarkerColor(ColorAwayEta);
    gAEtaZt1ZYAM2->SetMarkerStyle(MarkerAwayEta);
    gAEtaZt1ZYAM2->SetMarkerSize(MarkerSize);
    gAEtaZt1ZYAM2->SetLineColor(ColorAwayEta);
    gAEtaZt1ZYAM2->SetMarkerColor(ColorAwayEta);
    gAEtaZt1ZYAM3->SetMarkerStyle(MarkerAwayEta);
    gAEtaZt1ZYAM3->SetMarkerSize(MarkerSize);
    gAEtaZt1ZYAM3->SetLineColor(ColorAwayEta);
    gAEtaZt1ZYAM3->SetMarkerColor(ColorAwayEta);
    gAEtaZt2ZYAM3->SetMarkerStyle(MarkerAwayEta);
    gAEtaZt2ZYAM3->SetMarkerSize(MarkerSize);
    gAEtaZt2ZYAM1->SetLineColor(ColorAwayEta);
    gAEtaZt2ZYAM1->SetMarkerColor(ColorAwayEta);
    gAEtaZt2ZYAM2->SetMarkerStyle(MarkerAwayEta);
    gAEtaZt2ZYAM2->SetMarkerSize(MarkerSize);
    gAEtaZt2ZYAM2->SetLineColor(ColorAwayEta);
    gAEtaZt2ZYAM2->SetMarkerColor(ColorAwayEta);
    gAEtaZt2ZYAM3->SetMarkerStyle(MarkerAwayEta);
    gAEtaZt2ZYAM3->SetMarkerSize(MarkerSize);
    gAEtaZt2ZYAM3->SetLineColor(ColorAwayEta);
    gAEtaZt2ZYAM3->SetMarkerColor(ColorAwayEta);

    gNEtaSpectraFit1->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaSpectraFit1->SetMarkerSize(MarkerSize);
    gNEtaSpectraFit1->SetMarkerColor(ColorNearEtaFit);
    gNEtaSpectraFit1->SetLineColor(ColorNearEtaFit);
    gNEtaSpectraFit2->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaSpectraFit2->SetMarkerSize(MarkerSize);
    gNEtaSpectraFit2->SetMarkerColor(ColorNearEtaFit);
    gNEtaSpectraFit2->SetLineColor(ColorNearEtaFit);
    gNEtaSpectraFit3->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaSpectraFit3->SetMarkerSize(MarkerSize);
    gNEtaSpectraFit3->SetMarkerColor(ColorNearEtaFit);
    gNEtaSpectraFit3->SetLineColor(ColorNearEtaFit);
    gNEtaZt1Fit1->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaZt1Fit1->SetMarkerSize(MarkerSize);
    gNEtaZt1Fit1->SetLineColor(ColorNearEtaFit);
    gNEtaZt1Fit1->SetMarkerColor(ColorNearEtaFit);
    gNEtaZt1Fit2->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaZt1Fit2->SetMarkerSize(MarkerSize);
    gNEtaZt1Fit2->SetLineColor(ColorNearEtaFit);
    gNEtaZt1Fit2->SetMarkerColor(ColorNearEtaFit);
    gNEtaZt1Fit3->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaZt1Fit3->SetMarkerSize(MarkerSize);
    gNEtaZt1Fit3->SetLineColor(ColorNearEtaFit);
    gNEtaZt1Fit3->SetMarkerColor(ColorNearEtaFit);

    gNEtaWidth1->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaWidth1->SetMarkerSize(MarkerSize);
    gNEtaWidth1->SetMarkerColor(ColorNearEtaFit);
    gNEtaWidth1->SetLineColor(ColorNearEtaFit);
    gNEtaWidth2->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaWidth2->SetMarkerSize(MarkerSize);
    gNEtaWidth2->SetMarkerColor(ColorNearEtaFit);
    gNEtaWidth2->SetLineColor(ColorNearEtaFit);
    gNEtaWidth3->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaWidth3->SetMarkerSize(MarkerSize);
    gNEtaWidth3->SetMarkerColor(ColorNearEtaFit);
    gNEtaWidth3->SetLineColor(ColorNearEtaFit);

  
    sprintf(outName,"%s/DrawSpectra_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_%d%s.root",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,APtTPtMult,FitTit[DrawFit]);
    TFile *fout=new TFile(name,"recreate");


    cSpectraZYAM1=new TCanvas("cSpectraZYAM1","Spectra MidPoint",800,600);
    SetMargins1D(cSpectraZYAM1);
    sprintf(name,"Associated Particle Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hSpectraZYAM1=new TH2F("hSpectraZYAM1",name,10,0,TPt1,10,plotmin,plotmax);
    if(NoTitle)hSpectraZYAM1->SetTitle("");
    hSpectraZYAM1->GetXaxis()->SetTitle("p_{T}^{Associated} (GeV/c) ");
    hSpectraZYAM1->GetXaxis()->SetTitleColor(1);
    hSpectraZYAM1->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dp_{T}}    ");
    SetTitles1D(hSpectraZYAM1);
    hSpectraZYAM1->Draw();
    cSpectraZYAM1->SetLogy(1);
    gNPhiSpectraZYAM1->Draw("p");
    gAPhiSpectraZYAM1->Draw("p");
    if(DrawFit)gNPhiSpectraFit1->Draw("p");
    if(DrawFit)gAPhiSpectraFit1->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_SpectraMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cSpectraZYAM1->SaveAs(outName);

    cSpectraZYAM2=new TCanvas("cSpectraZYAM2","Spectra Average",800,600);
    SetMargins1D(cSpectraZYAM2);
    sprintf(name,"Associated Particle Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hSpectraZYAM2=new TH2F("hSpectraZYAM2",name,10,0,TPt1,10,plotmin,plotmax);
    if(NoTitle)hSpectraZYAM2->SetTitle("");
    hSpectraZYAM2->GetXaxis()->SetTitle("p_{T}^{Associated} (GeV/c) ");
    hSpectraZYAM2->GetXaxis()->SetTitleColor(1);
    hSpectraZYAM2->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dp_{T}}    ");
    SetTitles1D(hSpectraZYAM2);
    hSpectraZYAM2->Draw();
    cSpectraZYAM2->SetLogy(1);
    gNPhiSpectraZYAM2->Draw("p");
    gAPhiSpectraZYAM2->Draw("p");
    if(DrawFit)gNPhiSpectraFit2->Draw("p");
    if(DrawFit)gAPhiSpectraFit2->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_SpectraAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cSpectraZYAM2->SaveAs(outName);

    cSpectraZYAM3=new TCanvas("cSpectraZYAM3","Spectra Lowest Error",800,600);
    SetMargins1D(cSpectraZYAM3);
    sprintf(name,"Associated Particle Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hSpectraZYAM3=new TH2F("hSpectraZYAM3",name,10,0,TPt1,10,plotmin,plotmax);
    if(NoTitle)hSpectraZYAM3->SetTitle("");
    hSpectraZYAM3->GetXaxis()->SetTitle("p_{T}^{Associated} (GeV/c) ");
    hSpectraZYAM3->GetXaxis()->SetTitleColor(1);
    hSpectraZYAM3->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dp_{T}}    ");
    SetTitles1D(hSpectraZYAM3);
    hSpectraZYAM3->Draw();
    cSpectraZYAM3->SetLogy(1);
    gNPhiSpectraZYAM3->Draw("p");
    gAPhiSpectraZYAM3->Draw("p");
    if(DrawFit)gNPhiSpectraFit3->Draw("p");
    if(DrawFit)gAPhiSpectraFit3->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_SpectraLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cSpectraZYAM3->SaveAs(outName);

    cZt1ZYAM1=new TCanvas("cZt1ZYAM1","Zt Spectra MidPoint",800,600);
    SetMargins1D(cZt1ZYAM1);
    sprintf(name,"Associated Particle Z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt1ZYAM1=new TH2F("hZt1ZYAM1",name,10,0,1,10,plotminZt,plotmaxZt);
    if(NoTitle)hZt1ZYAM1->SetTitle("");
    hZt1ZYAM1->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig}  ");
    hZt1ZYAM1->GetXaxis()->SetTitleColor(1);
    hZt1ZYAM1->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt1ZYAM1);
    hZt1ZYAM1->Draw();
    cZt1ZYAM1->SetLogy(1);
    gNPhiZt1ZYAM1->Draw("p");
    gAPhiZt1ZYAM1->Draw("p");
    if(DrawFit)gNPhiZt1Fit1->Draw("p");
    if(DrawFit)gAPhiZt1Fit1->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cZt1ZYAM1->SaveAs(outName);

    cZt1ZYAM2=new TCanvas("cZt1ZYAM2","Zt Spectra Average",800,600);
    SetMargins1D(cZt1ZYAM2);
    sprintf(name,"Associated Particle Z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt1ZYAM2=new TH2F("hZt1ZYAM2",name,10,0,1,10,plotminZt,plotmaxZt);
    if(NoTitle)hZt1ZYAM2->SetTitle("");
    hZt1ZYAM2->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig}  ");
    hZt1ZYAM2->GetXaxis()->SetTitleColor(1);
    hZt1ZYAM2->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt1ZYAM2);
    hZt1ZYAM2->Draw();
    cZt1ZYAM2->SetLogy(1);
    gNPhiZt1ZYAM2->Draw("p");
    gAPhiZt1ZYAM2->Draw("p");
    if(DrawFit)gNPhiZt1Fit2->Draw("p");
    if(DrawFit)gAPhiZt1Fit2->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cZt1ZYAM2->SaveAs(outName);

    cZt1ZYAM3=new TCanvas("cZt1ZYAM3","Zt Spectra Lowest Error",800,600);
    SetMargins1D(cZt1ZYAM3);
    sprintf(name,"Associated Particle Z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt1ZYAM3=new TH2F("hZt1ZYAM3",name,10,0,1,10,plotminZt,plotmaxZt);
    if(NoTitle)hZt1ZYAM3->SetTitle("");
    hZt1ZYAM3->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig} ");
    hZt1ZYAM3->GetXaxis()->SetTitleColor(1);
    hZt1ZYAM3->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt1ZYAM3);
    hZt1ZYAM3->Draw();
    cZt1ZYAM3->SetLogy(1);
    gNPhiZt1ZYAM3->Draw("p");
    gAPhiZt1ZYAM3->Draw("p");
    if(DrawFit)gNPhiZt1Fit3->Draw("p");
    if(DrawFit)gAPhiZt1Fit3->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cZt1ZYAM3->SaveAs(outName);
  
    /*
      cZt2ZYAM1=new TCanvas("cZt2ZYAM1","Spectra MidPoint",800,600);
      SetMargins1D(cZt2ZYAM1);
      sprintf(name,"Associated Particle Z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
      if(NoTitle)name="";
      TH2F *hZt2ZYAM1=new TH2F("hZt2ZYAM1",name,10,0,1,10,plotmin,plotmax);
      hZt2ZYAM1->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Near}  ");
      hZt2ZYAM1->GetXaxis()->SetTitleColor(1);
      hZt2ZYAM1->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
      hZt2ZYAM1->Draw();
      cZt2ZYAM1->SetLogy(1);
      gNPhiZt2ZYAM1->Draw("p");
      gAPhiZt2ZYAM1->Draw("p");
      gNPhiZt2Fit1->Draw("p");
      gAPhiZt2Fit1->Draw("p");
      sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
      keySymbol(.65,.85,name,1,20,0.06,1.5);
      keySymbol(.65,.8,"Away",2,21,0.06,1.5);
      keySymbol(.65,.75,"Near Fit",kBlue+3,24,0.06,1.5);
      keySymbol(.65,.7,"Away Fit",kRed+2,25,0.06,1.5);
      sprintf(outName,"%s/DrawSpectra_ZtNearMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,APtTPtMult);
      cZt2ZYAM1->SaveAs(outName);

      cZt2ZYAM2=new TCanvas("cZt2ZYAM2","Spectra Average",800,600);
      SetMargins1D(cZt2ZYAM2);
      sprintf(name,"Associated Particle Z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
      if(NoTitle)name="";
      TH2F *hZt2ZYAM2=new TH2F("hZt2ZYAM2",name,10,0,1,10,plotmin,plotmax);
      hZt2ZYAM2->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Near} (GeV/c) ");
      hZt2ZYAM2->GetXaxis()->SetTitleColor(1);
      hZt2ZYAM2->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
      hZt2ZYAM2->Draw();
      cZt2ZYAM2->SetLogy(1);
      gNPhiZt2ZYAM2->Draw("p");
      gAPhiZt2ZYAM2->Draw("p");
      gNPhiZt2Fit2->Draw("p");
      gAPhiZt2Fit2->Draw("p");
      sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
      keySymbol(.65,.85,name,1,20,0.06,1.5);
      keySymbol(.65,.8,"Away",2,21,0.06,1.5);
      keySymbol(.65,.75,"Near Fit",kBlue+3,24,0.06,1.5);
      keySymbol(.65,.7,"Away Fit",kRed+2,25,0.06,1.5);
      sprintf(outName,"%s/DrawSpectra_ZtNearAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,APtTPtMult);
      cZt2ZYAM2->SaveAs(outName);

      cZt2ZYAM3=new TCanvas("cZt2ZYAM3","Spectra Lowest Error",800,600);
      SetMargins1D(cZt2ZYAM3);
      sprintf(name,"Associated Particle Z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
      if(NoTitle)name="";
      TH2F *hZt2ZYAM3=new TH2F("hZt2ZYAM3",name,10,0,1,10,plotmin,plotmax);
      hZt2ZYAM3->GetXaxis()->SetTitle("p_{T}^{Assoc}/p_{T}^{Near} (GeV/c) ");
      hZt2ZYAM3->GetXaxis()->SetTitleColor(1);
      hZt2ZYAM3->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
      hZt2ZYAM3->Draw();
      cZt2ZYAM3->SetLogy(1);
      gNPhiZt2ZYAM3->Draw("p");
      gAPhiZt2ZYAM3->Draw("p");
      gNPhiZt2Fit3->Draw("p");
      gAPhiZt2Fit3->Draw("p");
      sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
      keySymbol(.65,.85,name,1,20,0.06,1.5);
      keySymbol(.65,.8,"Away",2,21,0.06,1.5);
      keySymbol(.65,.75,"Near Fit",kBlue+3,24,0.06,1.5);
      keySymbol(.65,.7,"Away Fit",kRed+2,25,0.06,1.5);
      sprintf(outName,"%s/DrawSpectra_ZtNearLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,APtTPtMult);
      cZt2ZYAM3->SaveAs(outName);
  
    */

    cSpectraZYAM1E=new TCanvas("cSpectraZYAM1E","Spectra MidPoint",800,600);
    SetMargins1D(cSpectraZYAM1E);
    sprintf(name,"Associated Particle Spectra Eta %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hSpectraZYAM1E=new TH2F("hSpectraZYAM1E",name,10,0,TPt1,10,plotmin,plotmax);
    if(NoTitle)hSpectraZYAM1E->SetTitle("");
    hSpectraZYAM1E->GetXaxis()->SetTitle("p_{T}^{Associated} (GeV/c) ");
    hSpectraZYAM1E->GetXaxis()->SetTitleColor(1);
    hSpectraZYAM1E->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dp_{T}}    ");
    SetTitles1D(hSpectraZYAM1E);
    hSpectraZYAM1E->Draw();
    cSpectraZYAM1E->SetLogy(1);
    gNEtaSpectraZYAM1->Draw("p");
    gAEtaSpectraZYAM1->Draw("p");
    if(DrawFit)gNEtaSpectraFit1->Draw("p");
    // if(DrawFit)gAPhiSpectraFit1->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_SpectraEtaMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cSpectraZYAM1E->SaveAs(outName);

    cSpectraZYAM2E=new TCanvas("cSpectraZYAM2E","Spectra Average",800,600);
    SetMargins1D(cSpectraZYAM2E);
    sprintf(name,"Associated Particle Spectra Eta %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hSpectraZYAM2E=new TH2F("hSpectraZYAM2E",name,10,0,TPt1,10,plotmin,plotmax);
    if(NoTitle)hSpectraZYAM2E->SetTitle("");
    hSpectraZYAM2E->GetXaxis()->SetTitle("p_{T}^{Associated} (GeV/c) ");
    hSpectraZYAM2E->GetXaxis()->SetTitleColor(1);
    hSpectraZYAM2E->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dp_{T}}    ");
    SetTitles1D(hSpectraZYAM2E);
    hSpectraZYAM2E->Draw();
    cSpectraZYAM2E->SetLogy(1);
    gNEtaSpectraZYAM2->Draw("p");
    gAEtaSpectraZYAM2->Draw("p");
    if(DrawFit)gNEtaSpectraFit2->Draw("p");
    //if(DrawFit)gAEtaSpectraFit2->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthEta);
    keySymbol(.65,.85,name,ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayEtaFit,MarkerAwayEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_SpectraEtaAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d,%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cSpectraZYAM2E->SaveAs(outName);

    cSpectraZYAM3E=new TCanvas("cSpectraZYAM3E","Spectra Lowest Error",800,600);
    SetMargins1D(cSpectraZYAM3);
    sprintf(name,"Associated Particle Spectra Eta %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hSpectraZYAM3E=new TH2F("hSpectraZYAM3E",name,10,0,TPt1,10,plotmin,plotmax);
    if(NoTitle)hSpectraZYAM3E->SetTitle("");
    hSpectraZYAM3E->GetXaxis()->SetTitle("p_{T}^{Associated} (GeV/c) ");
    hSpectraZYAM3E->GetXaxis()->SetTitleColor(1);
    hSpectraZYAM3E->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dp_{T}}    ");
    SetTitles1D(hSpectraZYAM3E);
    hSpectraZYAM3E->Draw();
    cSpectraZYAM3E->SetLogy(1);
    gNEtaSpectraZYAM3->Draw("p");
    gAEtaSpectraZYAM3->Draw("p");
    if(DrawFit)gNEtaSpectraFit3->Draw("p");
    //if(DrawFit)gAEtaSpectraFit3->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthEta);
    keySymbol(.65,.85,name,ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    //if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayEtaFit,MarkerAwayEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_SpectraEtaLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cSpectraZYAM3E->SaveAs(outName);

    cZt1ZYAM1E=new TCanvas("cZt1ZYAM1E","Spectra MidPoint",800,600);
    SetMargins1D(cZt1ZYAM1E);
    sprintf(name,"Associated Particle Z_{T} Spectra Eta %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt1ZYAM1E=new TH2F("hZt1ZYAM1E",name,10,0,1,10,plotmin,plotmax);
    if(NoTitle)hZt1ZYAM1E->SetTitle("");
    hZt1ZYAM1E->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig}  ");
    hZt1ZYAM1E->GetXaxis()->SetTitleColor(1);
    hZt1ZYAM1E->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt1ZYAM1E);
    hZt1ZYAM1E->Draw();
    cZt1ZYAM1E->SetLogy(1);
    gNEtaZt1ZYAM1->Draw("p");
    gAEtaZt1ZYAM1->Draw("p");
    if(DrawFit)gNEtaZt1Fit1->Draw("p");
    //if(DrawFit)gAEtaZt1Fit1->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthEta);
    keySymbol(.65,.85,name,ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    //if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayEtaFit,MarkerAwayEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigEtaMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cZt1ZYAM1E->SaveAs(outName);

    cZt1ZYAM2E=new TCanvas("cZt1ZYAM2E","Spectra Average",800,600);
    SetMargins1D(cZt1ZYAM2E);
    sprintf(name,"Associated Particle Z_{T} Spectra Eta %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt1ZYAM2E=new TH2F("hZt1ZYAM2E",name,10,0,1,10,plotmin,plotmax);
    if(NoTitle)hZt1ZYAM2E->SetTitle("");
    hZt1ZYAM2E->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig}  ");
    hZt1ZYAM2E->GetXaxis()->SetTitleColor(1);
    hZt1ZYAM2E->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt1ZYAM2E);
    hZt1ZYAM2E->Draw();
    cZt1ZYAM2E->SetLogy(1);
    gNEtaZt1ZYAM2->Draw("p");
    gAEtaZt1ZYAM2->Draw("p");
    if(DrawFit)gNEtaZt1Fit2->Draw("p");
    //if(DrawFit)gAEtaZt1Fit2->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthEta);
    keySymbol(.65,.85,name,ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    //if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayEtaFit,MarkerAwayEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigEtaAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cZt1ZYAM2E->SaveAs(outName);

    cZt1ZYAM3E=new TCanvas("cZt1ZYAM3E","Spectra Lowest Error",800,600);
    SetMargins1D(cZt1ZYAM3E);
    sprintf(name,"Associated Particle Z_{T} Spectra Eta %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt1ZYAM3E=new TH2F("hZt1ZYAM3E",name,10,0,1,10,plotmin,plotmax);
    if(NoTitle)hZt1ZYAM3->SetTitle("");
    hZt1ZYAM3E->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig} ");
    hZt1ZYAM3E->GetXaxis()->SetTitleColor(1);
    hZt1ZYAM3E->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt1ZYAM3E);
    hZt1ZYAM3E->Draw();
    cZt1ZYAM3E->SetLogy(1);
    gNEtaZt1ZYAM3->Draw("p");
    gAEtaZt1ZYAM3->Draw("p");
    if(DrawFit)gNEtaZt1Fit3->Draw("p");
    //if(DrawFit)gAEtaZt1Fit3->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthEta);
    keySymbol(.65,.85,name,ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    //if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayEtaFit,MarkerAwayEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigEtaLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cZt1ZYAM3E->SaveAs(outName);


    cSpectraZYAM1C=new TCanvas("cSpectraZYAM1C","Spectra MidPoint",800,600);
    SetMargins1D(cSpectraZYAM1C);
    sprintf(name,"Associated Particle Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hSpectraZYAM1C=new TH2F("hSpectraZYAM1C",name,10,0,TPt1,10,plotmin,plotmax);
    if(NoTitle)hSpectraZYAM1C->SetTitle("");
    hSpectraZYAM1C->GetXaxis()->SetTitle("p_{T}^{Associated} (GeV/c) ");
    hSpectraZYAM1C->GetXaxis()->SetTitleColor(1);
    hSpectraZYAM1C->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dp_{T}}    ");
    SetTitles1D(hSpectraZYAM1C);
    hSpectraZYAM1C->Draw();
    cSpectraZYAM1C->SetLogy(1);
    gNPhiSpectraZYAM1->Draw("p");
    gAPhiSpectraZYAM1->Draw("p");
    gNEtaSpectraZYAM1->Draw("p");
    gAEtaSpectraZYAM1->Draw("p");
    if(DrawFit)gNPhiSpectraFit1->Draw("p");
    if(DrawFit)gAPhiSpectraFit1->Draw("p");
    if(DrawFit)gNEtaSpectraFit1->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.65,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.6,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.7,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.7,"Away #Delta#eta Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_SpectraCombMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cSpectraZYAM1C->SaveAs(outName);

    cSpectraZYAM2C=new TCanvas("cSpectraZYAM2C","Spectra MidPoint",800,600);
    SetMargins1D(cSpectraZYAM2C);
    sprintf(name,"Associated Particle Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hSpectraZYAM2C=new TH2F("hSpectraZYAM2C",name,10,0,TPt1,10,plotmin,plotmax);
    if(NoTitle)hSpectraZYAM2C->SetTitle("");
    hSpectraZYAM2C->GetXaxis()->SetTitle("p_{T}^{Associated} (GeV/c) ");
    hSpectraZYAM2C->GetXaxis()->SetTitleColor(1);
    hSpectraZYAM2C->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dp_{T}}    ");

    SetTitles1D(hSpectraZYAM2C);
    hSpectraZYAM2C->Draw();
    cSpectraZYAM2C->SetLogy(1);
    gNPhiSpectraZYAM2->Draw("p");
    gAPhiSpectraZYAM2->Draw("p");
    gNEtaSpectraZYAM2->Draw("p");
    gAEtaSpectraZYAM2->Draw("p");
    if(DrawFit)gNPhiSpectraFit2->Draw("p");
    if(DrawFit)gAPhiSpectraFit2->Draw("p");
    if(DrawFit)gNEtaSpectraFit2->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.65,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.6,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.7,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.7,"Away #Delta#eta Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_SpectraCombAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cSpectraZYAM2C->SaveAs(outName);

    cSpectraZYAM3C=new TCanvas("cSpectraZYAM3C","Spectra MidPoint",800,600);
    SetMargins1D(cSpectraZYAM3C);
    sprintf(name,"Associated Particle Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hSpectraZYAM3C=new TH2F("hSpectraZYAM3C",name,10,0,TPt1,10,plotmin,plotmax);
    if(NoTitle)hSpectraZYAM3C->SetTitle("");
    hSpectraZYAM3C->GetXaxis()->SetTitle("p_{T}^{Associated} (GeV/c) ");
    hSpectraZYAM3C->GetXaxis()->SetTitleColor(1);
    hSpectraZYAM3C->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dp_{T}}    ");
    SetTitles1D(hSpectraZYAM3C);
    hSpectraZYAM3C->Draw();
    cSpectraZYAM3C->SetLogy(1);
    gNPhiSpectraZYAM3->Draw("p");
    gAPhiSpectraZYAM3->Draw("p");
    gNEtaSpectraZYAM3->Draw("p");
    gAEtaSpectraZYAM3->Draw("p");
    if(DrawFit)gNPhiSpectraFit3->Draw("p");
    if(DrawFit)gAPhiSpectraFit3->Draw("p");
    if(DrawFit)gNEtaSpectraFit3->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.65,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.6,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.7,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_SpectraCombLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cSpectraZYAM3C->SaveAs(outName);


    cZt1ZYAM1C=new TCanvas("cZt1ZYAM1C","ZT MidPoint",800,600);
    SetMargins1D(cZt1ZYAM1C);
    sprintf(name,"Associated Particle z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt1ZYAM1C=new TH2F("hZt1ZYAM1C",name,10,0,1,10,plotmin,plotmax);
    if(NoTitle)hZt1ZYAM1C->SetTitle("");
    hZt1ZYAM1C->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig} ");
    hZt1ZYAM1C->GetXaxis()->SetTitleColor(1);
    hZt1ZYAM1C->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt1ZYAM1C);
    hZt1ZYAM1C->Draw();
    cZt1ZYAM1C->SetLogy(1);
    gNPhiZt1ZYAM1->Draw("p");
    gAPhiZt1ZYAM1->Draw("p");
    gNEtaZt1ZYAM1->Draw("p");
    gAEtaZt1ZYAM1->Draw("p");
    if(DrawFit)gNPhiZt1Fit1->Draw("p");
    if(DrawFit)gAPhiZt1Fit1->Draw("p");
    if(DrawFit)gNEtaZt1Fit1->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.65,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.6,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.7,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigCombMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cZt1ZYAM1C->SaveAs(outName);

    cZt1ZYAM2C=new TCanvas("cZt1ZYAM2C","Zt MidPoint",800,600);
    SetMargins1D(cZt1ZYAM2C);
    sprintf(name,"Associated Particle z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt1ZYAM2C=new TH2F("hZt1ZYAM2C",name,10,0,1,10,plotmin,plotmax);
    if(NoTitle)hZt1ZYAM2C->SetTitle("");
    hZt1ZYAM2C->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig} ");
    hZt1ZYAM2C->GetXaxis()->SetTitleColor(1);
    hZt1ZYAM2C->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
  
    SetTitles1D(hZt1ZYAM2C);
    hZt1ZYAM2C->Draw();
    cZt1ZYAM2C->SetLogy(1);
    gNPhiZt1ZYAM2->Draw("p");
    gAPhiZt1ZYAM2->Draw("p");
    gNEtaZt1ZYAM2->Draw("p");
    gAEtaZt1ZYAM2->Draw("p");
    if(DrawFit)gNPhiZt1Fit2->Draw("p");
    if(DrawFit)gAPhiZt1Fit2->Draw("p");
    if(DrawFit)gNEtaZt1Fit2->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.65,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.6,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.7,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigCombAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cZt1ZYAM2C->SaveAs(outName);

    cZt1ZYAM3C=new TCanvas("cZt1ZYAM3C","Zt1_3 MidPoint",800,600);
    SetMargins1D(cZt1ZYAM3C);
    sprintf(name,"Associated Particle z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt1ZYAM3C=new TH2F("hZt1ZYAM3C",name,10,0,1,10,plotmin,plotmax);
    if(NoTitle)hZt1ZYAM3C->SetTitle("");
    hZt1ZYAM3C->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig} ");
    hZt1ZYAM3C->GetXaxis()->SetTitleColor(1);
    hZt1ZYAM3C->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt1ZYAM3C);
    hZt1ZYAM3C->Draw();
    cZt1ZYAM3C->SetLogy(1);
    gNPhiZt1ZYAM3->Draw("p");
    gAPhiZt1ZYAM3->Draw("p");
    gNEtaZt1ZYAM3->Draw("p");
    gAEtaZt1ZYAM3->Draw("p");
    if(DrawFit)gNPhiZt1Fit3->Draw("p");
    if(DrawFit)gAPhiZt1Fit3->Draw("p");
    if(DrawFit)gNEtaZt1Fit3->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.65,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.6,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.7,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigCombLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cZt1ZYAM3C->SaveAs(outName);
    /*
    cZt2ZYAM1C=new TCanvas("cZt2ZYAM1C"," MidPoint",800,600);
    SetMargins1D(cZt2ZYAM1C);
    sprintf(name,"Associated Particle z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt2ZYAM1C=new TH2F("hZt2ZYAM1C",name,10,0,1,10,plotmin,plotmax);
    if(NoTitle)hZt2ZYAM1C->SetTitle("");
 
    hZt2ZYAM1C->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig} ");
    hZt2ZYAM1C->GetXaxis()->SetTitleColor(1);
    hZt2ZYAM1C->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt2ZYAM1C);
    hZt2ZYAM1C->Draw();
    cZt2ZYAM1C->SetLogy(1);
    gNPhiZt2ZYAM1->Draw("p");
    gAPhiZt2ZYAM1->Draw("p");
    gNEtaZt2ZYAM1->Draw("p");
    gAEtaZt2ZYAM1->Draw("p");
    //if(DrawFit)gNPhiZt2Fit1->Draw("p");
    //if(DrawFit)gAPhiZt2Fit1->Draw("p");
    //if(DrawFit)gNEtaZt2Fit1->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.75,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.7,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.85,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.75,"Near #Delta#etaFit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigCombMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,APtTPtMult,FitTit[DrawFit],filetype);
    cZt2ZYAM1C->SaveAs(outName);

    cZt2ZYAM2C=new TCanvas("cZt2ZYAM2C"," MidPoint",800,600);
    SetMargins1D(cZt2ZYAM2C);
    sprintf(name,"Associated Particle z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt2ZYAM2C=new TH2F("hZt2ZYAM2C",name,10,0,1,10,plotmin,plotmax);
    if(NoTitle)hZt2ZYAM2C->SetTitle("");
  
    hZt2ZYAM2C->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig} ");
    hZt2ZYAM2C->GetXaxis()->SetTitleColor(1);
    hZt2ZYAM2C->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
    SetTitles1D(hZt2ZYAM2C);
    hZt2ZYAM2C->Draw();
    cZt2ZYAM2C->SetLogy(1);
    gNPhiZt2ZYAM2->Draw("p");
    gAPhiZt2ZYAM2->Draw("p");
    gNEtaZt2ZYAM2->Draw("p");
    gAEtaZt2ZYAM2->Draw("p");
    //if(DrawFit)gNPhiZt2Fit2->Draw("p");
    //if(DrawFit)gAPhiZt2Fit2->Draw("p");
    // if(DrawFit)gNEtaZt2Fit2->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.75,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    //  if(DrawFit)keySymbol(.65,.7,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.85,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.75,"Near #Delta#etaFit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigCombAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,APtTPtMult,FitTit[DrawFit],filetype);
    cZt2ZYAM2C->SaveAs(outName);

    cZt2ZYAM3C=new TCanvas("cZt2ZYAM3C"," MidPoint",800,600);
    SetMargins1D(cZt2ZYAM3C);
    sprintf(name,"Associated Particle z_{T} Spectra %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hZt2ZYAM3C=new TH2F("hZt2ZYAM3C",name,10,0,1,10,plotmin,plotmax);
    if(NoTitle)hZt2ZYAM3C->SetTitle("");
    hZt2ZYAM3C->GetXaxis()->SetTitle("z_{T}=p_{T}^{Assoc}/p_{T}^{Trig} ");
    hZt2ZYAM3C->GetXaxis()->SetTitleColor(1);
    hZt2ZYAM3C->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}#frac{dN}{dz_{T}}    ");
  
    SetTitles1D(hZt2ZYAM3C);
    hZt2ZYAM3C->Draw();
    cZt2ZYAM3C->SetLogy(1);
    gNPhiZt2ZYAM3->Draw("p");
    gAPhiZt2ZYAM3->Draw("p");
    gNEtaZt2ZYAM3->Draw("p");
    gAEtaZt2ZYAM3->Draw("p");
    // if(DrawFit)gNPhiZt2Fit3->Draw("p");
    // if(DrawFit)gAPhiZt2Fit3->Draw("p");
    //  if(DrawFit)gNEtaZt2Fit3->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.75,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    //  if(DrawFit)keySymbol(.65,.7,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.85,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    // if(DrawFit)keySymbol(.65,.75,"Near #Delta#etaFit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_ZtTrigCombLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,APtTPtMult,FitTit[DrawFit],filetype);
    cZt2ZYAM3C->SaveAs(outName);
    */

    cWidth1=new TCanvas("cWidth1","Width MidPoint",800,600);
    SetMargins1D(cWidth1);
    sprintf(name,"Peak Width %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hWidth1=new TH2F("hWidth1",name,10,0,TPt1,10,plotminW,plotmaxW);
    if(NoTitle)hWidth1->SetTitle("");
    hWidth1->GetXaxis()->SetTitle("p_{T}^{Assoc} (GeV/c) ");
    hWidth1->GetXaxis()->SetTitleColor(1);
    hWidth1->GetYaxis()->SetTitle("#sigma (radians)   ");
    SetTitles1D(hWidth1);
    hWidth1->Draw();
    gNPhiWidth1->Draw("p");
    gAPhiWidth1->Draw("p");
    gNEtaWidth1->Draw("p");
    keySymbol(.65,.85,"Near #Delta#phi",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_WidthMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_C%d_%dM%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,filetype);
    cWidth1->SaveAs(outName);

    cWidth2=new TCanvas("cWidth2","Width MidPoint",800,600);
    SetMargins1D(cWidth2);
    sprintf(name,"Peak Width %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hWidth2=new TH2F("hWidth2",name,10,0,TPt1,10,plotminW,plotmaxW);
    if(NoTitle)hWidth2->SetTitle("");
    hWidth2->GetXaxis()->SetTitle("p_{T}^{Assoc} (GeV/c) ");
    hWidth2->GetXaxis()->SetTitleColor(1);
    hWidth2->GetYaxis()->SetTitle("#sigma (radians)   ");
    SetTitles1D(hWidth2);
    hWidth2->Draw();
    gNPhiWidth2->Draw("p");
    gAPhiWidth2->Draw("p");
    gNEtaWidth2->Draw("p");
    keySymbol(.65,.85,"Near #Delta#phi",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_WidthAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,filetype);
    cWidth2->SaveAs(outName);

    cWidth3=new TCanvas("cWidth3","Width MidPoint",800,600);
    SetMargins1D(cWidth3);
    sprintf(name,"Peak Width %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hWidth3=new TH2F("hWidth3",name,10,0,TPt1,10,plotminW,plotmaxW);
    if(NoTitle)hWidth3->SetTitle("");
    hWidth3->GetXaxis()->SetTitle("p_{T}^{Assoc} (GeV/c) ");
    hWidth3->GetXaxis()->SetTitleColor(1);
    hWidth3->GetYaxis()->SetTitle("#sigma (radians)   ");
    SetTitles1D(hWidth3);
    hWidth3->Draw();
    gNPhiWidth3->Draw("p");
    gAPhiWidth3->Draw("p");
    gNEtaWidth3->Draw("p");
    keySymbol(.65,.85,"Near #Delta#phi",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_WidthLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffCorr,APtTPtMult,filetype);
    cWidth3->SaveAs(outName);

  }

  if(APtTPtMult==1){
    if(HorizontalErrors==0){
      for(int mloop=0;mloop<nMainLoop;mloop++){
	NXPhi[1][mloop]=0;
	NXPhi[3][mloop]=0;
	NXPhi[5][mloop]=0;
	NXPhi[7][mloop]=0;
	NXPhi[9][mloop]=0;
	NXPhi[11][mloop]=0;
	NXPhi[13][mloop]=0;
	NXPhi[15][mloop]=0;
	NXPhi[17][mloop]=0;
	AXPhi[1][mloop]=0;
	AXPhi[3][mloop]=0;
	AXPhi[5][mloop]=0;
	AXPhi[7][mloop]=0;
	AXPhi[9][mloop]=0;
	AXPhi[11][mloop]=0;
	AXPhi[13][mloop]=0;
	AXPhi[15][mloop]=0;
	AXPhi[17][mloop]=0;
	NXPhiFit[1][mloop]=0;
	NXPhiFit[3][mloop]=0;
	NXPhiFit[5][mloop]=0;
	NXPhiFit[7][mloop]=0;
	NXPhiFit[9][mloop]=0;
	NXPhiFit[11][mloop]=0;
	NXPhiFit[13][mloop]=0;
	NXPhiFit[15][mloop]=0;
	NXPhiFit[17][mloop]=0;
	AXPhiFit[1][mloop]=0;
	AXPhiFit[3][mloop]=0;
	AXPhiFit[5][mloop]=0;
	AXPhiFit[7][mloop]=0;
	AXPhiFit[9][mloop]=0;
	AXPhiFit[11][mloop]=0;
	AXPhiFit[13][mloop]=0;
	AXPhiFit[15][mloop]=0;
	AXPhiFit[17][mloop]=0;
	NXEta[1][mloop]=0;
	NXEta[3][mloop]=0;
	NXEta[5][mloop]=0;
	NXEta[7][mloop]=0;
	NXEta[9][mloop]=0;
	NXEta[11][mloop]=0;
	NXEta[13][mloop]=0;
	NXEta[15][mloop]=0;
	NXEta[17][mloop]=0;
	AXEta[1][mloop]=0;
	AXEta[3][mloop]=0;
	AXEta[5][mloop]=0;
	AXEta[7][mloop]=0;
	AXEta[9][mloop]=0;
	AXEta[11][mloop]=0;
	AXEta[13][mloop]=0;
	AXEta[15][mloop]=0;
	AXEta[17][mloop]=0;
	NXEtaFit[1][mloop]=0;
	NXEtaFit[3][mloop]=0;
	NXEtaFit[5][mloop]=0;
	NXEtaFit[7][mloop]=0;
	NXEtaFit[9][mloop]=0;
	NXEtaFit[11][mloop]=0;
	NXEtaFit[13][mloop]=0;
	NXEtaFit[15][mloop]=0;
	NXEtaFit[17][mloop]=0;
	AXEtaFit[1][mloop]=0;
	AXEtaFit[3][mloop]=0;
	AXEtaFit[5][mloop]=0;
	AXEtaFit[7][mloop]=0;
	AXEtaFit[9][mloop]=0;
	AXEtaFit[11][mloop]=0;
	AXEtaFit[13][mloop]=0;
	AXEtaFit[15][mloop]=0;
	AXEtaFit[17][mloop]=0;
      }
    }
    gNYieldPhi=new TGraphErrors(nMainLoop,NXPhi[0],NYieldPhiZYAM[0],NXPhi[1],NYieldPhiZYAM[1]);
    gAYieldPhi=new TGraphErrors(nMainLoop,NXPhi[0],AYieldPhiZYAM[0],NXPhi[1],AYieldPhiZYAM[1]);
    gMYieldPhi=new TGraphErrors(nMainLoop,NXPhi[0],MYieldPhi[0],NXPhi[1],MYieldPhi[1]);
    gNYieldPhiPt=new TGraphErrors(nMainLoop,NXPhi[0],NYieldPhiZYAM[2],NXPhi[1],NYieldPhiZYAM[3]);
    gAYieldPhiPt=new TGraphErrors(nMainLoop,NXPhi[0],AYieldPhiZYAM[2],NXPhi[1],AYieldPhiZYAM[3]);
    gMYieldPhiPt=new TGraphErrors(nMainLoop,NXPhi[0],MYieldPhi[2],NXPhi[1],MYieldPhi[3]);

    gNYieldPhi2=new TGraphErrors(nMainLoop,NXPhi[2],NYieldPhiZYAM[0],NXPhi[3],NYieldPhiZYAM[1]);
    gAYieldPhi2=new TGraphErrors(nMainLoop,NXPhi[2],AYieldPhiZYAM[0],NXPhi[3],AYieldPhiZYAM[1]);
    gMYieldPhi2=new TGraphErrors(nMainLoop,NXPhi[2],MYieldPhi[0],NXPhi[3],MYieldPhi[1]);
    gNYieldPhiPt2=new TGraphErrors(nMainLoop,NXPhi[2],NYieldPhiZYAM[2],NXPhi[3],NYieldPhiZYAM[3]);
    gAYieldPhiPt2=new TGraphErrors(nMainLoop,NXPhi[2],AYieldPhiZYAM[2],NXPhi[3],AYieldPhiZYAM[3]);
    gMYieldPhiPt2=new TGraphErrors(nMainLoop,NXPhi[2],MYieldPhi[2],NXPhi[3],MYieldPhi[3]);

    gNYieldPhi3=new TGraphErrors(nMainLoop,NXPhi[4],NYieldPhiZYAM[0],NXPhi[5],NYieldPhiZYAM[1]);
    gAYieldPhi3=new TGraphErrors(nMainLoop,NXPhi[4],AYieldPhiZYAM[0],NXPhi[5],AYieldPhiZYAM[1]);
    gMYieldPhi3=new TGraphErrors(nMainLoop,NXPhi[4],MYieldPhi[0],NXPhi[5],MYieldPhi[1]);
    gNYieldPhiPt3=new TGraphErrors(nMainLoop,NXPhi[4],NYieldPhiZYAM[2],NXPhi[5],NYieldPhiZYAM[3]);
    gAYieldPhiPt3=new TGraphErrors(nMainLoop,NXPhi[4],AYieldPhiZYAM[2],NXPhi[5],AYieldPhiZYAM[3]);
    gMYieldPhiPt3=new TGraphErrors(nMainLoop,NXPhi[4],MYieldPhi[2],NXPhi[5],MYieldPhi[3]);

    gNYieldPhiFit=new TGraphErrors(nMainLoop,NXPhiFit[0],NYieldPhiFit[0],NXPhiFit[1],NYieldPhiFit[1]);
    gAYieldPhiFit=new TGraphErrors(nMainLoop,NXPhiFit[0],AYieldPhiFit[0],NXPhiFit[1],AYieldPhiFit[1]);
    gMYieldPhiFit=new TGraphErrors(nMainLoop,NXPhiFit[0],MYieldPhiFit[0],NXPhiFit[1],MYieldPhiFit[1]);
    gNYieldPhiFitPt=new TGraphErrors(nMainLoop,NXPhiFit[0],NYieldPhiFit[2],NXPhiFit[1],NYieldPhiFit[3]);
    gAYieldPhiFitPt=new TGraphErrors(nMainLoop,NXPhiFit[0],AYieldPhiFit[2],NXPhiFit[1],AYieldPhiFit[3]);
    gMYieldPhiFitPt=new TGraphErrors(nMainLoop,NXPhiFit[0],MYieldPhiFit[2],NXPhiFit[1],MYieldPhiFit[3]);

    gNYieldPhiFit2=new TGraphErrors(nMainLoop,NXPhiFit[2],NYieldPhiFit[0],NXPhiFit[3],NYieldPhiFit[1]);
    gAYieldPhiFit2=new TGraphErrors(nMainLoop,NXPhiFit[2],AYieldPhiFit[0],NXPhiFit[3],AYieldPhiFit[1]);
    gMYieldPhiFit2=new TGraphErrors(nMainLoop,NXPhiFit[2],MYieldPhiFit[0],NXPhiFit[3],MYieldPhiFit[1]);
    gNYieldPhiFitPt2=new TGraphErrors(nMainLoop,NXPhiFit[2],NYieldPhiFit[2],NXPhiFit[3],NYieldPhiFit[3]);
    gAYieldPhiFitPt2=new TGraphErrors(nMainLoop,NXPhiFit[2],AYieldPhiFit[2],NXPhiFit[3],AYieldPhiFit[3]);
    gMYieldPhiFitPt2=new TGraphErrors(nMainLoop,NXPhiFit[2],MYieldPhiFit[2],NXPhiFit[3],MYieldPhiFit[3]);

    gNYieldPhiFit3=new TGraphErrors(nMainLoop,NXPhiFit[4],NYieldPhiFit[0],NXPhiFit[5],NYieldPhiFit[1]);
    gAYieldPhiFit3=new TGraphErrors(nMainLoop,NXPhiFit[4],AYieldPhiFit[0],NXPhiFit[5],AYieldPhiFit[1]);
    gMYieldPhiFit3=new TGraphErrors(nMainLoop,NXPhiFit[4],MYieldPhiFit[0],NXPhiFit[5],MYieldPhiFit[1]);
    gNYieldPhiFitPt3=new TGraphErrors(nMainLoop,NXPhiFit[4],NYieldPhiFit[2],NXPhiFit[5],NYieldPhiFit[3]);
    gAYieldPhiFitPt3=new TGraphErrors(nMainLoop,NXPhiFit[4],AYieldPhiFit[2],NXPhiFit[5],AYieldPhiFit[3]);
    gMYieldPhiFitPt3=new TGraphErrors(nMainLoop,NXPhiFit[4],MYieldPhiFit[2],NXPhiFit[5],MYieldPhiFit[3]);
   
    gNPhiWidth1=new TGraphErrors(nMainLoop,NXPhiFit[0],NWidthPhi[0],NXPhiFit[1],NWidthPhi[1]);
    gNPhiWidth2=new TGraphErrors(nMainLoop,NXPhiFit[2],NWidthPhi[0],NXPhiFit[3],NWidthPhi[1]);
    gNPhiWidth3=new TGraphErrors(nMainLoop,NXPhiFit[4],NWidthPhi[0],NXPhiFit[5],NWidthPhi[1]);
    gAPhiWidth1=new TGraphErrors(nMainLoop,AXPhiFit[0],AWidthPhi[0],AXPhiFit[1],AWidthPhi[1]);
    gAPhiWidth2=new TGraphErrors(nMainLoop,AXPhiFit[2],AWidthPhi[0],AXPhiFit[3],AWidthPhi[1]);
    gAPhiWidth3=new TGraphErrors(nMainLoop,AXPhiFit[4],AWidthPhi[0],AXPhiFit[5],AWidthPhi[1]);

    gNPhiJt1=new TGraphErrors(nMainLoop,NXPhiFit[0],NJtPhi[0],NXPhiFit[1],NJtPhi[1]);
    gNPhiJt2=new TGraphErrors(nMainLoop,NXPhiFit[2],NJtPhi[0],NXPhiFit[3],NJtPhi[1]);
    gNPhiJt3=new TGraphErrors(nMainLoop,NXPhiFit[4],NJtPhi[0],NXPhiFit[5],NJtPhi[1]);

    gNYieldEta=new TGraphErrors(nMainLoop,NXEta[0],NYieldEtaZYAM[0],NXEta[1],NYieldEtaZYAM[1]);
    gAYieldEta=new TGraphErrors(nMainLoop,NXEta[0],AYieldEtaZYAM[0],NXEta[1],AYieldEtaZYAM[1]);
    gNYieldEtaPt=new TGraphErrors(nMainLoop,NXEta[0],NYieldEtaZYAM[2],NXEta[1],NYieldEtaZYAM[3]);
    gAYieldEtaPt=new TGraphErrors(nMainLoop,NXEta[0],AYieldEtaZYAM[2],NXEta[1],AYieldEtaZYAM[3]);
  
    gNYieldEta2=new TGraphErrors(nMainLoop,NXEta[2],NYieldEtaZYAM[0],NXEta[3],NYieldEtaZYAM[1]);
    gAYieldEta2=new TGraphErrors(nMainLoop,NXEta[2],AYieldEtaZYAM[0],NXEta[3],AYieldEtaZYAM[1]);
    gNYieldEtaPt2=new TGraphErrors(nMainLoop,NXEta[2],NYieldEtaZYAM[2],NXEta[3],NYieldEtaZYAM[3]);
    gAYieldEtaPt2=new TGraphErrors(nMainLoop,NXEta[2],AYieldEtaZYAM[2],NXEta[3],AYieldEtaZYAM[3]);

    gNYieldEta3=new TGraphErrors(nMainLoop,NXEta[4],NYieldEtaZYAM[0],NXEta[5],NYieldEtaZYAM[1]);
    gAYieldEta3=new TGraphErrors(nMainLoop,NXEta[4],AYieldEtaZYAM[0],NXEta[5],AYieldEtaZYAM[1]);
    gNYieldEtaPt3=new TGraphErrors(nMainLoop,NXEta[4],NYieldEtaZYAM[2],NXEta[5],NYieldEtaZYAM[3]);
    gAYieldEtaPt3=new TGraphErrors(nMainLoop,NXEta[4],AYieldEtaZYAM[2],NXEta[5],AYieldEtaZYAM[3]);
 
    gNYieldEtaFit=new TGraphErrors(nMainLoop,NXEtaFit[0],NYieldEtaFit[0],NXEtaFit[1],NYieldEtaFit[1]);
    // gAYieldEtaFit=new TGraphErrors(nMainLoop,NXEtaFit[0],AYieldEtaFit[0],NXEtaFit[1],AYieldEtaFit[1]);
    gNYieldEtaFitPt=new TGraphErrors(nMainLoop,NXEtaFit[0],NYieldEtaFit[2],NXEtaFit[1],NYieldEtaFit[3]);
    // gAYieldEtaFitPt=new TGraphErrors(nMainLoop,NXEtaFit[0],AYieldEtaFit[2],NXEtaFit[1],AYieldEtaFit[3]);

    gNYieldEtaFit2=new TGraphErrors(nMainLoop,NXEtaFit[2],NYieldEtaFit[0],NXEtaFit[3],NYieldEtaFit[1]);
    //gAYieldEtaFit2=new TGraphErrors(nMainLoop,NXEtaFit[2],AYieldEtaFit[0],NXEtaFit[3],AYieldEtaFit[1]);
    gNYieldEtaFitPt2=new TGraphErrors(nMainLoop,NXEtaFit[2],NYieldEtaFit[2],NXEtaFit[3],NYieldEtaFit[3]);
    // gAYieldEtaFitPt2=new TGraphErrors(nMainLoop,NXEtaFit[2],AYieldEtaFit[2],NXEtaFit[3],AYieldEtaFit[3]);

    gNYieldEtaFit3=new TGraphErrors(nMainLoop,NXEtaFit[4],NYieldEtaFit[0],NXEtaFit[5],NYieldEtaFit[1]);
    // gAYieldEtaFit3=new TGraphErrors(nMainLoop,NXEtaFit[4],AYieldEtaFit[0],NXEtaFit[5],AYieldEtaFit[1]);
    gNYieldEtaFitPt3=new TGraphErrors(nMainLoop,NXEtaFit[4],NYieldEtaFit[2],NXEtaFit[5],NYieldEtaFit[3]);
    // gAYieldEtaFitPt3=new TGraphErrors(nMainLoop,NXEtaFit[4],AYieldEtaFit[2],NXEtaFit[5],AYieldEtaFit[3]);

    gNEtaWidth1=new TGraphErrors(nMainLoop,NXPhiFit[0],NWidthEta[0],NXPhiFit[1],NWidthEta[1]);
    gNEtaWidth2=new TGraphErrors(nMainLoop,NXPhiFit[2],NWidthEta[0],NXPhiFit[3],NWidthEta[1]);
    gNEtaWidth3=new TGraphErrors(nMainLoop,NXPhiFit[4],NWidthEta[0],NXPhiFit[5],NWidthEta[1]);

    gNEtaJt1=new TGraphErrors(nMainLoop,NXEtaFit[0],NJtEta[0],NXEtaFit[1],NJtEta[1]);
    gNEtaJt2=new TGraphErrors(nMainLoop,NXEtaFit[2],NJtEta[0],NXEtaFit[3],NJtEta[1]);
    gNEtaJt3=new TGraphErrors(nMainLoop,NXEtaFit[4],NJtEta[0],NXEtaFit[5],NJtEta[1]);

    gNAvePhi1=new TGraphErrors(nMainLoop,NXPhi[0],NAvePhi[0],NXPhi[1],NAvePhi[1]);
    gNAvePhi2=new TGraphErrors(nMainLoop,NXPhi[2],NAvePhi[0],NXPhi[3],NAvePhi[1]);
    gNAvePhi3=new TGraphErrors(nMainLoop,NXPhi[4],NAvePhi[0],NXPhi[5],NAvePhi[1]);

    gAAvePhi1=new TGraphErrors(nMainLoop,NXPhi[0],AAvePhi[0],NXPhi[1],AAvePhi[1]);
    gAAvePhi2=new TGraphErrors(nMainLoop,NXPhi[2],AAvePhi[0],NXPhi[3],AAvePhi[1]);
    gAAvePhi3=new TGraphErrors(nMainLoop,NXPhi[4],AAvePhi[0],NXPhi[5],AAvePhi[1]);

    gMAvePhi1=new TGraphErrors(nMainLoop,NXPhi[0],MAvePhi[0],NXPhi[1],MAvePhi[1]);
    gMAvePhi2=new TGraphErrors(nMainLoop,NXPhi[2],MAvePhi[0],NXPhi[3],MAvePhi[1]);
    gMAvePhi3=new TGraphErrors(nMainLoop,NXPhi[4],MAvePhi[0],NXPhi[5],MAvePhi[1]);

    gNAveEta1=new TGraphErrors(nMainLoop,NXEta[0],NAveEta[0],NXEta[1],NAveEta[1]);
    gNAveEta2=new TGraphErrors(nMainLoop,NXEta[2],NAveEta[0],NXEta[3],NAveEta[1]);
    gNAveEta3=new TGraphErrors(nMainLoop,NXEta[4],NAveEta[0],NXEta[5],NAveEta[1]);

    gAAveEta1=new TGraphErrors(nMainLoop,NXEta[0],AAveEta[0],NXEta[1],AAveEta[1]);
    gAAveEta2=new TGraphErrors(nMainLoop,NXEta[2],AAveEta[0],NXEta[3],AAveEta[1]);
    gAAveEta3=new TGraphErrors(nMainLoop,NXEta[4],AAveEta[0],NXEta[5],AAveEta[1]);

    gNAvePhiFit1=new TGraphErrors(nMainLoop,NXPhi[0],NAvePhiFit[0],NXPhi[1],NAvePhiFit[1]);
    gNAvePhiFit2=new TGraphErrors(nMainLoop,NXPhi[2],NAvePhiFit[0],NXPhi[3],NAvePhiFit[1]);
    gNAvePhiFit3=new TGraphErrors(nMainLoop,NXPhi[4],NAvePhiFit[0],NXPhi[5],NAvePhiFit[1]);

    gAAvePhiFit1=new TGraphErrors(nMainLoop,NXPhi[0],AAvePhiFit[0],NXPhi[1],AAvePhiFit[1]);
    gAAvePhiFit2=new TGraphErrors(nMainLoop,NXPhi[2],AAvePhiFit[0],NXPhi[3],AAvePhiFit[1]);
    gAAvePhiFit3=new TGraphErrors(nMainLoop,NXPhi[4],AAvePhiFit[0],NXPhi[5],AAvePhiFit[1]);

    gMAvePhiFit1=new TGraphErrors(nMainLoop,NXPhi[0],MAvePhiFit[0],NXPhi[1],MAvePhiFit[1]);
    gMAvePhiFit2=new TGraphErrors(nMainLoop,NXPhi[2],MAvePhiFit[0],NXPhi[3],MAvePhiFit[1]);
    gMAvePhiFit3=new TGraphErrors(nMainLoop,NXPhi[4],MAvePhiFit[0],NXPhi[5],MAvePhiFit[1]);

    gNAveEtaFit1=new TGraphErrors(nMainLoop,NXEta[0],NAveEtaFit[0],NXEta[1],NAveEtaFit[1]);
    gNAveEtaFit2=new TGraphErrors(nMainLoop,NXEta[2],NAveEtaFit[0],NXEta[3],NAveEtaFit[1]);
    gNAveEtaFit3=new TGraphErrors(nMainLoop,NXEta[4],NAveEtaFit[0],NXEta[5],NAveEtaFit[1]);
   
    gNYieldPhi->SetMarkerStyle(MarkerNearPhi);
    gAYieldPhi->SetMarkerStyle(MarkerAwayPhi);
    gMYieldPhi->SetMarkerStyle(MarkerMin);
    gNYieldPhiPt->SetMarkerStyle(MarkerNearPhi);
    gAYieldPhiPt->SetMarkerStyle(MarkerAwayPhi);
    gMYieldPhiPt->SetMarkerStyle(MarkerMin);
    gNYieldPhi2->SetMarkerStyle(MarkerNearPhi);
    gAYieldPhi2->SetMarkerStyle(MarkerAwayPhi);
    gMYieldPhi2->SetMarkerStyle(MarkerMin);
    gNYieldPhiPt2->SetMarkerStyle(MarkerNearPhi);
    gAYieldPhiPt2->SetMarkerStyle(MarkerAwayPhi);
    gMYieldPhiPt2->SetMarkerStyle(MarkerMin);
    gNYieldPhi3->SetMarkerStyle(MarkerNearPhi);
    gAYieldPhi3->SetMarkerStyle(MarkerAwayPhi);
    gMYieldPhi3->SetMarkerStyle(MarkerMin);
    gNYieldPhiPt3->SetMarkerStyle(MarkerNearPhi);
    gAYieldPhiPt3->SetMarkerStyle(MarkerAwayPhi);
    gMYieldPhiPt3->SetMarkerStyle(MarkerMin);
    gNYieldPhi->SetMarkerColor(ColorNearPhi);
    gAYieldPhi->SetMarkerColor(ColorAwayPhi);
    gMYieldPhi->SetMarkerColor(ColorMin);
    gNYieldPhi->SetLineColor(ColorNearPhi);
    gAYieldPhi->SetLineColor(ColorAwayPhi);
    gMYieldPhi->SetLineColor(ColorMin);
    gNYieldPhiPt->SetMarkerColor(ColorNearPhi);
    gAYieldPhiPt->SetMarkerColor(ColorAwayPhi);
    gMYieldPhiPt->SetMarkerColor(ColorMin);
    gNYieldPhiPt->SetLineColor(ColorNearPhi);
    gAYieldPhiPt->SetLineColor(ColorAwayPhi);
    gMYieldPhiPt->SetLineColor(ColorMin);
    gNYieldPhi2->SetMarkerColor(ColorNearPhi);
    gAYieldPhi2->SetMarkerColor(ColorAwayPhi);
    gMYieldPhi2->SetMarkerColor(ColorMin);
    gNYieldPhi2->SetLineColor(ColorNearPhi);
    gAYieldPhi2->SetLineColor(ColorAwayPhi);
    gMYieldPhi2->SetLineColor(ColorMin);
    gNYieldPhiPt2->SetMarkerColor(ColorNearPhi);
    gAYieldPhiPt2->SetMarkerColor(ColorAwayPhi);
    gMYieldPhiPt2->SetMarkerColor(ColorMin);
    gNYieldPhiPt2->SetLineColor(ColorNearPhi);
    gAYieldPhiPt2->SetLineColor(ColorAwayPhi);
    gMYieldPhiPt2->SetLineColor(ColorMin);
    gNYieldPhi3->SetMarkerColor(ColorNearPhi);
    gAYieldPhi3->SetMarkerColor(ColorAwayPhi);
    gMYieldPhi3->SetMarkerColor(ColorMin);
    gNYieldPhi3->SetLineColor(ColorNearPhi);
    gAYieldPhi3->SetLineColor(ColorAwayPhi);
    gMYieldPhi3->SetLineColor(ColorMin);
    gNYieldPhiPt3->SetMarkerColor(ColorNearPhi);
    gAYieldPhiPt3->SetMarkerColor(ColorAwayPhi);
    gMYieldPhiPt3->SetMarkerColor(ColorMin);
    gNYieldPhiPt3->SetLineColor(ColorNearPhi);
    gAYieldPhiPt3->SetLineColor(ColorAwayPhi);
    gMYieldPhiPt3->SetLineColor(ColorMin);
    gNYieldPhi->SetMarkerSize(MarkerSize);
    gAYieldPhi->SetMarkerSize(MarkerSize);
    gMYieldPhi->SetMarkerSize(MarkerSize);
    gNYieldPhiPt->SetMarkerSize(MarkerSize);
    gAYieldPhiPt->SetMarkerSize(MarkerSize);
    gMYieldPhiPt->SetMarkerSize(MarkerSize);
    gNYieldPhi2->SetMarkerSize(MarkerSize);
    gAYieldPhi2->SetMarkerSize(MarkerSize);
    gMYieldPhi2->SetMarkerSize(MarkerSize);
    gNYieldPhiPt2->SetMarkerSize(MarkerSize);
    gAYieldPhiPt2->SetMarkerSize(MarkerSize);
    gMYieldPhiPt2->SetMarkerSize(MarkerSize);
    gNYieldPhi3->SetMarkerSize(MarkerSize);
    gAYieldPhi3->SetMarkerSize(MarkerSize);
    gMYieldPhi3->SetMarkerSize(MarkerSize);
    gNYieldPhiPt3->SetMarkerSize(MarkerSize);
    gAYieldPhiPt3->SetMarkerSize(MarkerSize);
    gNYieldPhiPt3->SetMarkerSize(MarkerSize);

    gNYieldPhiFit->SetMarkerStyle(MarkerNearPhiFit);
    gNYieldPhiFit->SetMarkerSize(MarkerSize);
    gNYieldPhiFit->SetMarkerColor(ColorNearPhiFit);
    gNYieldPhiFit->SetLineColor(ColorNearPhiFit);
    gNYieldPhiFit2->SetMarkerStyle(MarkerNearPhiFit);
    gNYieldPhiFit2->SetMarkerSize(MarkerSize);
    gNYieldPhiFit2->SetMarkerColor(ColorNearPhiFit);
    gNYieldPhiFit2->SetLineColor(ColorNearPhiFit);
    gNYieldPhiFit3->SetMarkerStyle(MarkerNearPhiFit);
    gNYieldPhiFit3->SetMarkerSize(MarkerSize);
    gNYieldPhiFit3->SetMarkerColor(ColorNearPhiFit);
    gNYieldPhiFit3->SetLineColor(ColorNearPhiFit);

    gAYieldPhiFit->SetMarkerStyle(MarkerAwayPhiFit);
    gAYieldPhiFit->SetMarkerSize(MarkerSize);
    gAYieldPhiFit->SetMarkerColor(ColorAwayPhiFit);
    gAYieldPhiFit->SetLineColor(ColorAwayPhiFit);
    gAYieldPhiFit2->SetMarkerStyle(MarkerAwayPhiFit);
    gAYieldPhiFit2->SetMarkerSize(MarkerSize);
    gAYieldPhiFit2->SetMarkerColor(ColorAwayPhiFit);
    gAYieldPhiFit2->SetLineColor(ColorAwayPhiFit);
    gAYieldPhiFit3->SetMarkerStyle(MarkerAwayPhiFit);
    gAYieldPhiFit3->SetMarkerSize(MarkerSize);
    gAYieldPhiFit3->SetMarkerColor(ColorAwayPhiFit);
    gAYieldPhiFit3->SetLineColor(ColorAwayPhiFit);

    gMYieldPhiFit->SetMarkerStyle(MarkerMinFit);
    gMYieldPhiFit->SetMarkerSize(MarkerSize);
    gMYieldPhiFit->SetMarkerColor(ColorMinFit);
    gMYieldPhiFit->SetLineColor(ColorMinFit);
    gMYieldPhiFit2->SetMarkerStyle(MarkerMinFit);
    gMYieldPhiFit2->SetMarkerSize(MarkerSize);
    gMYieldPhiFit2->SetMarkerColor(ColorMinFit);
    gMYieldPhiFit2->SetLineColor(ColorMinFit);
    gMYieldPhiFit3->SetMarkerStyle(MarkerMinFit);
    gMYieldPhiFit3->SetMarkerSize(MarkerSize);
    gMYieldPhiFit3->SetMarkerColor(ColorMinFit);
    gMYieldPhiFit3->SetLineColor(ColorMinFit);

    gNYieldPhiFitPt->SetMarkerStyle(MarkerNearPhiFit);
    gNYieldPhiFitPt->SetMarkerSize(MarkerSize);
    gNYieldPhiFitPt->SetMarkerColor(ColorNearPhiFit);
    gNYieldPhiFitPt->SetLineColor(ColorNearPhiFit);
    gNYieldPhiFitPt2->SetMarkerStyle(MarkerNearPhiFit);
    gNYieldPhiFitPt2->SetMarkerSize(MarkerSize);
    gNYieldPhiFitPt2->SetMarkerColor(ColorNearPhiFit);
    gNYieldPhiFitPt2->SetLineColor(ColorNearPhiFit);
    gNYieldPhiFitPt3->SetMarkerStyle(MarkerNearPhiFit);
    gNYieldPhiFitPt3->SetMarkerSize(MarkerSize);
    gNYieldPhiFitPt3->SetMarkerColor(ColorNearPhiFit);
    gNYieldPhiFitPt3->SetLineColor(ColorNearPhiFit);

    gAYieldPhiFitPt->SetMarkerStyle(MarkerAwayPhiFit);
    gAYieldPhiFitPt->SetMarkerSize(MarkerSize);
    gAYieldPhiFitPt->SetMarkerColor(ColorAwayPhiFit);
    gAYieldPhiFitPt->SetLineColor(ColorAwayPhiFit);
    gAYieldPhiFitPt2->SetMarkerStyle(MarkerAwayPhiFit);
    gAYieldPhiFitPt2->SetMarkerSize(MarkerSize);
    gAYieldPhiFitPt2->SetMarkerColor(ColorAwayPhiFit);
    gAYieldPhiFitPt2->SetLineColor(ColorAwayPhiFit);
    gAYieldPhiFitPt3->SetMarkerStyle(MarkerAwayPhiFit);
    gAYieldPhiFitPt3->SetMarkerSize(MarkerSize);
    gAYieldPhiFitPt3->SetMarkerColor(ColorAwayPhiFit);
    gAYieldPhiFitPt3->SetLineColor(ColorAwayPhiFit);

    gMYieldPhiFitPt->SetMarkerStyle(MarkerMinFit);
    gMYieldPhiFitPt->SetMarkerSize(MarkerSize);
    gMYieldPhiFitPt->SetMarkerColor(ColorMinFit);
    gMYieldPhiFitPt->SetLineColor(ColorMinFit);
    gMYieldPhiFitPt2->SetMarkerStyle(MarkerMinFit);
    gMYieldPhiFitPt2->SetMarkerSize(MarkerSize);
    gMYieldPhiFitPt2->SetMarkerColor(ColorMinFit);
    gMYieldPhiFitPt2->SetLineColor(ColorMinFit);
    gMYieldPhiFitPt3->SetMarkerStyle(MarkerMinFit);
    gMYieldPhiFitPt3->SetMarkerSize(MarkerSize);
    gMYieldPhiFitPt3->SetMarkerColor(ColorMinFit);
    gMYieldPhiFitPt3->SetLineColor(ColorMinFit);

    gNPhiWidth1->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiWidth1->SetMarkerSize(MarkerSize);
    gNPhiWidth1->SetMarkerColor(ColorNearPhiFit);
    gNPhiWidth1->SetLineColor(ColorNearPhiFit);
    gNPhiWidth2->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiWidth2->SetMarkerSize(MarkerSize);
    gNPhiWidth2->SetMarkerColor(ColorNearPhiFit);
    gNPhiWidth2->SetLineColor(ColorNearPhiFit);
    gNPhiWidth3->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiWidth3->SetMarkerSize(MarkerSize);
    gNPhiWidth3->SetMarkerColor(ColorNearPhiFit);
    gNPhiWidth3->SetLineColor(ColorNearPhiFit);

    gAPhiWidth1->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiWidth1->SetMarkerSize(MarkerSize);
    gAPhiWidth1->SetMarkerColor(ColorAwayPhiFit);
    gAPhiWidth1->SetLineColor(ColorAwayPhiFit);
    gAPhiWidth2->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiWidth2->SetMarkerSize(MarkerSize);
    gAPhiWidth2->SetMarkerColor(ColorAwayPhiFit);
    gAPhiWidth2->SetLineColor(ColorAwayPhiFit);
    gAPhiWidth3->SetMarkerStyle(MarkerAwayPhiFit);
    gAPhiWidth3->SetMarkerSize(MarkerSize);
    gAPhiWidth3->SetMarkerColor(ColorAwayPhiFit);
    gAPhiWidth3->SetLineColor(ColorAwayPhiFit);

    gNPhiJt1->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiJt1->SetMarkerSize(MarkerSize);
    gNPhiJt1->SetMarkerColor(ColorNearPhiFit);
    gNPhiJt1->SetLineColor(ColorNearPhiFit);
    gNPhiJt2->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiJt2->SetMarkerSize(MarkerSize);
    gNPhiJt2->SetMarkerColor(ColorNearPhiFit);
    gNPhiJt2->SetLineColor(ColorNearPhiFit);
    gNPhiJt3->SetMarkerStyle(MarkerNearPhiFit);
    gNPhiJt3->SetMarkerSize(MarkerSize);
    gNPhiJt3->SetMarkerColor(ColorNearPhiFit);
    gNPhiJt3->SetLineColor(ColorNearPhiFit);


    gNYieldEta->SetMarkerStyle(MarkerNearEta);
    gAYieldEta->SetMarkerStyle(MarkerAwayEta);
    gNYieldEtaPt->SetMarkerStyle(MarkerNearEta);
    gAYieldEtaPt->SetMarkerStyle(MarkerAwayEta);
    gNYieldEta2->SetMarkerStyle(MarkerNearEta);
    gAYieldEta2->SetMarkerStyle(MarkerAwayEta);
    gNYieldEtaPt2->SetMarkerStyle(MarkerNearEta);
    gAYieldEtaPt2->SetMarkerStyle(MarkerAwayEta);
    gNYieldEta3->SetMarkerStyle(MarkerNearEta);
    gAYieldEta3->SetMarkerStyle(MarkerAwayEta);
    gNYieldEtaPt3->SetMarkerStyle(MarkerNearEta);
    gAYieldEtaPt3->SetMarkerStyle(MarkerAwayEta);
    gNYieldEta->SetMarkerColor(ColorNearEta);
    gAYieldEta->SetMarkerColor(ColorAwayEta);
    gNYieldEta->SetLineColor(ColorNearEta);
    gAYieldEta->SetLineColor(ColorAwayEta);
    gNYieldEtaPt->SetMarkerColor(ColorNearEta);
    gAYieldEtaPt->SetMarkerColor(ColorAwayEta);
    gNYieldEtaPt->SetLineColor(ColorNearEta);
    gAYieldEtaPt->SetLineColor(ColorAwayEta);
    gNYieldEta2->SetMarkerColor(ColorNearEta);
    gAYieldEta2->SetMarkerColor(ColorAwayEta);
    gNYieldEta2->SetLineColor(ColorNearEta);
    gAYieldEta2->SetLineColor(ColorAwayEta);
    gNYieldEtaPt2->SetMarkerColor(ColorNearEta);
    gAYieldEtaPt2->SetMarkerColor(ColorAwayEta);
    gNYieldEtaPt2->SetLineColor(ColorNearEta);
    gAYieldEtaPt2->SetLineColor(ColorAwayEta);
    gNYieldEta3->SetMarkerColor(ColorNearEta);
    gAYieldEta3->SetMarkerColor(ColorAwayEta);
    gNYieldEta3->SetLineColor(ColorNearEta);
    gAYieldEta3->SetLineColor(ColorAwayEta);
    gNYieldEtaPt3->SetMarkerColor(ColorNearEta);
    gAYieldEtaPt3->SetMarkerColor(ColorAwayEta);
    gNYieldEtaPt3->SetLineColor(ColorNearEta);
    gAYieldEtaPt3->SetLineColor(ColorAwayEta);
    gNYieldEta->SetMarkerSize(MarkerSize);
    gAYieldEta->SetMarkerSize(MarkerSize);
    gNYieldEtaPt->SetMarkerSize(MarkerSize);
    gAYieldEtaPt->SetMarkerSize(MarkerSize);
    gNYieldEta2->SetMarkerSize(MarkerSize);
    gAYieldEta2->SetMarkerSize(MarkerSize);
    gNYieldEtaPt2->SetMarkerSize(MarkerSize);
    gAYieldEtaPt2->SetMarkerSize(MarkerSize);
    gNYieldEta3->SetMarkerSize(MarkerSize);
    gAYieldEta3->SetMarkerSize(MarkerSize);
    gNYieldEtaPt3->SetMarkerSize(MarkerSize);
    gAYieldEtaPt3->SetMarkerSize(MarkerSize);
    gNYieldEtaPt3->SetMarkerSize(MarkerSize);

    gNYieldEtaFit->SetMarkerStyle(MarkerNearEtaFit);
    gNYieldEtaFit->SetMarkerSize(MarkerSize);
    gNYieldEtaFit->SetMarkerColor(ColorNearEtaFit);
    gNYieldEtaFit->SetLineColor(ColorNearEtaFit);
    gNYieldEtaFit2->SetMarkerStyle(MarkerNearEtaFit);
    gNYieldEtaFit2->SetMarkerSize(MarkerSize);
    gNYieldEtaFit2->SetMarkerColor(ColorNearEtaFit);
    gNYieldEtaFit2->SetLineColor(ColorNearEtaFit);
    gNYieldEtaFit3->SetMarkerStyle(MarkerNearEtaFit);
    gNYieldEtaFit3->SetMarkerSize(MarkerSize);
    gNYieldEtaFit3->SetMarkerColor(ColorNearEtaFit);
    gNYieldEtaFit3->SetLineColor(ColorNearEtaFit);

    gNYieldEtaFitPt->SetMarkerStyle(MarkerNearEtaFit);
    gNYieldEtaFitPt->SetMarkerSize(MarkerSize);
    gNYieldEtaFitPt->SetMarkerColor(ColorNearEtaFit);
    gNYieldEtaFitPt->SetLineColor(ColorNearEtaFit);
    gNYieldEtaFitPt2->SetMarkerStyle(MarkerNearEtaFit);
    gNYieldEtaFitPt2->SetMarkerSize(MarkerSize);
    gNYieldEtaFitPt2->SetMarkerColor(ColorNearEtaFit);
    gNYieldEtaFitPt2->SetLineColor(ColorNearEtaFit);
    gNYieldEtaFitPt3->SetMarkerStyle(MarkerNearEtaFit);
    gNYieldEtaFitPt3->SetMarkerSize(MarkerSize);
    gNYieldEtaFitPt3->SetMarkerColor(ColorNearEtaFit);
    gNYieldEtaFitPt3->SetLineColor(ColorNearEtaFit);

    gNEtaWidth1->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaWidth1->SetMarkerSize(MarkerSize);
    gNEtaWidth1->SetMarkerColor(ColorNearEtaFit);
    gNEtaWidth1->SetLineColor(ColorNearEtaFit);
    gNEtaWidth2->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaWidth2->SetMarkerSize(MarkerSize);
    gNEtaWidth2->SetMarkerColor(ColorNearEtaFit);
    gNEtaWidth2->SetLineColor(ColorNearEtaFit);
    gNEtaWidth3->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaWidth3->SetMarkerSize(MarkerSize);
    gNEtaWidth3->SetMarkerColor(ColorNearEtaFit);
    gNEtaWidth3->SetLineColor(ColorNearEtaFit);

    gNEtaJt1->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaJt1->SetMarkerSize(MarkerSize);
    gNEtaJt1->SetMarkerColor(ColorNearEtaFit);
    gNEtaJt1->SetLineColor(ColorNearEtaFit);
    gNEtaJt2->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaJt2->SetMarkerSize(MarkerSize);
    gNEtaJt2->SetMarkerColor(ColorNearEtaFit);
    gNEtaJt2->SetLineColor(ColorNearEtaFit);
    gNEtaJt3->SetMarkerStyle(MarkerNearEtaFit);
    gNEtaJt3->SetMarkerSize(MarkerSize);
    gNEtaJt3->SetMarkerColor(ColorNearEtaFit);
    gNEtaJt3->SetLineColor(ColorNearEtaFit);
 
    gNAvePhi1->SetMarkerStyle(MarkerNearPhi);
    gNAvePhi1->SetMarkerSize(MarkerSize);
    gNAvePhi1->SetMarkerColor(ColorNearPhi);
    gNAvePhi1->SetLineColor(ColorNearPhi);
    gAAvePhi1->SetMarkerStyle(MarkerAwayPhi);
    gAAvePhi1->SetMarkerSize(MarkerSize);
    gAAvePhi1->SetMarkerColor(ColorAwayPhi);
    gAAvePhi1->SetLineColor(ColorAwayPhi);
    gMAvePhi1->SetMarkerStyle(MarkerMin);
    gMAvePhi1->SetMarkerSize(MarkerSize);
    gMAvePhi1->SetMarkerColor(ColorMin);
    gMAvePhi1->SetLineColor(ColorMin);
 
    gNAvePhi2->SetMarkerStyle(MarkerNearPhi);
    gNAvePhi2->SetMarkerSize(MarkerSize);
    gNAvePhi2->SetMarkerColor(ColorNearPhi);
    gNAvePhi2->SetLineColor(ColorNearPhi);
    gAAvePhi2->SetMarkerStyle(MarkerAwayPhi);
    gAAvePhi2->SetMarkerSize(MarkerSize);
    gAAvePhi2->SetMarkerColor(ColorAwayPhi);
    gAAvePhi2->SetLineColor(ColorAwayPhi);
    gMAvePhi2->SetMarkerStyle(MarkerMin);
    gMAvePhi2->SetMarkerSize(MarkerSize);
    gMAvePhi2->SetMarkerColor(ColorMin);
    gMAvePhi2->SetLineColor(ColorMin);


    gNAvePhiFit1->SetMarkerStyle(MarkerNearPhiFit);
    gNAvePhiFit1->SetMarkerSize(MarkerSize);
    gNAvePhiFit1->SetMarkerColor(ColorNearPhiFit);
    gNAvePhiFit1->SetLineColor(ColorNearPhiFit);
    gAAvePhiFit1->SetMarkerStyle(MarkerAwayPhiFit);
    gAAvePhiFit1->SetMarkerSize(MarkerSize);
    gAAvePhiFit1->SetMarkerColor(ColorAwayPhiFit);
    gAAvePhiFit1->SetLineColor(ColorAwayPhiFit);
    gMAvePhiFit1->SetMarkerStyle(MarkerMinFit);
    gMAvePhiFit1->SetMarkerSize(MarkerSize);
    gMAvePhiFit1->SetMarkerColor(ColorMinFit);
    gMAvePhiFit1->SetLineColor(ColorMinFit);
  
    gNAvePhiFit2->SetMarkerStyle(MarkerNearPhiFit);
    gNAvePhiFit2->SetMarkerSize(MarkerSize);
    gNAvePhiFit2->SetMarkerColor(ColorNearPhiFit);
    gNAvePhiFit2->SetLineColor(ColorNearPhiFit);
    gAAvePhiFit2->SetMarkerStyle(MarkerAwayPhiFit);
    gAAvePhiFit2->SetMarkerSize(MarkerSize);
    gAAvePhiFit2->SetMarkerColor(ColorAwayPhiFit);
    gAAvePhiFit2->SetLineColor(ColorAwayPhiFit);
    gMAvePhiFit2->SetMarkerStyle(MarkerMinFit);
    gMAvePhiFit2->SetMarkerSize(MarkerSize);
    gMAvePhiFit2->SetMarkerColor(ColorMinFit);
    gMAvePhiFit2->SetLineColor(ColorMinFit);
   
    char outName[200];
    sprintf(outName,"%s/DrawSpectra_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d.root",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult);
    TFile *fout=new TFile(name,"recreate");
   
   
    float plotmin=100, plotmax=-100;
    plotmin2=100; plotmax2=-100;
    plotmin3=100; plotmax3=-100;
    float plotminW=100, plotmaxW=-100;
    float plotminJ=100, plotmaxJ=-100;
    float maxheight=100;
    float plotminAve=100, plotmaxAve=-100;
    for(int i=0;i<=2;i+=2){
      for(int j=0;j<nMainLoop;j++){
	if(NYieldPhiZYAM[i][j]<plotmin)plotmin=NYieldPhiZYAM[i][j];
	if(NYieldPhiZYAM[i][j]>plotmax&&NYieldPhiZYAM[i][j]<maxheight)plotmax=NYieldPhiZYAM[i][j];
	if(AYieldPhiZYAM[i][j]<plotmin)plotmin=AYieldPhiZYAM[i][j];
	if(AYieldPhiZYAM[i][j]>plotmax&&AYieldPhiZYAM[i][j]<maxheight)plotmax=AYieldPhiZYAM[i][j];
	if(NYieldEtaZYAM[i][j]<plotmin)plotmin=NYieldEtaZYAM[i][j];
	if(NYieldEtaZYAM[i][j]>plotmax&&NYieldEtaZYAM[i][j]<maxheight)plotmax=NYieldEtaZYAM[i][j];
	if(AYieldEtaZYAM[i][j]<plotmin)plotmin=AYieldEtaZYAM[i][j];
	if(AYieldEtaZYAM[i][j]>plotmax&&AYieldEtaZYAM[i][j]<maxheight)plotmax=AYieldEtaZYAM[i][j];
	if(MYieldPhi[i][j]<plotmin2)plotmin2=MYieldPhi[0][j];
	if(MYieldPhi[i][j]>plotmax2&&MYieldPhi[0][j]<maxheight)plotmax2=MYieldPhi[0][j];
	if(MYieldPhi[2][j]<plotmin3)plotmin3=MYieldPhi[2][j];
	if(MYieldPhi[2][j]>plotmax3&&MYieldPhi[2][j]<maxheight)plotmax3=MYieldPhi[2][j];
	if(DrawFit){
	  if(NYieldPhiFit[i][j]<plotmin)plotmin=NYieldPhiFit[i][j];
	  if(NYieldPhiFit[i][j]>plotmax&&NYieldPhiFit[i][j]<maxheight)plotmax=NYieldPhiFit[i][j];
	  if(AYieldPhiFit[i][j]<plotmin)plotmin=AYieldPhiFit[i][j];
	  if(AYieldPhiFit[i][j]>plotmax&&AYieldPhiFit[i][j]<maxheight)plotmax=AYieldPhiFit[i][j];
	  if(MYieldPhiFit[0][j]<plotmin2)plotmin2=MYieldPhiFit[0][j];
	  if(MYieldPhiFit[0][j]>plotmax2&&MYieldPhiFit[0][j]<maxheight)plotmax2=MYieldPhiFit[0][j];
	  if(MYieldPhiFit[2][j]<plotmin3)plotmin3=MYieldPhiFit[2][j];
	  if(MYieldPhiFit[2][j]>plotmax3&&MYieldPhiFit[2][j]<maxheight)plotmax3=MYieldPhiFit[2][j];
	}
	if(i==0){
	  if(NWidthPhi[0][j]<plotminW)plotminW=NWidthPhi[0][j];
	  if(NWidthPhi[0][j]>plotmaxW&&NWidthPhi[0][j]<5)plotmaxW=NWidthPhi[0][j];
	  if(AWidthPhi[0][j]<plotminW)plotminW=AWidthPhi[0][j];
	  if(AWidthPhi[0][j]>plotmaxW&&AWidthPhi[0][j]<6)plotmaxW=AWidthPhi[0][j];
	  if(NJtPhi[0][j]<plotminJ)plotminJ=NJtPhi[0][j];
	  if(NJtPhi[0][j]>plotmaxJ&&NJtPhi[0][j]<5)plotmaxJ=NJtPhi[0][j];
	  if(AJtPhi[0][j]<plotminJ)plotminJ=AJtPhi[0][j];
	  if(AJtPhi[0][j]>plotmaxJ&&AJtPhi[0][j]<6)plotmaxJ=AJtPhi[0][j];
	  if(NAvePhi[0][j]<plotminAve)plotminAve=NAvePhi[0][j];
	  if(NAvePhi[0][j]>plotmaxAve&&NAvePhi[0][j]<5)plotmaxW=NAvePhi[0][j];
	  if(AAvePhi[0][j]<plotminAve)plotminAve=AAvePhi[0][j];
	  if(AAvePhi[0][j]>plotmaxAve&&AAvePhi[0][j]<6)plotmaxAve=AAvePhi[0][j];
	}
      }
    }
    plotmin/=1.2-0.2;
    plotmax*=1.2+0.2;
    plotmin2/=1.5-0.2;
    plotmax2*=1.5+0.2;
    plotmin3/=1.5-0.2;
    plotmax3*=1.5+0.2;
    if(plotmax>maxheight)plotmax=maxheight;
    if(plotmax2>maxheight)plotmax2=maxheight;
    if(plotmax3>maxheight)plotmax3=maxheight;
    plotminW/1.5-0.2;
    plotmaxW*=1.5+0.2;
    plotmaxJ*=1.5+0.2;
    plotmaxAve*=1.5+1;
    plotminAve/=1.6-0.2;
    plotmin=0;
    plotmin2=0;
    plotmin3=0;
    plotminW=0;
    plotminJ=0;
    plotminAve=0;
    //plotminAve=0;
    if(plotmin<=0)plotmin=-0.00001;
    if(plotmin2<=0)plotmin2=-0.00001;
    if(plotmin3<=0)plotmin3=-0.00001;
    if(plotminW<=0)plotminW=-0.00001;
    if(plotminAve<=0)plotminAve=-0.00001;
    plotmaxW=1;

    cTYield=new TCanvas("cTYield","P_{T}^{Trig} Dependence of Yields",800,600);
    SetMargins1D(cTYield);
    sprintf(name,"Yields %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYield=new TH2F("hTYield",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYield->SetTitle("");
    hTYield->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYield->GetXaxis()->SetTitleColor(1);
    hTYield->GetYaxis()->SetTitle("#frac{N_{ch}}{N_{Trigger}}    ");
    SetTitles1D(hTYield);
    hTYield->Draw();
    //cTYield->SetLogy(1);
    gNYieldPhi->Draw("p");
    gAYieldPhi->Draw("p");
    if(DrawFit)gNYieldPhiFit->Draw("p");
    if(DrawFit)gAYieldPhiFit->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_YieldMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYield->SaveAs(outName);

    cTYieldPt=new TCanvas("cTYieldPt","P_{T}^{Trig} Dependence of Peak P_{T}",800,600);
    SetMargins1D(cTYieldPt);
    sprintf(name,"Peak P_{T} %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYieldPt=new TH2F("hTYieldPt",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYieldPt->SetTitle("");
    hTYieldPt->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYieldPt->GetXaxis()->SetTitleColor(1);
    hTYieldPt->GetYaxis()->SetTitle("#frac{p_{T}}{N_{Trigger}}    ");
    SetTitles1D(hTYieldPt);
    hTYieldPt->Draw();
    // cTYieldPt->SetLogy(1);
    gNYieldPhiPt->Draw("p");
    gAYieldPhiPt->Draw("p");
    if(DrawFit)gNYieldPhiFitPt->Draw("p");
    if(DrawFit)gAYieldPhiFitPt->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_PtMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYieldPt->SaveAs(outName);

    cTYield2=new TCanvas("cTYield2","P_{T}^{Trig} Dependence of Yields",800,600);
    SetMargins1D(cTYield2);
    sprintf(name,"Yields %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYield2=new TH2F("hTYield2",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYield2->SetTitle("");
    hTYield2->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYield2->GetXaxis()->SetTitleColor(1);
    hTYield2->GetYaxis()->SetTitle("#frac{N_{ch}}{N_{Trigger}}    ");
    SetTitles1D(hTYield2);
    hTYield2->Draw();
    //cTYield2->SetLogy(1);
    gNYieldPhi2->Draw("p");
    gAYieldPhi2->Draw("p");
    if(DrawFit)gNYieldPhiFit2->Draw("p");
    if(DrawFit)gAYieldPhiFit2->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_YieldAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYield2->SaveAs(outName);

    cTYieldPt2=new TCanvas("cTYieldPt2","P_{T}^{Trig} Dependence of Peak P_{T}",800,600);
    SetMargins1D(cTYieldPt2);
    sprintf(name,"Peak P_{T} %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYieldPt2=new TH2F("hTYieldPt2",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYieldPt2->SetTitle("");
    hTYieldPt2->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYieldPt2->GetXaxis()->SetTitleColor(1);
    hTYieldPt2->GetYaxis()->SetTitle("#frac{p_{T}}{N_{Trigger}}    ");
    SetTitles1D(hTYieldPt2);
    hTYieldPt2->Draw();
    //cTYieldPt2->SetLogy(1);
    gNYieldPhiPt2->Draw("p");
    gAYieldPhiPt2->Draw("p");
    if(DrawFit)gNYieldPhiFitPt2->Draw("p");
    if(DrawFit)gAYieldPhiFitPt2->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_PtAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYieldPt2->SaveAs(outName);

    //low
    cTYield3=new TCanvas("cTYield3","P_{T}^{Trig} Dependence of Yields",800,600);
    SetMargins1D(cTYield3);
    sprintf(name,"Yields %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYield3=new TH2F("hTYield3",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYield3->SetTitle("");
    hTYield3->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYield3->GetXaxis()->SetTitleColor(1);
    hTYield3->GetYaxis()->SetTitle("#frac{N_{ch}}{N_{Trigger}}    ");
    SetTitles1D(hTYield2);
    hTYield3->Draw();
    //cTYield2->SetLogy(1);
    gNYieldPhi3->Draw("p");
    gAYieldPhi3->Draw("p");
    if(DrawFit)gNYieldPhiFit3->Draw("p");
    if(DrawFit)gAYieldPhiFit3->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_YieldLowErr_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYield3->SaveAs(outName);

    cTYieldPt3=new TCanvas("cTYieldPt3","P_{T}^{Trig} Dependence of Peak P_{T}",800,600);
    SetMargins1D(cTYieldPt3);
    sprintf(name,"Peak P_{T} %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYieldPt3=new TH2F("hTYieldPt3",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYieldPt3->SetTitle("");
    hTYieldPt3->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYieldPt3->GetXaxis()->SetTitleColor(1);
    hTYieldPt3->GetYaxis()->SetTitle("#frac{p_{T}}{N_{Trigger}}    ");
    SetTitles1D(hTYieldPt3);
    hTYieldPt3->Draw();
    //cTYieldPt2->SetLogy(1);
    gNYieldPhiPt3->Draw("p");
    gAYieldPhiPt3->Draw("p");
    if(DrawFit)gNYieldPhiFitPt3->Draw("p");
    if(DrawFit)gAYieldPhiFitPt3->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_PtLowErr_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYieldPt3->SaveAs(outName);



    cTYieldE=new TCanvas("cTYieldE","P_{T}^{Trig} Dependence of Eta Yields",800,600);
    SetMargins1D(cTYieldE);
    sprintf(name,"Yields %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYield=new TH2F("hTYieldE",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYieldE->SetTitle("");
    hTYieldE->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYieldE->GetXaxis()->SetTitleColor(1);
    hTYieldE->GetYaxis()->SetTitle("#frac{N_{ch}}{N_{Trigger}}    ");
    SetTitles1D(hTYieldE);
    hTYieldE->Draw();
    //cTYield->SetLogy(1);
    gNYieldEta->Draw("p");
    gAYieldEta->Draw("p");
    if(DrawFit)gNYieldEtaFit->Draw("p");
    // if(DrawFit)gAYieldEtaFit->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    //if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayEtaFit,MarkerAwayEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_YieldEtaMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYieldE->SaveAs(outName);

    cTYieldPtE=new TCanvas("cTYieldPtE","P_{T}^{Trig} Dependence of Peak P_{T} Eta",800,600);
    SetMargins1D(cTYieldPtE);
    sprintf(name,"Peak P_{T} %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYieldPtE=new TH2F("hTYieldPtE",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYieldPtE->SetTitle("");
    hTYieldPtE->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYieldPtE->GetXaxis()->SetTitleColor(1);
    hTYieldPtE->GetYaxis()->SetTitle("#frac{p_{T}}{N_{Trigger}}    ");
    SetTitles1D(hTYieldPtE);
    hTYieldPtE->Draw();
    // cTYieldPt->SetLogy(1);
    gNYieldEtaPt->Draw("p");
    gAYieldEtaPt->Draw("p");
    if(DrawFit)gNYieldEtaFitPt->Draw("p");
    //if(DrawFit)gAYieldEtaFitPt->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    //if(DrawFit)keySymbol(.65,.7,"Away Fit",ColorAwayEtaFit,MarkerAwayEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_PtEtaMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYieldPt->SaveAs(outName);

  
    cTYieldC=new TCanvas("cTYieldC","P_{T}^{Trig} Dependence of Yields",800,600);
    SetMargins1D(cTYieldC);
    sprintf(name,"Yields %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYieldC=new TH2F("hTYieldC",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYieldC->SetTitle("");
    hTYieldC->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYieldC->GetXaxis()->SetTitleColor(1);
    hTYieldC->GetYaxis()->SetTitle("#frac{N_{ch}}{N_{Trigger}}    ");
    SetTitles1D(hTYieldC);
    hTYieldC->Draw();
    //cTYield->SetLogy(1);
    gNYieldPhi->Draw("p");
    gAYieldPhi->Draw("p");
    if(DrawFit)gNYieldPhiFit->Draw("p");
    if(DrawFit)gAYieldPhiFit->Draw("p");
    gNYieldEta->Draw("p");
    gAYieldEta->Draw("p");
    if(DrawFit)gNYieldEtaFit->Draw("p");
    // if(DrawFit)gAYieldEtaFit->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.65,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.6,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_YieldCombMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYieldC->SaveAs(outName);

    cTYieldPtC=new TCanvas("cTYieldPtC","P_{T}^{Trig} Dependence of Yields",800,600);
    SetMargins1D(cTYieldPtC);
    sprintf(name,"Peak p_{T} %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYieldPtC=new TH2F("hTYieldPtC",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYieldPtC->SetTitle("");
    hTYieldPtC->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYieldPtC->GetXaxis()->SetTitleColor(1);
    hTYieldPtC->GetYaxis()->SetTitle("#frac{p_{T}^{ch}}{N_{Trigger}}    ");
    SetTitles1D(hTYieldPtC);
    hTYieldPtC->Draw();
    //cTYield->SetLogy(1);
    gNYieldPhiPt->Draw("p");
    gAYieldPhiPt->Draw("p");
    if(DrawFit)gNYieldPhiFitPt->Draw("p");
    if(DrawFit)gAYieldPhiFitPt->Draw("p");
    gNYieldEtaPt->Draw("p");
    gAYieldEtaPt->Draw("p");
    if(DrawFit)gNYieldEtaFitPt->Draw("p");
    // if(DrawFit)gAYieldEtaFit->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.65,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.6,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_PtCombMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYieldPtC->SaveAs(outName);

    //
    cTYieldC2=new TCanvas("cTYieldC2","P_{T}^{Trig} Dependence of Yields 2",800,600);
    SetMargins1D(cTYieldC2);
    sprintf(name,"Yields %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYieldC2=new TH2F("hTYieldC2",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYieldC2->SetTitle("");
    hTYieldC2->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYieldC2->GetXaxis()->SetTitleColor(1);
    hTYieldC2->GetYaxis()->SetTitle("#frac{N_{ch}}{N_{Trigger}}    ");
    SetTitles1D(hTYieldC2);
    hTYieldC2->Draw();
    //cTYield->SetLogy(1);
    gNYieldPhi2->Draw("p");
    gAYieldPhi2->Draw("p");
    if(DrawFit)gNYieldPhiFit2->Draw("p");
    if(DrawFit)gAYieldPhiFit2->Draw("p");
    gNYieldEta2->Draw("p");
    gAYieldEta2->Draw("p");
    if(DrawFit)gNYieldEtaFit2->Draw("p");
    // if(DrawFit)gAYieldEtaFit->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.65,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.6,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_YieldCombAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYieldC2->SaveAs(outName);

    cTYieldPtC2=new TCanvas("cTYieldPtC2","P_{T}^{Trig} Dependence of Yields",800,600);
    SetMargins1D(cTYieldPtC2);
    sprintf(name,"Peak p_{T} %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTYieldPtC2=new TH2F("hTYieldPtC2",name,10,0,TPt2+1,10,plotmin,plotmax);
    if(NoTitle)hTYieldPtC2->SetTitle("");
    hTYieldPtC2->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTYieldPtC2->GetXaxis()->SetTitleColor(1);
    hTYieldPtC2->GetYaxis()->SetTitle("#frac{p_{T}^{ch}}{N_{Trigger}}    ");
    SetTitles1D(hTYieldPtC2);
    hTYieldPtC2->Draw();
    //cTYield->SetLogy(1);
    gNYieldPhiPt2->Draw("p");
    gAYieldPhiPt2->Draw("p");
    if(DrawFit)gNYieldPhiFitPt2->Draw("p");
    if(DrawFit)gAYieldPhiFitPt2->Draw("p");
    gNYieldEtaPt2->Draw("p");
    gAYieldEtaPt2->Draw("p");
    if(DrawFit)gNYieldEtaFitPt2->Draw("p");
    // if(DrawFit)gAYieldEtaFit->Draw("p");
    sprintf(name,"Near #Delta#phi |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.65,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.75,"Near #Delta#phi Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.7,"Away #Delta#phi Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.65,"Near #Delta#eta",ColorNearEta,MarkerNearEta,0.06,1.2*MarkerSize);
    keySymbol(.65,.6,"Away #Delta#eta",ColorAwayEta,MarkerAwayEta,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.65,.55,"Near #Delta#eta Fit",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_PtCombAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTYieldPtC2->SaveAs(outName);


    //
    cTMin=new TCanvas("cTMin","P_{T}^{Trig} Dependence of Underlying Event",800,600);
    SetMargins1D(cTMin);
    sprintf(name,"Underlying Event %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTMin=new TH2F("hTMin",name,10,0,TPt2+1,10,plotmin2,plotmax2);
    if(NoTitle)hTMin->SetTitle("");
    hTMin->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTMin->GetXaxis()->SetTitleColor(1);
    hTMin->GetYaxis()->SetTitle("#frac{p_{T}^{ch}}{N_{Trigger}}    ");
    SetTitles1D(hTMin);
    hTMin->Draw();
    // cTMin->SetLogy(1);
    gMYieldPhi->Draw("p");
    if(DrawFit)gMYieldPhiFit->Draw("p");
    sprintf(name,"ZYAM |#Delta#phi#pm%1.1f|<%1.1f",ZYAMCent,ZYAMWidth);
    keySymbol(.5,.85,name,ColorMin,MarkerMin,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.8,"Fit",ColorMinFit,MarkerMinFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_MinMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTMin->SaveAs(outName);

    cTMinPt=new TCanvas("cTMinPt","P_{T}^{Trig} Dependence of Underlying Event P_{T}",800,600);
    SetMargins1D(cTMinPt);
    sprintf(name,"Min P_{T} %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTMinPt=new TH2F("hTMinPt",name,10,0,TPt2+1,10,plotmin3,plotmax3);
    if(NoTitle)hTMinPt->SetTitle("");
    hTMinPt->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTMinPt->GetXaxis()->SetTitleColor(1);
    hTMinPt->GetYaxis()->SetTitle("#frac{p_{T}}{N_{Trigger}}    ");
    SetTitles1D(hTMinPt);
    hTMinPt->Draw();
    //cTMinPt->SetLogy(1);
    gMYieldPhiPt->Draw("p");
    if(DrawFit)gMYieldPhiFitPt->Draw("p");
    sprintf(name,"ZYAM |#Delta#phi#pm%1.1f|<%1.1f",ZYAMCent,ZYAMWidth);
    keySymbol(.5,.85,name,ColorMin,MarkerMin,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.8,"Fit",ColorMinFit,MarkerMinFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_MinPtMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTMinPt->SaveAs(outName);

    cTMin2=new TCanvas("cTMin2","P_{T}^{Trig} Dependence of Underlying Event",800,600);
    SetMargins1D(cTMin2);
    sprintf(name,"Underlying Event %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTMin2=new TH2F("hTMin2",name,10,0,TPt2+1,10,plotmin2,plotmax2);
    if(NoTitle)hTMin2->SetTitle("");
    hTMin2->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTMin2->GetXaxis()->SetTitleColor(1);
    hTMin2->GetYaxis()->SetTitle("#frac{1}{N_{Trigger}}    ");
    SetTitles1D(hTMin2);
    hTMin2->Draw();
    //cTMin2->SetLogy(1);
    gMYieldPhi2->Draw("p");
    if(DrawFit)gMYieldPhiFit2->Draw("p");
    sprintf(name,"ZYAM |#Delta#phi#pm%1.1f|<%1.1f",ZYAMCent,ZYAMWidth);
     keySymbol(.5,.85,name,ColorMin,MarkerMin,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.8,"Fit",ColorMinFit,MarkerMinFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_MinAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTMin2->SaveAs(outName);

    cTMinPt2=new TCanvas("cTMinPt2","P_{T}^{Trig} Dependence of Underlying Event P_{T}",800,600);
    SetMargins1D(cTMinPt2);
    sprintf(name,"Min P_{T} %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTMinPt2=new TH2F("hTMinPt2",name,10,0,TPt2+1,10,plotmin3,plotmax3);
    if(NoTitle)hTMinPt2->SetTitle("");
    hTMinPt2->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTMinPt2->GetXaxis()->SetTitleColor(1);
    hTMinPt2->GetYaxis()->SetTitle("#frac{p_{T}}{N_{Trigger}}    ");
    SetTitles1D(hTMinPt2);
    hTMinPt2->Draw();
    //cTMinPt2->SetLogy(1);
    gMYieldPhiPt2->Draw("p");
    if(DrawFit)gMYieldPhiFitPt2->Draw("p");
    sprintf(name,"ZYAM |#Delta#phi#pm%1.1f|<%1.1f",ZYAMCent,ZYAMWidth);
    keySymbol(.5,.85,name,ColorMin,MarkerMin,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.8,"Fit",ColorMinFit,MarkerMinFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_MinPtAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffCorr,APtTPtMult,FitTit[DrawFit],filetype);
    cTMinPt2->SaveAs(outName); 

    cTMin3=new TCanvas("cTMin3","P_{T}^{Trig} Dependence of Underlying Event",800,600);
    SetMargins1D(cTMin3);
    sprintf(name,"Underlying Event %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTMin3=new TH2F("hTMin3",name,10,0,TPt2+1,10,plotmin2,plotmax2);
    if(NoTitle)hTMin3->SetTitle("");
    hTMin3->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTMin3->GetXaxis()->SetTitleColor(1);
    hTMin3->GetYaxis()->SetTitle("#frac{N}{N_{Trigger}}    ");
    SetTitles1D(hTMin3);
    hTMin3->Draw();
    //cTMin2->SetLogy(1);
    gMYieldPhi3->Draw("p");
    if(DrawFit)gMYieldPhiFit3->Draw("p");
    sprintf(name,"ZYAM |#Delta#phi#pm%1.1f|<%1.1f",ZYAMCent,ZYAMWidth);
    keySymbol(.5,.85,name,ColorMin,MarkerMin,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.8,"Fit",ColorMinFit,MarkerMinFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_MinLowErr_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffCorr,APtTPtMult,FitTit[DrawFit],filetype);
    cTMin3->SaveAs(outName);

    cTMinPt3=new TCanvas("cTMinPt3","P_{T}^{Trig} Dependence of Underlying Event P_{T}",800,600);
    SetMargins1D(cTMinPt3);
    sprintf(name,"Min P_{T} %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTMinPt3=new TH2F("hTMinPt3",name,10,0,TPt2+1,10,plotmin3,plotmax3);
    if(NoTitle)hTMinPt3->SetTitle("");
    hTMinPt3->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTMinPt3->GetXaxis()->SetTitleColor(1);
    hTMinPt3->GetYaxis()->SetTitle("#frac{p_{T}}{N_{Trigger}}    ");
    SetTitles1D(hTMinPt3);
    hTMinPt3->Draw();
    //cTMinPt2->SetLogy(1);
    gMYieldPhiPt3->Draw("p");
    if(DrawFit)gMYieldPhiFitPt3->Draw("p");
    sprintf(name,"ZYAM |#Delta#phi#pm%1.1f|<%1.1f",ZYAMCent,ZYAMWidth);
    keySymbol(.5,.85,name,ColorMin,MarkerMin,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.8,"Fit",ColorMinFit,MarkerMinFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_MinPtLowErr_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffCorr,APtTPtMult,filetype);
    cTMinPt3->SaveAs(outName); 
    cWidth1=new TCanvas("cWidth1","Width MidPoint",800,600);
    SetMargins1D(cWidth1);
    sprintf(name,"Peak Width %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hWidth1=new TH2F("hWidth1",name,10,0,TPt2+1,10,plotminW,plotmaxW);
    if(NoTitle)hWidth1->SetTitle("");
    hWidth1->GetXaxis()->SetTitle("p_{T}^{Trig} (GeV/c) ");
    hWidth1->GetXaxis()->SetTitleColor(1);
    hWidth1->GetYaxis()->SetTitle("#sigma (radians)   ");
    SetTitles1D(hWidth1);
    hWidth1->Draw();
    gNPhiWidth1->Draw("p");
    gAPhiWidth1->Draw("p");
    gNEtaWidth1->Draw("p");
    keySymbol(.65,.85,"Near #Delta#phi",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_WidthTrigMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,filetype);
    cWidth1->SaveAs(outName);

    cWidth2=new TCanvas("cWidth2","Width Average",800,600);
    SetMargins1D(cWidth2);
    sprintf(name,"Peak Width %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hWidth2=new TH2F("hWidth2",name,10,0,TPt2+1,10,plotminW,plotmaxW);
    if(NoTitle)hWidth2->SetTitle("");
    hWidth2->GetXaxis()->SetTitle("p_{T}^{Trig} (GeV/c) ");
    hWidth2->GetXaxis()->SetTitleColor(1);
    hWidth2->GetYaxis()->SetTitle("#sigma (radians)   ");
    SetTitles1D(hWidth2);
    hWidth2->Draw();
    gNPhiWidth2->Draw("p");
    gAPhiWidth2->Draw("p");
    gNEtaWidth2->Draw("p");
    keySymbol(.65,.85,"Near #Delta#phi",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_WidthTrigAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,filetype);
    cWidth2->SaveAs(outName);

    cWidth3=new TCanvas("cWidth3","Width Low Err",800,600);
    SetMargins1D(cWidth3);
    sprintf(name,"Peak Width %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hWidth3=new TH2F("hWidth3",name,10,0,TPt2+1,10,plotminW,plotmaxW);
    if(NoTitle)hWidth3->SetTitle("");
    hWidth3->GetXaxis()->SetTitle("p_{T}^{Trig} (GeV/c) ");
    hWidth3->GetXaxis()->SetTitleColor(1);
    hWidth3->GetYaxis()->SetTitle("#sigma (radians)   ");
    SetTitles1D(hWidth3);
    hWidth3->Draw();
    gNPhiWidth3->Draw("p");
    gAPhiWidth3->Draw("p");
    gNEtaWidth3->Draw("p");
    keySymbol(.65,.85,"Near #Delta#phi",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.75,"Near #Delta#eta",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_WidthTrigLowError_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,filetype);
    cWidth3->SaveAs(outName);

    cJt1=new TCanvas("cJt1","Jt MidPoint",800,600);
    SetMargins1D(cJt1);
    sprintf(name,"Peak Jt %3.1f<Pt^{Trig}<%3.1f %dM%d",TPt1,TPt2,Mult1,Mult2);
    TH2F *hJt1=new TH2F("hJt1",name,10,0,TPt2+1,10,plotminJ,plotmaxJ);
    if(NoTitle)hJt1->SetTitle("");
    hJt1->GetXaxis()->SetTitle("p_{T}^{Trig} (GeV/c) ");
    hJt1->GetXaxis()->SetTitleColor(1);
    hJt1->GetYaxis()->SetTitle("#sqrt{<j_{T}^{2}>}   ");
    SetTitles1D(hJt1);
    hJt1->Draw();
    gNPhiJt1->Draw("p");
    //  gAPhiJt1->Draw("p");
    gNEtaJt1->Draw("p");
    keySymbol(.65,.85,"Near #Delta#phi",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    // keySymbol(.65,.8,"Away #Delta#phi",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    keySymbol(.65,.8,"Near #Delta#eta",ColorNearEtaFit,MarkerNearEtaFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_JtTrigMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,filetype);
    cJt1->SaveAs(outName);

    cTAve=new TCanvas("cTAve","P_{T}^{Trig} Dependence of <p_{T}>",800,600);
    SetMargins1D(cTAve);
    sprintf(name,"<p_{T}^{Assoc}> %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTAve=new TH2F("hTAve",name,10,0,TPt2+1,10,plotminAve,plotmaxAve);
    if(NoTitle)hTAve->SetTitle("");
    hTAve->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTAve->GetXaxis()->SetTitleColor(1);
    hTAve->GetYaxis()->SetTitle("<p_{T}^{Assoc}>    ");
    SetTitles1D(hTAve);
    hTAve->Draw();
    //cTYield->SetLogy(1);
    gNAvePhi1->Draw("p");
    gAAvePhi1->Draw("p");
    gMAvePhi1->Draw("p");
    if(DrawFit)gNAvePhiFit1->Draw("p");
    if(DrawFit)gAAvePhiFit1->Draw("p");
    if(DrawFit)gMAvePhiFit1->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.50,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.50,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    keySymbol(.50,.75,"Underlying Event",ColorMin,MarkerMin,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.70,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.65,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.50,.6,"Underlying Event Fit",ColorMinFit,MarkerMinFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_AveMidPoint_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTAve->SaveAs(outName);

    cTAve2=new TCanvas("cTAve2","P_{T}^{Trig} Dependence of <p_{T}>",800,600);
    SetMargins1D(cTAve2);
    sprintf(name,"<p_{T}^{Assoc}> %3.2f<Pt^{Assoc}<%3.1f %dM%d",APt1,APt2,Mult1,Mult2);
    TH2F *hTAve2=new TH2F("hTAve2",name,10,0,TPt2+1,10,plotminAve,plotmaxAve);
    if(NoTitle)hTAve2->SetTitle("");
    hTAve2->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c) ");
    hTAve2->GetXaxis()->SetTitleColor(1);
    hTAve2->GetYaxis()->SetTitle("<p_{T}^{Assoc}>    ");
    SetTitles1D(hTAve2);
    hTAve2->Draw();
    //cTYield->SetLogy(1);
    gNAvePhi2->Draw("p");
    gAAvePhi2->Draw("p");
    gMAvePhi2->Draw("p");
    if(DrawFit)gNAvePhiFit2->Draw("p");
    if(DrawFit)gAAvePhiFit2->Draw("p");
    if(DrawFit)gMAvePhiFit2->Draw("p");
    sprintf(name,"Near |#Delta#phi|<%1.1f",NearWidthPhi);
    keySymbol(.50,.85,name,ColorNearPhi,MarkerNearPhi,0.06,1.2*MarkerSize);
    keySymbol(.50,.8,"Away",ColorAwayPhi,MarkerAwayPhi,0.06,1.2*MarkerSize);
    keySymbol(.50,.75,"Underlying Event",ColorMin,MarkerMin,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.70,"Near Fit",ColorNearPhiFit,MarkerNearPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.5,.65,"Away Fit",ColorAwayPhiFit,MarkerAwayPhiFit,0.06,1.2*MarkerSize);
    if(DrawFit)keySymbol(.50,.6,"Underlying Event Fit",ColorMinFit,MarkerMinFit,0.06,1.2*MarkerSize);
    sprintf(outName,"%s/DrawSpectra_AveAve_%3.1fPT%3.1f_%3.1fpt%3.1f_%dM%d_C%d_%d%s%s",Folder,TPt1,TPt2,APt1,APt2,Mult1,Mult2,EffMethod,APtTPtMult,FitTit[DrawFit],filetype);
    cTAve2->SaveAs(outName);
  }
}
