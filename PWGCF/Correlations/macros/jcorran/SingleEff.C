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
//Obtain single particle efficiencies from output AliAnalysisTaskDiHadron (when ran over MC events and MC histograms filled)
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

void SingleEff(int Cent){
  gROOT->Reset();
  gROOT->LoadMacro("Util9.C");
  Style(1);
  
  
  // char *EffFolder="2010-08-13/7Pythia_LHC10b5";
     char *EffFolder="2010-08-13/LHC10c6_900Pythia";

  char effName[100];
  sprintf(effName,"%s/julery_DiHadron.root",EffFolder);
  TFile *effFile=new TFile(effName);
  TList *effList=effFile->Get("julery_DiHadron");

  char name[100];
  sprintf(name,"fHistPtEff_C%d",Cent);
  TH1F *fHistPt=(TH1F*)effList->FindObject(name);
  sprintf(name,"fHistPtEff_C%d_MC",Cent);
  TH1F *fHistPtMC=(TH1F*)effList->FindObject(name);

  TH1F *fHistEff=(TH1F*)fHistPt->Clone();
  fHistEff->SetName("fHistEff");
  fHistEff->GetYaxis()->SetTitle("Efficiency+Contamination");
  fHistEff->SetTitle("");
  fHistEff->Divide(fHistPtMC);

  c1=new TCanvas("c1","",800,600);
  fHistEff->Draw();
  fHistEff->SetMaximum(1);
  fHistEff->SetMinimum(0);
  fHistEff->GetXaxis()->SetRange(0,100);
  //c1->SetLogx(1);
  //TF1 *fit1=new TF1("fit1","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5",0,100);
  //fit1->SetParameters(0.8,0.001,0.001,0.001,0.001,0.001)
  fit2=new TF1("fit2","[0]",3,100);
  fit2->SetParameter(0,0.8);
    TF1 *fit3=new TF1("fit3","[0]/[1]*exp(-0.5*pow(x/[1],2))+[2]+[3]*x",0,3);
    fit3->SetParameters(0.1,0.1,1,0.1,0.1,0.1);
    fit3->SetParLimits(0,-2,0);
    fit3->SetLineColor(4);
  fHistEff->Fit(fit3,"r");
  fit2->SetLineColor(2);
  fHistEff->Fit(fit2,"r");
  fit3->Draw("same");
  sprintf(name,"%s/SingleEff1_%d.gif",EffFolder,Cent);
  c1->SaveAs(name);
   c2=new TCanvas("c2","",800,600);
  fHistEff->GetXaxis()->SetRange(0,30);
  fHistEff->Draw();
  fit2->Draw("same");
  fit3->Draw("same");
  sprintf(name,"%s/SingleEff2_%d.gif",EffFolder,Cent);
  c2->SaveAs(name);

  cout << fit3->GetParameter(0) << ", " << fit3->GetParameter(1) << ", " << fit3->GetParameter(2) << ", " << fit3->GetParameter(3) << endl;
  cout << fit2->GetParameter(0) << endl;

}
