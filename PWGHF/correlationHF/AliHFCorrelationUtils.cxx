/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: $ */

/////////////////////////////////////////////////////////////
//
// Collections of methods used in several points of correlation analyses
//
// Author: A. Rossi, andrea.rossi@cern.ch
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>

#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>

#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

#include "AliHFCorrelationUtils.h"


using std::cout;
using std::endl;

ClassImp(AliHFCorrelationUtils)

//________________________________________________________________________________________________
TH1D* AliHFCorrelationUtils::ReflectHisto(TH1D *h,Double_t scale){
  
  TH1D *h2=new TH1D(Form("%sReflected",h->GetName()),Form("%sReflected",h->GetName()),h->GetNbinsX()/2.,0.,TMath::Pi());
  for(Int_t j=1;j<=h->GetNbinsX();j++){
    Double_t x=h->GetBinCenter(j);
    Double_t y0=h->GetBinContent(j);
    Double_t ey0=h->GetBinError(j);
    Int_t j2;
    if(x>0&&x<TMath::Pi())j2=h2->FindBin(x);
    else if(x<0)j2=h2->FindBin(-1.*x);
    else if(x>TMath::Pi())j2=h2->FindBin(2.*TMath::Pi()-x);
    else {
      printf("Point %d excluded \n",j);
      continue;
    }
    Double_t y=h2->GetBinContent(j2);
    Double_t ey=h2->GetBinError(j2);
    h2->SetBinContent(j2,(y+y0));
    h2->SetBinError(j2,TMath::Sqrt(ey0*ey0+ey*ey));  
  }
  h2->Scale(scale);


  return h2;
}

//________________________________________________________________________________________________
TH1D* AliHFCorrelationUtils::DuplicateHistoTo2piRange(TH1D *h,Double_t scale){
  if(h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(h->GetNbinsX())-h->GetBinLowEdge(1)>1.01*TMath::Pi()){
    Printf("AliHFCorrelationUtils: DuplicateHistoTo2PiRange works for histos with x range between 0 and pi");
    return 0x0;
  }

  if(h->GetBinLowEdge(h->GetNbinsX())>1.01*TMath::Pi()){
    Printf("AliHFCorrelationUtils: DuplicateHistoTo2PiRange works for histos with x range between 0 and pi");
    return 0x0;
  }
  if(h->GetBinLowEdge(1)<-0.01){
    Printf("AliHFCorrelationUtils: DuplicateHistoTo2PiRange works for histos with x range between 0 and pi");
    return 0x0;
  }


  TH1D *h2=new TH1D(Form("%sTwoPiRange",h->GetName()),Form("%sTwoPiRange",h->GetName()),h->GetNbinsX()*2.,-TMath::Pi()/2.,1.5*TMath::Pi());
  for(Int_t j=1;j<=h2->GetNbinsX();j++){
    Double_t x=h2->GetBinCenter(j);
    Int_t j2;
    if(x>0&&x<TMath::Pi())j2=h->FindBin(x);
    else if(x<0)j2=h->FindBin(-1.*x);
    else if(x>TMath::Pi())j2=h->FindBin(2.*TMath::Pi()-x);
    else {
      printf("Point %d excluded \n",j);
      continue;
    }
    Double_t y=h->GetBinContent(j2);
    Double_t ey=h->GetBinError(j2);
    h2->SetBinContent(j,y*scale);
    h2->SetBinError(j,ey*scale);
  }
  
  return h2;
}

//________________________________________________________________________________________________
void AliHFCorrelationUtils::GetMCClosureModulation(Double_t ptD, Double_t ptTrmin, Double_t ptTrmax, Double_t mod[]) {
printf("prova\n");
  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9936;
        mod[1]=0.9979;
        mod[2]=1.;
        mod[3]=1.0015;
        mod[4]=1.0013;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9979;
        mod[1]=0.9983;
        mod[2]=0.999;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9809;
        mod[1]=0.9967;
        mod[2]=1.0035;
        mod[3]=1.0045;
        mod[4]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9860;
        mod[1]=0.9963;
        mod[2]=1.;
        mod[3]=1.004;
        mod[4]=1.003;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.971;
        mod[1]=0.9915;
        mod[2]=1.0073;
        mod[3]=1.0114;
        mod[4]=1.0105;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.971;
        mod[1]=0.9915;
        mod[2]=1.0073;
        mod[3]=1.0114;
        mod[4]=1.0105;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9561;
        mod[1]=1.;
        mod[2]=1.002;
        mod[3]=1.002;
        mod[4]=1.0012;
      }                        
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9962;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9871;
        mod[1]=1.0022;
        mod[2]=1.0073;
        mod[3]=1.0061;
        mod[4]=1.0033;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9919;
        mod[1]=0.9993;
        mod[2]=1.0033;
        mod[3]=1.0042;
        mod[4]=1.0028;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9814;
        mod[1]=1.001;
        mod[2]=1.013;
        mod[3]=1.011;
        mod[4]=1.0075;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9814;
        mod[1]=1.001;
        mod[2]=1.013;
        mod[3]=1.011;
        mod[4]=1.0075;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9659;
        mod[1]=1.0254;
        mod[2]=1.0335;
        mod[3]=1.0141;
        mod[4]=1.;
      }  
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9922;
        mod[1]=1.0089;
        mod[2]=1.0072;
        mod[3]=1.0035;
        mod[4]=1.001;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.996;
        mod[1]=1.0037;
        mod[2]=1.0033;
        mod[3]=1.0024;
        mod[4]=1.002;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.990;
        mod[1]=1.008;
        mod[2]=1.013;
        mod[3]=1.009;
        mod[4]=1.;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.990;
        mod[1]=1.008;
        mod[2]=1.013;
        mod[3]=1.009;
        mod[4]=1.;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9787;
        mod[1]=1.0378;
        mod[2]=1.0262;
        mod[3]=1.004;
        mod[4]=1.;
      }  
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
      }  
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
  }
  



}
