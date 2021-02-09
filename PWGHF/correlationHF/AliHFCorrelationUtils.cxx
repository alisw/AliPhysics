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
void AliHFCorrelationUtils::GetMCClosureModulation(Double_t ptD, Double_t ptTrmin, Double_t ptTrmax, Double_t mod[], Int_t system, Int_t centbin) {

printf("Getting MC closure modulations -> Choosen system = %d\n",system);

if(system==0) { //start system 0 (usually pp2010 - if so, the modulations are at 0 since we didnt' apply a correction to data)
  
  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
    break; 

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
    break; 


    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
    break;    

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;

    }

} //end system 0

/***********************************************/
if(system==1) { //start system 1 (usually pPb2013 - if so, the modulations are at 0 since we didnt' apply a correction to data)
  
  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
    break; 

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
    break; 


    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
    break;    

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;

  }

} //end system 1

/***********************************************/
if(system==2 && centbin == 0) { //start system 2, 0-100%

  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9954;
        mod[1]=0.9984;
        mod[2]=1.;
        mod[3]=1.0014;
        mod[4]=1.009;
        mod[5]=1.0010;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9984;
        mod[1]=0.9986;
        mod[2]=0.9993;
        mod[3]=1.0004;
        mod[4]=1.0003;
        mod[5]=1.0005;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9851;
        mod[1]=0.9977;
        mod[2]=1.0028;
        mod[3]=1.0052;
        mod[4]=1.0032;
        mod[5]=1.0028;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9889;
        mod[1]=0.9974;
        mod[2]=1.0005;
        mod[3]=1.0035;
        mod[4]=1.0024;
        mod[5]=1.0028;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9763;
        mod[1]=0.9927;
        mod[2]=1.0063;
        mod[3]=1.0088;
        mod[4]=1.0065;
        mod[5]=1.0025;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9763;
        mod[1]=0.9927;
        mod[2]=1.0063;
        mod[3]=1.0088;
        mod[4]=1.0065;
        mod[5]=1.0025;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9663;
        mod[1]=1.0093;
        mod[2]=1.0179;
        mod[3]=1.0151;
        mod[4]=1.0040;
        mod[5]=1.0035;
      }                        
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9975;
        mod[1]=1.0004;
        mod[2]=1.0019;
        mod[3]=1.0015;
        mod[4]=1.0020;
        mod[5]=0.9996;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.0004;
        mod[1]=0.9999;
        mod[2]=1.0007;
        mod[3]=1.0007;
        mod[4]=1.0007;
        mod[5]=0.9994;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9892;
        mod[1]=1.0020;
        mod[2]=1.0063;
        mod[3]=1.0051;
        mod[4]=1.0027;
        mod[5]=1.0002;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9932;
        mod[1]=0.9995;
        mod[2]=1.0035;
        mod[3]=1.0035;
        mod[4]=1.0025;
        mod[5]=1.0009;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9837;
        mod[1]=1.0020;
        mod[2]=1.0102;
        mod[3]=1.0089;
        mod[4]=1.0069;
        mod[5]=1.0007;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9837;
        mod[1]=1.0020;
        mod[2]=1.0102;
        mod[3]=1.0089;
        mod[4]=1.0069;
        mod[5]=1.0007;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9690;
        mod[1]=1.0200;
        mod[2]=1.0265;
        mod[3]=1.0138;
        mod[4]=0.9968;
        mod[5]=0.9920;
      }  
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9944;
        mod[1]=1.0067;
        mod[2]=1.0056;
        mod[3]=1.0019;
        mod[4]=1.0006;
        mod[5]=0.9992;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9970;
        mod[1]=1.0029;
        mod[2]=1.0029;
        mod[3]=1.0016;
        mod[4]=1.0018;
        mod[5]=0.9991;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9928;
        mod[1]=1.0052;
        mod[2]=1.0100;
        mod[3]=1.0058;
        mod[4]=0.9992;
        mod[5]=0.9998;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9928;
        mod[1]=1.0052;
        mod[2]=1.0100;
        mod[3]=1.0058;
        mod[4]=0.9992;
        mod[5]=0.9998;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9820;
        mod[1]=1.0305;
        mod[2]=1.0192;
        mod[3]=0.9992;
        mod[4]=0.9951;
        mod[5]=1.0007;
      }  
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9983;
        mod[1]=1.0036;
        mod[2]=0.9994;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }  
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 2

/***********************************************/
if(system==2 && centbin == 1) { //start system 2, 0-20%

  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9960;
        mod[1]=0.9984;
        mod[2]=1.;
        mod[3]=1.0015;
        mod[4]=1.0006;
        mod[5]=1.0006;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9989;
        mod[1]=0.9986;
        mod[2]=0.9993;
        mod[3]=1.0005;
        mod[4]=1.0003;
        mod[5]=1.0003;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9864;
        mod[1]=0.9979;
        mod[2]=1.0026;
        mod[3]=1.0045;
        mod[4]=1.0022;
        mod[5]=1.0018;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9981;
        mod[1]=1.0001;
        mod[2]=1.0018;
        mod[3]=1.0017;
        mod[4]=1.0008;
        mod[5]=0.9996;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.0007;
        mod[1]=0.9997;
        mod[2]=1.0009;
        mod[3]=1.0006;
        mod[4]=1.0006;
        mod[5]=0.9992;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9905;
        mod[1]=1.0016;
        mod[2]=1.0054;
        mod[3]=1.0041;
        mod[4]=1.0016;
        mod[5]=1.0010;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }      
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9948;
        mod[1]=1.0054;
        mod[2]=1.0047;
        mod[3]=1.0011;
        mod[4]=1.0008;
        mod[5]=0.9988;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 2 - 0-20%


/***********************************************/
if(system==2 && centbin == 2) { //start system 2, 20-60%

  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9934;
        mod[1]=0.9980;
        mod[2]=1.0006;
        mod[3]=1.0016;
        mod[4]=1.0017;
        mod[5]=1.0014;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9974;
        mod[1]=0.9982;
        mod[2]=0.9997;
        mod[3]=1.0003;
        mod[4]=1.0010;
        mod[5]=1.0009;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9813;
        mod[1]=0.9971;
        mod[2]=1.0039;
        mod[3]=1.0066;
        mod[4]=1.0044;
        mod[5]=1.0034;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9957;
        mod[1]=1.0004;
        mod[2]=1.0021;
        mod[3]=1.0018;
        mod[4]=1.0016;
        mod[5]=0.9973;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9867;
        mod[1]=1.0021;
        mod[2]=1.0071;
        mod[3]=1.0062;
        mod[4]=1.0050;
        mod[5]=0.9967;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9976;
        mod[1]=1.0031;
        mod[2]=1.0024;
        mod[3]=1.0008;
        mod[4]=0.9992;
        mod[5]=0.9998;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9930;
        mod[1]=1.0093;
        mod[2]=1.0078;
        mod[3]=1.0041;
        mod[4]=0.9998;
        mod[5]=0.9994;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 2 - 20-60%

/***********************************************/
if(system==2 && centbin == 3) { //start system 2, 60-100%

  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9891;
        mod[1]=0.9966;
        mod[2]=1.0004;
        mod[3]=1.0022;
        mod[4]=1.0034;
        mod[5]=1.0029;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9955;
        mod[1]=0.9973;
        mod[2]=0.9989;
        mod[3]=1.;
        mod[4]=1.0018;
        mod[5]=1.0015;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9709;
        mod[1]=0.9938;
        mod[2]=1.0054;
        mod[3]=1.0100;
        mod[4]=1.0096;
        mod[5]=1.0079;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9942;
        mod[1]=1.0008;
        mod[2]=1.0030;
        mod[3]=1.0033;
        mod[4]=1.0022;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9828;
        mod[1]=1.0036;
        mod[2]=1.0102;
        mod[3]=1.0102;
        mod[4]=1.0072;
        mod[5]=1.0009;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9977;
        mod[1]=1.0034;
        mod[2]=1.0022;
        mod[3]=1.0003;
        mod[4]=0.9992;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9929;
        mod[1]=1.0099;
        mod[2]=1.0081;
        mod[3]=1.0045;
        mod[4]=0.9993;
        mod[5]=0.9983;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      }
      else { //dummy
        mod[0]=1.;
        mod[1]=1.;
        mod[2]=1.;
        mod[3]=1.;
        mod[4]=1.;
        mod[5]=1.;
      } 
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 2 - 60-100%

/***********************************************/
if(system==3) { //start system 3, pp at 5 TeV (2017)

  Int_t ptbinD = 0;
  if(ptD>2 && ptD<3) ptbinD = 1;
  if(ptD>3 && ptD<5) ptbinD = 2;
  if(ptD>5 && ptD<8) ptbinD = 3;
  if(ptD>8 && ptD<16) ptbinD = 4;
  if(ptD>16 && ptD<24) ptbinD = 5;

  switch(ptbinD) {

    case(1): //ptD 2 to 3
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9841;
        mod[1]=0.9920;
        mod[2]=0.9975;
        mod[3]=1.0027;
        mod[4]=1.0035;
        mod[5]=1.0051;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9902;
        mod[1]=0.9938;
        mod[2]=0.9966;
        mod[3]=1.0005;
        mod[4]=1.0012;
        mod[5]=1.0034;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9683;
        mod[1]=0.9871;
        mod[2]=0.9998;
        mod[3]=1.0087;
        mod[4]=1.0103;
        mod[5]=1.0100;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9746;
        mod[1]=0.9871;
        mod[2]=0.9967;
        mod[3]=1.0044;
        mod[4]=1.0074;
        mod[5]=1.0074;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9468;
        mod[1]=0.9803;
        mod[2]=1.0041;
        mod[3]=1.0187;
        mod[4]=1.0145;
        mod[5]=1.0196;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9459;
        mod[1]=0.9820;
        mod[2]=1.0097;
        mod[3]=1.0225;
        mod[4]=1.0197;
        mod[5]=1.0190;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9435;
        mod[1]=0.9847;
        mod[2]=1.0214;
        mod[3]=1.0290;
        mod[4]=1.0306;
        mod[5]=1.0176;
      }                        
    break;


    case(2): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9877;
        mod[1]=0.9970;
        mod[2]=1.0015;
        mod[3]=1.0038;
        mod[4]=1.0033;
        mod[5]=1.0031;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9943;
        mod[1]=0.9967;
        mod[2]=0.9991;
        mod[3]=1.0014;
        mod[4]=1.0024;
        mod[5]=1.0018;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9736;
        mod[1]=0.9974;
        mod[2]=1.0076;
        mod[3]=1.0110;
        mod[4]=1.0099;
        mod[5]=1.0071;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9797;
        mod[1]=0.9951;
        mod[2]=1.0029;
        mod[3]=1.0064;
        mod[4]=1.0078;
        mod[5]=1.0063;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9627;
        mod[1]=0.9936;
        mod[2]=1.0132;
        mod[3]=1.0231;
        mod[4]=1.0178;
        mod[5]=1.0121;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9553;
        mod[1]=1.0027;
        mod[2]=1.0213;
        mod[3]=1.0252;
        mod[4]=1.0170;
        mod[5]=1.0097;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9627;
        mod[1]=0.9935;
        mod[2]=1.0132;
        mod[3]=1.0231;
        mod[4]=1.0178;
        mod[5]=1.0121;
      }                        
    break;

    case(3): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9925;
        mod[1]=1.0019;
        mod[2]=1.0045;
        mod[3]=1.0043;
        mod[4]=1.0023;
        mod[5]=1.0012;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9983;
        mod[1]=0.9991;
        mod[2]=1.0006;
        mod[3]=1.0017;
        mod[4]=1.0009;
        mod[5]=1.0009;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9818;
        mod[1]=1.0066;
        mod[2]=1.0133;
        mod[3]=1.0111;
        mod[4]=1.0068;
        mod[5]=1.0022;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9863;
        mod[1]=1.0006;
        mod[2]=1.0073;
        mod[3]=1.0069;
        mod[4]=1.0070;
        mod[5]=1.0041;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9771;
        mod[1]=1.0080;
        mod[2]=1.0206;
        mod[3]=1.0206;
        mod[4]=1.0104;
        mod[5]=1.0005;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9704;
        mod[1]=1.0184;
        mod[2]=1.0296;
        mod[3]=1.0239;
        mod[4]=1.0063;
        mod[5]=0.9973;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9583;
        mod[1]=1.0317;
        mod[2]=1.0421;
        mod[3]=1.0276;
        mod[4]=1.0005;
        mod[5]=0.9931;
      }  
    break;

    case(4): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9950;
        mod[1]=1.0051;
        mod[2]=1.0036;
        mod[3]=1.0019;
        mod[4]=1.0004;
        mod[5]=0.9992;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9987;
        mod[1]=0.9993;
        mod[2]=1.0001;
        mod[3]=1.0008;
        mod[4]=1.0005;
        mod[5]=1.0002;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9897;
        mod[1]=1.0014;
        mod[2]=1.0114;
        mod[3]=1.0052;
        mod[4]=1.0009;
        mod[5]=0.9970;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9918;
        mod[1]=1.0046;
        mod[2]=1.0075;
        mod[3]=1.0050;
        mod[4]=1.0012;
        mod[5]=0.9916;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9876;
        mod[1]=1.0174;
        mod[2]=1.0156;
        mod[3]=1.0135;
        mod[4]=1.0006;
        mod[5]=0.9923;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9847;
        mod[1]=1.0291;
        mod[2]=1.0211;
        mod[3]=1.0071;
        mod[4]=1.0016;
        mod[5]=0.9935;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9784;
        mod[1]=1.0436;
        mod[2]=1.0287;
        mod[3]=1.0022;
        mod[4]=0.9971;
        mod[5]=0.9953;
      }  
    break;

    case(5): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.0009;
        mod[1]=1.0045;
        mod[2]=1.0002;
        mod[3]=0.9960;
        mod[4]=0.9950;
        mod[5]=1.0003;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.0005;
        mod[1]=1.0025;
        mod[2]=1.0017;
        mod[3]=0.9963;
        mod[4]=0.9946;
        mod[5]=1.0009;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.0010;
        mod[1]=1.0067;
        mod[2]=1.0032;
        mod[3]=0.9953;
        mod[4]=0.9966;
        mod[5]=0.9954;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9981;
        mod[1]=1.0009;
        mod[2]=1.0048;
        mod[3]=0.9931;
        mod[4]=0.9988;
        mod[5]=1.0093;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9967;
        mod[1]=1.0126;
        mod[2]=1.0035;
        mod[3]=0.9926;
        mod[4]=0.9893;
        mod[5]=1.0056;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=1.0015;
        mod[1]=1.0144;
        mod[2]=1.0004;
        mod[3]=1.0015;
        mod[4]=0.9933;
        mod[5]=0.9978;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=1.0031;
        mod[1]=1.0135;
        mod[2]=0.9966;
        mod[3]=1.0110;
        mod[4]=0.9991;
        mod[5]=0.9860;
      }  
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 3


if(system==4 && centbin==0) { //pp at 13 TeV (2016-17-18)

  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;
  if(ptD>24 && ptD<36) ptbinD = 5;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9877;
        mod[1]=0.9970;
        mod[2]=1.0015;
        mod[3]=1.0038;
        mod[4]=1.0033;
        mod[5]=1.0031;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9943;
        mod[1]=0.9967;
        mod[2]=0.9991;
        mod[3]=1.0014;
        mod[4]=1.0024;
        mod[5]=1.0018;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9736;
        mod[1]=0.9974;
        mod[2]=1.0076;
        mod[3]=1.0110;
        mod[4]=1.0099;
        mod[5]=1.0071;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9797;
        mod[1]=0.9951;
        mod[2]=1.0029;
        mod[3]=1.0064;
        mod[4]=1.0078;
        mod[5]=1.0063;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9627;
        mod[1]=0.9936;
        mod[2]=1.0132;
        mod[3]=1.0231;
        mod[4]=1.0178;
        mod[5]=1.0121;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9553;
        mod[1]=1.0027;
        mod[2]=1.0213;
        mod[3]=1.0252;
        mod[4]=1.0170;
        mod[5]=1.0097;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9627;
        mod[1]=0.9935;
        mod[2]=1.0132;
        mod[3]=1.0231;
        mod[4]=1.0178;
        mod[5]=1.0121;
      }                        
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9925;
        mod[1]=1.0019;
        mod[2]=1.0045;
        mod[3]=1.0043;
        mod[4]=1.0023;
        mod[5]=1.0012;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9983;
        mod[1]=0.9991;
        mod[2]=1.0006;
        mod[3]=1.0017;
        mod[4]=1.0009;
        mod[5]=1.0009;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9818;
        mod[1]=1.0066;
        mod[2]=1.0133;
        mod[3]=1.0111;
        mod[4]=1.0068;
        mod[5]=1.0022;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9863;
        mod[1]=1.0006;
        mod[2]=1.0073;
        mod[3]=1.0069;
        mod[4]=1.0070;
        mod[5]=1.0041;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9771;
        mod[1]=1.0080;
        mod[2]=1.0206;
        mod[3]=1.0206;
        mod[4]=1.0104;
        mod[5]=1.0005;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9704;
        mod[1]=1.0184;
        mod[2]=1.0296;
        mod[3]=1.0239;
        mod[4]=1.0063;
        mod[5]=0.9973;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9583;
        mod[1]=1.0317;
        mod[2]=1.0421;
        mod[3]=1.0276;
        mod[4]=1.0005;
        mod[5]=0.9931;
      }  
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.9950;
        mod[1]=1.0051;
        mod[2]=1.0036;
        mod[3]=1.0019;
        mod[4]=1.0004;
        mod[5]=0.9992;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.9987;
        mod[1]=0.9993;
        mod[2]=1.0001;
        mod[3]=1.0008;
        mod[4]=1.0005;
        mod[5]=1.0002;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=0.9897;
        mod[1]=1.0014;
        mod[2]=1.0114;
        mod[3]=1.0052;
        mod[4]=1.0009;
        mod[5]=0.9970;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9918;
        mod[1]=1.0046;
        mod[2]=1.0075;
        mod[3]=1.0050;
        mod[4]=1.0012;
        mod[5]=0.9916;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9876;
        mod[1]=1.0174;
        mod[2]=1.0156;
        mod[3]=1.0135;
        mod[4]=1.0006;
        mod[5]=0.9923;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=0.9847;
        mod[1]=1.0291;
        mod[2]=1.0211;
        mod[3]=1.0071;
        mod[4]=1.0016;
        mod[5]=0.9935;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.9784;
        mod[1]=1.0436;
        mod[2]=1.0287;
        mod[3]=1.0022;
        mod[4]=0.9971;
        mod[5]=0.9953;
      }  
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.0009;
        mod[1]=1.0045;
        mod[2]=1.0002;
        mod[3]=0.9960;
        mod[4]=0.9950;
        mod[5]=1.0003;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.0005;
        mod[1]=1.0025;
        mod[2]=1.0017;
        mod[3]=0.9963;
        mod[4]=0.9946;
        mod[5]=1.0009;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.0010;
        mod[1]=1.0067;
        mod[2]=1.0032;
        mod[3]=0.9953;
        mod[4]=0.9966;
        mod[5]=0.9954;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9981;
        mod[1]=1.0009;
        mod[2]=1.0048;
        mod[3]=0.9931;
        mod[4]=0.9988;
        mod[5]=1.0093;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9967;
        mod[1]=1.0126;
        mod[2]=1.0035;
        mod[3]=0.9926;
        mod[4]=0.9893;
        mod[5]=1.0056;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=1.0015;
        mod[1]=1.0144;
        mod[2]=1.0004;
        mod[3]=1.0015;
        mod[4]=0.9933;
        mod[5]=0.9978;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=1.0031;
        mod[1]=1.0135;
        mod[2]=0.9966;
        mod[3]=1.0110;
        mod[4]=0.9991;
        mod[5]=0.9860;
      }  
    break;

    case(5): //ptD 24 to 36
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.0009;
        mod[1]=1.0045;
        mod[2]=1.0002;
        mod[3]=0.9960;
        mod[4]=0.9950;
        mod[5]=1.0003;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=1.0005;
        mod[1]=1.0025;
        mod[2]=1.0017;
        mod[3]=0.9963;
        mod[4]=0.9946;
        mod[5]=1.0009;
      }
      else if(ptTrmin==1. && ptTrmax==99.) {
        mod[0]=1.0010;
        mod[1]=1.0067;
        mod[2]=1.0032;
        mod[3]=0.9953;
        mod[4]=0.9966;
        mod[5]=0.9954;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.9981;
        mod[1]=1.0009;
        mod[2]=1.0048;
        mod[3]=0.9931;
        mod[4]=0.9988;
        mod[5]=1.0093;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.9967;
        mod[1]=1.0126;
        mod[2]=1.0035;
        mod[3]=0.9926;
        mod[4]=0.9893;
        mod[5]=1.0056;
      }
      else if(ptTrmin==2. && ptTrmax==99.) {
        mod[0]=1.0015;
        mod[1]=1.0144;
        mod[2]=1.0004;
        mod[3]=1.0015;
        mod[4]=0.9933;
        mod[5]=0.9978;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=1.0031;
        mod[1]=1.0135;
        mod[2]=0.9966;
        mod[3]=1.0110;
        mod[4]=0.9991;
        mod[5]=0.9860;
      }  
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 4 - 0-100%

if(system==4 && centbin==1) { //pp at 13 TeV (2016-17-18), 0-0.1% V0M

  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.99239;
        mod[1]=0.99838;
        mod[2]=1.00185;
        mod[3]=1.00245;
        mod[4]=1.00207;
        mod[5]=1.00173;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99561;
        mod[1]=0.99783;
        mod[2]=0.99992;
        mod[3]=1.00084;
        mod[4]=1.00147;
        mod[5]=1.00148;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.98858;
        mod[1]=0.99767;
        mod[2]=1.00245;
        mod[3]=1.00411;
        mod[4]=1.00295;
        mod[5]=1.00262;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.98363;
        mod[1]=0.99802;
        mod[2]=1.00825;
        mod[3]=1.00569;
        mod[4]=1.00614;
        mod[5]=1.00344;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.98448;
        mod[1]=1.00920;
        mod[2]=1.01670;
        mod[3]=1.01324;
        mod[4]=1.00077;
        mod[5]=0.99827;
      }                        
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.99226;
        mod[1]=1.00158;
        mod[2]=1.00314;
        mod[3]=1.00176;
        mod[4]=1.00291;
        mod[5]=1.00057;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99799;
        mod[1]=0.99830;
        mod[2]=1.00133;
        mod[3]=1.00091;
        mod[4]=1.00267;
        mod[5]=0.99986;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.98766;
        mod[1]=1.00218;
        mod[2]=1.00371;
        mod[3]=1.00188;
        mod[4]=1.00199;
        mod[5]=1.00177;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.98100;
        mod[1]=1.00683;
        mod[2]=1.01136;
        mod[3]=1.00582;
        mod[4]=1.00362;
        mod[5]=1.00155;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.97719;
        mod[1]=1.02189;
        mod[2]=1.00972;
        mod[3]=1.00632;
        mod[4]=1.00895;
        mod[5]=1.00279;
      }  
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.99574;
        mod[1]=1.00457;
        mod[2]=1.00100;
        mod[3]=0.99948;
        mod[4]=1.00192;
        mod[5]=0.99713;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99792;
        mod[1]=0.99964;
        mod[2]=0.99937;
        mod[3]=0.99916;
        mod[4]=1.00155;
        mod[5]=0.99787;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.99513;
        mod[1]=1.00480;
        mod[2]=1.00303;
        mod[3]=1.00075;
        mod[4]=1.00241;
        mod[5]=0.99673;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.98683;
        mod[1]=1.01060;
        mod[2]=1.00836;
        mod[3]=1.00075;
        mod[4]=1.00034;
        mod[5]=0.98986;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.99065;
        mod[1]=1.02478;
        mod[2]=1.00012;
        mod[3]=0.99814;
        mod[4]=1.00703;
        mod[5]=0.99949;
      }  
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.99963;
        mod[1]=1.00303;
        mod[2]=1.00193;
        mod[3]=1.00233;
        mod[4]=0.99972;
        mod[5]=1.00211;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99943;
        mod[1]=0.99992;
        mod[2]=1.02045;
        mod[3]=1.00260;
        mod[4]=0.99948;
        mod[5]=1.00274;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.99702;
        mod[1]=1.00315;
        mod[2]=1.00124;
        mod[3]=1.00128;
        mod[4]=1.00191;
        mod[5]=1.00257;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.99021;
        mod[1]=1.00690;
        mod[2]=1.00693;
        mod[3]=1.00185;
        mod[4]=1.00620;
        mod[5]=0.99884;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.99339;
        mod[1]=1.00812;
        mod[2]=1.00293;
        mod[3]=1.01048;
        mod[4]=0.99855;
        mod[5]=1.00686;
      }  
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 4 - 0-0.1%

if(system==4 && centbin==2) { //pp at 13 TeV (2016-17-18), 0.1-10% V0M

  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.98787;
        mod[1]=0.99573;
        mod[2]=1.00115;
        mod[3]=1.00229;
        mod[4]=1.00333;
        mod[5]=1.00247;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99300;
        mod[1]=0.99540;
        mod[2]=1.00018;
        mod[3]=1.00016;
        mod[4]=1.00205;
        mod[5]=1.00203;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.98286;
        mod[1]=0.99457;
        mod[2]=0.99995;
        mod[3]=1.00425;
        mod[4]=1.00489;
        mod[5]=1.00278;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.97291;
        mod[1]=0.99319;
        mod[2]=1.00465;
        mod[3]=1.01279;
        mod[4]=1.01091;
        mod[5]=1.00662;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.96763;
        mod[1]=1.01013;
        mod[2]=1.02067;
        mod[3]=1.01253;
        mod[4]=1.00683;
        mod[5]=1.00370;
      }                        
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.98941;
        mod[1]=1.00102;
        mod[2]=1.00415;
        mod[3]=1.00456;
        mod[4]=1.00145;
        mod[5]=1.00140;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99767;
        mod[1]=0.99797;
        mod[2]=0.99997;
        mod[3]=1.00239;
        mod[4]=1.00036;
        mod[5]=1.00109;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.98387;
        mod[1]=1.00084;
        mod[2]=1.00524;
        mod[3]=1.00736;
        mod[4]=1.00333;
        mod[5]=1.00352;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.97078;
        mod[1]=1.00302;
        mod[2]=1.02058;
        mod[3]=1.01602;
        mod[4]=1.00682;
        mod[5]=0.99869;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.95768;
        mod[1]=1.03394;
        mod[2]=1.03572;
        mod[3]=1.00724;
        mod[4]=1.00100;
        mod[5]=0.99799;
      }  
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.99352;
        mod[1]=1.00380;
        mod[2]=1.00226;
        mod[3]=1.00200;
        mod[4]=0.99896;
        mod[5]=0.99809;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99705;
        mod[1]=0.99940;
        mod[2]=0.99961;
        mod[3]=1.00260;
        mod[4]=1.00003;
        mod[5]=0.99765;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.99238;
        mod[1]=1.00223;
        mod[2]=1.00299;
        mod[3]=1.00230;
        mod[4]=0.99683;
        mod[5]=1.00076;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.98377;
        mod[1]=1.00636;
        mod[2]=1.01313;
        mod[3]=1.00585;
        mod[4]=1.00383;
        mod[5]=0.99805;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.98122;
        mod[1]=1.03577;
        mod[2]=1.01683;
        mod[3]=0.99095;
        mod[4]=0.99076;
        mod[5]=0.99200;
      }  
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.00011;
        mod[1]=1.00478;
        mod[2]=1.00001;
        mod[3]=1.00427;
        mod[4]=1.00218;
        mod[5]=0.99985;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99702;
        mod[1]=1.00325;
        mod[2]=1.00061;
        mod[3]=1.00209;
        mod[4]=0.99858;
        mod[5]=1.00444;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=1.00167;
        mod[1]=1.00182;
        mod[2]=0.99607;
        mod[3]=1.00592;
        mod[4]=1.00802;
        mod[5]=0.99947;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=1.00123;
        mod[1]=1.00316;
        mod[2]=1.00701;
        mod[3]=0.99745;
        mod[4]=1.01092;
        mod[5]=1.00749;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.98774;
        mod[1]=1.01569;
        mod[2]=1.00069;
        mod[3]=1.01980;
        mod[4]=1.00773;
        mod[5]=0.97071;
      }  
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 4 - 0.1-10%

if(system==4 && centbin==3) { //pp at 13 TeV (2016-17-18), 10-30% V0M

  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.98563;
        mod[1]=0.99515;
        mod[2]=1.00186;
        mod[3]=1.00361;
        mod[4]=1.00406;
        mod[5]=1.00347;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99249;
        mod[1]=0.99496;
        mod[2]=0.99960;
        mod[3]=1.00081;
        mod[4]=1.00194;
        mod[5]=1.00248;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.97855;
        mod[1]=0.99325;
        mod[2]=1.00204;
        mod[3]=1.00622;
        mod[4]=1.00716;
        mod[5]=1.00532;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.96369;
        mod[1]=0.99037;
        mod[2]=1.01090;
        mod[3]=1.01798;
        mod[4]=1.01588;
        mod[5]=1.01009;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.94703;
        mod[1]=1.01513;
        mod[2]=1.03407;
        mod[3]=1.02629;
        mod[4]=1.01333;
        mod[5]=1.00220;
      }                        
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.98727;
        mod[1]=1.00116;
        mod[2]=1.00533;
        mod[3]=1.00457;
        mod[4]=1.00270;
        mod[5]=1.00137;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99631;
        mod[1]=0.99784;
        mod[2]=1.00042;
        mod[3]=1.00154;
        mod[4]=1.00129;
        mod[5]=1.00112;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.97996;
        mod[1]=0.99994;
        mod[2]=1.00776;
        mod[3]=1.00827;
        mod[4]=1.00577;
        mod[5]=1.00283;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.96685;
        mod[1]=1.00490;
        mod[2]=1.02422;
        mod[3]=1.02094;
        mod[4]=1.01075;
        mod[5]=1.00288;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.94692;
        mod[1]=1.04107;
        mod[2]=1.04941;
        mod[3]=1.01859;
        mod[4]=1.00242;
        mod[5]=0.99620;
      }  
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.99291;
        mod[1]=1.00530;
        mod[2]=1.00418;
        mod[3]=1.00133;
        mod[4]=1.00006;
        mod[5]=0.99928;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99723;
        mod[1]=0.99884;
        mod[2]=1.00025;
        mod[3]=1.00100;
        mod[4]=1.00044;
        mod[5]=0.99991;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.99035;
        mod[1]=1.00390;
        mod[2]=1.00686;
        mod[3]=1.00247;
        mod[4]=1.00107;
        mod[5]=0.99956;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.98265;
        mod[1]=1.01568;
        mod[2]=1.01956;
        mod[3]=1.00738;
        mod[4]=1.00096;
        mod[5]=0.99922;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.97334;
        mod[1]=1.05114;
        mod[2]=1.02965;
        mod[3]=0.99944;
        mod[4]=0.99434;
        mod[5]=0.99375;
      }  
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.99964;
        mod[1]=1.00407;
        mod[2]=0.99909;
        mod[3]=0.99984;
        mod[4]=0.99961;
        mod[5]=0.99981;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99717;
        mod[1]=1.00074;
        mod[2]=0.99915;
        mod[3]=1.00000;
        mod[4]=0.99991;
        mod[5]=1.00099;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.99771;
        mod[1]=1.00227;
        mod[2]=0.99902;
        mod[3]=1.00115;
        mod[4]=1.00149;
        mod[5]=1.00149;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.99875;
        mod[1]=1.00743;
        mod[2]=1.00506;
        mod[3]=1.00137;
        mod[4]=1.00086;
        mod[5]=0.99765;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.99526;
        mod[1]=1.02510;
        mod[2]=1.00186;
        mod[3]=1.00402;
        mod[4]=1.00085;
        mod[5]=0.99381;
      }  
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 4 - 10-30%

if(system==4 && centbin==4) { //pp at 13 TeV (2016-17-18), 30-100% V0M

  Int_t ptbinD = 0;
  if(ptD>3 && ptD<5) ptbinD = 1;
  if(ptD>5 && ptD<8) ptbinD = 2;
  if(ptD>8 && ptD<16) ptbinD = 3;
  if(ptD>16 && ptD<24) ptbinD = 4;

  switch(ptbinD) {

    case(1): //ptD 3 to 5
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.98358;
        mod[1]=0.99444;
        mod[2]=1.00210;
        mod[3]=1.00460;
        mod[4]=1.00502;
        mod[5]=1.00436;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99166;
        mod[1]=0.99463;
        mod[2]=0.99928;
        mod[3]=1.00119;
        mod[4]=1.00227;
        mod[5]=1.00266;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.97569;
        mod[1]=0.99144;
        mod[2]=1.00211;
        mod[3]=1.00758;
        mod[4]=1.00914;
        mod[5]=1.00820;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.95435;
        mod[1]=0.98594;
        mod[2]=1.01483;
        mod[3]=1.02457;
        mod[4]=1.02213;
        mod[5]=1.01421;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.92741;
        mod[1]=1.01431;
        mod[2]=1.04607;
        mod[3]=1.03815;
        mod[4]=1.02106;
        mod[5]=1.00369;
      }                        
    break;

    case(2): //ptD 5 to 8
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.98717;
        mod[1]=1.00215;
        mod[2]=1.00695;
        mod[3]=1.00602;
        mod[4]=1.00340;
        mod[5]=1.00147;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99648;
        mod[1]=0.99833;
        mod[2]=1.00116;
        mod[3]=1.00207;
        mod[4]=1.00157;
        mod[5]=1.00097;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.97947;
        mod[1]=1.00039;
        mod[2]=1.00941;
        mod[3]=1.01080;
        mod[4]=1.00777;
        mod[5]=1.00433;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.96577;
        mod[1]=1.00533;
        mod[2]=1.02987;
        mod[3]=1.02912;
        mod[4]=1.01691;
        mod[5]=1.00488;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.93792;
        mod[1]=1.04501;
        mod[2]=1.06461;
        mod[3]=1.03129;
        mod[4]=1.00231;
        mod[5]=0.99348;
      }  
    break;

    case(3): //ptD 8 to 16
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=0.99396;
        mod[1]=1.00678;
        mod[2]=1.00542;
        mod[3]=1.00169;
        mod[4]=0.99946;
        mod[5]=0.99885;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99783;
        mod[1]=0.99947;
        mod[2]=1.00074;
        mod[3]=1.00093;
        mod[4]=0.99979;
        mod[5]=0.99990;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.99095;
        mod[1]=1.00516;
        mod[2]=1.00769;
        mod[3]=1.00346;
        mod[4]=1.00172;
        mod[5]=0.99909;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.98406;
        mod[1]=1.01627;
        mod[2]=1.02625;
        mod[3]=1.01337;
        mod[4]=1.00185;
        mod[5]=0.99803;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.97160;
        mod[1]=1.05781;
        mod[2]=1.04215;
        mod[3]=1.00227;
        mod[4]=0.99054;
        mod[5]=0.99023;
      }  
    break;

    case(4): //ptD 16 to 24
      if(ptTrmin==0.3 && ptTrmax==99.) {
        mod[0]=1.00104;
        mod[1]=1.00440;
        mod[2]=0.99820;
        mod[3]=0.99913;
        mod[4]=0.99825;
        mod[5]=0.99845;
      }
      else if(ptTrmin==0.3 && ptTrmax==1.) {
        mod[0]=0.99842;
        mod[1]=0.99969;
        mod[2]=0.99815;
        mod[3]=1.00020;
        mod[4]=0.99975;
        mod[5]=0.99956;
      }
      else if(ptTrmin==1. && ptTrmax==2.) {
        mod[0]=0.99887;
        mod[1]=1.00198;
        mod[2]=0.99901;
        mod[3]=1.00093;
        mod[4]=0.99879;
        mod[5]=1.00079;
      }
      else if(ptTrmin==2. && ptTrmax==3.) {
        mod[0]=0.99866;
        mod[1]=1.00894;
        mod[2]=1.00358;
        mod[3]=1.00077;
        mod[4]=0.99630;
        mod[5]=0.99354;
      }
      else if(ptTrmin==3. && ptTrmax==99.) {
        mod[0]=0.99586;
        mod[1]=1.03353;
        mod[2]=1.00304;
        mod[3]=0.99325;
        mod[4]=0.99551;
        mod[5]=0.99621;
      }  
    break;

    default:
      printf("Error! Incorrect definition of pT(D) range!\n");
      mod[0] = -999;
      return;
      
  }
  
} //end system 4 - 30-100%

}
