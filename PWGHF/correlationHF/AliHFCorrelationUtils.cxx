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
if(system==3) { //start system 2, 0-100%

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


}
