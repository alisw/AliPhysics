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
