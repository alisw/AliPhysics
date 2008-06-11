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

/*
$Log$
*/ 

#include "Riostream.h"
#include "AliFlowLYZHist2.h"
#include "AliFlowCommonConstants.h" 
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h" 
#include "TComplex.h"

class TH1D;


// Class to organize the histograms in the second run
// in the Lee Yang Zeros Flow analysis.
// Also contains methods to get values from the histograms
// which are called in AliFlowLeeYandZerosMaker::Finish().
// author: N. van der Kolk (kolk@nikhef.nl)

 

ClassImp(AliFlowLYZHist2)

  

//-----------------------------------------------------------------------

  AliFlowLYZHist2::AliFlowLYZHist2(Int_t theta):
    fHistProReNumer(0),
    fHistProImNumer(0),
    fHistProReNumerPt(0),
    fHistProImNumerPt(0),
    fHistProReNumer2D(0),
    fHistProImNumer2D(0)
{

  //constructor creating histograms 
  TString title, name;
  Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Int_t iNbinsEta = AliFlowCommonConstants::GetNbinsEta();

  Double_t  dPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  dPtMax = AliFlowCommonConstants::GetPtMax();
  Double_t  dEtaMin = AliFlowCommonConstants::GetEtaMin();	     
  Double_t  dEtaMax = AliFlowCommonConstants::GetEtaMax();
    
  //fHistProReNumer
  name = "Second_FlowPro_ReNumer";
  name +=theta;
  name +="_LYZ";
  title = "Second_FlowPro_ReNumer";
  title +=theta;
  title +="_LYZ";
  fHistProReNumer = new TProfile(name.Data(),title.Data(),iNbinsEta,dEtaMin,dEtaMax); 
  fHistProReNumer->SetXTitle("eta");
  fHistProReNumer->SetYTitle("v (%)");

  //fHistProImNumer
  name = "Second_FlowPro_ImNumer";
  name +=theta;
  name +="_LYZ";
  title = "Second_FlowPro_ImNumer";
  title +=theta;
  title +="_LYZ";
  fHistProImNumer = new TProfile(name.Data(),title.Data(),iNbinsEta,dEtaMin,dEtaMax);  
  fHistProImNumer->SetXTitle("eta");
  fHistProImNumer->SetYTitle("v (%)");

  //fHistProReNumerPt
  name = "Second_FlowPro_ReNumerPt";
  name +=theta;
  name +="_LYZ";
  title = "Second_FlowPro_ReNumerPt";
  title +=theta;
  title +="_LYZ";
  fHistProReNumerPt = new TProfile(name.Data(),title.Data(),iNbinsPt,dPtMin,dPtMax); 
  fHistProReNumerPt->SetXTitle("Pt");
  fHistProReNumerPt->SetYTitle("v (%)");

  //fHistProImNumerPt
  name = "Second_FlowPro_ImNumerPt";
  name +=theta;
  name +="_LYZ";
  title = "Second_FlowPro_ImNumerPt";
  title +=theta;
  title +="_LYZ";
  fHistProImNumerPt = new TProfile(name.Data(),title.Data(),iNbinsPt,dPtMin,dPtMax);  
  fHistProImNumerPt->SetXTitle("Pt");
  fHistProImNumerPt->SetYTitle("v (%)");

  //fHistProReNumer2D
  name = "Second_FlowPro_ReNumer2D";
  name +=theta;
  name +="_LYZ";
  title = "Second_FlowPro_ReNumer2D";
  title +=theta;
  title +="_LYZ";
  fHistProReNumer2D = new TProfile2D(name.Data(),title.Data(),iNbinsEta,dEtaMin,dEtaMax,iNbinsPt,dPtMin,dPtMax);  
  fHistProReNumer2D->SetXTitle("eta");
  fHistProReNumer2D->SetYTitle("Pt (GeV/c)");

  //fHistProImNumer2D 
  name = "Second_FlowPro_ImNumer2D";
  name +=theta;
  name +="_LYZ";
  title = "Second_FlowPro_ImNumer2D";
  title +=theta;
  title +="_LYZ";
  fHistProImNumer2D = new TProfile2D(name.Data(),title.Data(),iNbinsEta,dEtaMin,dEtaMax,iNbinsPt,dPtMin,dPtMax);  
  fHistProImNumer2D->SetXTitle("eta");
  fHistProImNumer2D->SetYTitle("Pt (GeV/c)");
}

//----------------------------------------------------------------------- 

AliFlowLYZHist2::~AliFlowLYZHist2()
{
  //deletes histograms
  delete fHistProReNumer;
  delete fHistProImNumer;
  delete fHistProReNumerPt;
  delete fHistProImNumerPt;
  delete fHistProReNumer2D;
  delete fHistProImNumer2D;
}

//----------------------------------------------------------------------- 
void AliFlowLYZHist2::Fill(Double_t f1, Double_t f2, TComplex c)
{
  //fill the real and imaginary part of fNumer

  fHistProReNumer->Fill(f1, c.Re());  
  fHistProImNumer->Fill(f1, c.Im());
   
  fHistProReNumerPt->Fill(f2, c.Re());  
  fHistProImNumerPt->Fill(f2, c.Im());
  
  fHistProReNumer2D->Fill(f1, f2, c.Re());          
  fHistProImNumer2D->Fill(f1, f2, c.Im());           
}

//-----------------------------------------------------------------------
TComplex AliFlowLYZHist2::GetNumerEta(Int_t i)
{
  //get the real and imaginary part of fNumer
  Double_t dReNumer = fHistProReNumer->GetBinContent(i);
  Double_t dImNumer = fHistProImNumer->GetBinContent(i);
  TComplex cNumer(dReNumer,dImNumer);
  //if (dNumer.Rho()==0) {cerr<<"modulus of dNumer is zero in AliFlowLYZHist2::GetNumer(Int_t i)"<<endl;}
  return cNumer;
}

//----------------------------------------------------------------------- 
TComplex AliFlowLYZHist2::GetNumerPt(Int_t i)
{
  //get the real and imaginary part of fNumer
  Double_t dReNumer = fHistProReNumerPt->GetBinContent(i);
  Double_t dImNumer = fHistProImNumerPt->GetBinContent(i);
  TComplex cNumer(dReNumer,dImNumer);
  return cNumer;
}

