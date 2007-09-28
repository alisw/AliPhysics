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
#include "AliFlowLYZConstants.h" //??
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

  AliFlowLYZHist2::AliFlowLYZHist2(Int_t theta, Int_t har)
{
  //constructor creating histograms 
  TString title, name;
  Int_t fEtaBins = AliFlowLYZConstants::kEtaBins;
  Int_t fPtBins = AliFlowLYZConstants::kPtBins;
  Float_t fEtaMin = AliFlowLYZConstants::fgEtaMin;
  Float_t fEtaMax = AliFlowLYZConstants::fgEtaMax;
  Float_t fPtMin = AliFlowLYZConstants::fgPtMin;
  Float_t fPtMax = AliFlowLYZConstants::fgPtMax;

  
  //fHistProReNumer
  name = "SecondHist_FlowProLYZ_ReNumer";
  name +=theta;
  name +="_Har";
  name +=har;
  title = "SecondHist_FlowProLYZ_ReNumer";
  title +=theta;
  title +="_Har";
  title +=har;
  fHistProReNumer = new TProfile(name.Data(),title.Data(),fEtaBins,fEtaMin,fEtaMax); 
  fHistProReNumer->SetXTitle("eta");
  fHistProReNumer->SetYTitle("v (%)");

  //fHistProImNumer
  name = "SecondHist_FlowProLYZ_ImNumer";
  name +=theta;
  name +="_Har";
  name +=har;
  title = "SecondHist_FlowProLYZ_ImNumer";
  title +=theta;
  title +="_Har";
  title +=har;
  fHistProImNumer = new TProfile(name.Data(),title.Data(),fEtaBins,fEtaMin,fEtaMax);  
  fHistProImNumer->SetXTitle("eta");
  fHistProImNumer->SetYTitle("v (%)");

  //fHistProReNumerPt
  name = "SecondHist_FlowProLYZ_ReNumerPt";
  name +=theta;
  name +="_Har";
  name +=har;
  title = "SecondHist_FlowProLYZ_ReNumerPt";
  title +=theta;
  title +="_Har";
  title +=har;
  fHistProReNumerPt = new TProfile(name.Data(),title.Data(),fPtBins,fPtMin,fPtMax); 
  fHistProReNumerPt->SetXTitle("Pt");
  fHistProReNumerPt->SetYTitle("v (%)");

  //fHistProImNumerPt
  name = "SecondHist_FlowProLYZ_ImNumerPt";
  name +=theta;
  name +="_Har";
  name +=har;
  title = "SecondHist_FlowProLYZ_ImNumerPt";
  title +=theta;
  title +="_Har";
  title +=har;
  fHistProImNumerPt = new TProfile(name.Data(),title.Data(),fPtBins,fPtMin,fPtMax);  
  fHistProImNumerPt->SetXTitle("Pt");
  fHistProImNumerPt->SetYTitle("v (%)");

  //fHistProReNumer2D
  name = "SecondHist_FlowProLYZ_ReNumer2D";
  name +=theta;
  name +="_Har";
  name +=har;
  title = "SecondHist_FlowProLYZ_ReNumer2D";
  title +=theta;
  title +="_Har";
  title +=har;
  fHistProReNumer2D = new TProfile2D(name.Data(),title.Data(),fEtaBins,fEtaMin,fEtaMax,fPtBins,fPtMin,fPtMax);  
  fHistProReNumer2D->SetXTitle("eta");
  fHistProReNumer2D->SetYTitle("Pt (GeV/c)");

  //fHistProImNumer2D 
  name = "SecondHist_FlowProLYZ_ImNumer2D";
  name +=theta;
  name +="_Har";
  name +=har;
  title = "SecondHist_FlowProLYZ_ImNumer2D";
  title +=theta;
  title +="_Har";
  title +=har;
  fHistProImNumer2D = new TProfile2D(name.Data(),title.Data(),fEtaBins,fEtaMin,fEtaMax,fPtBins,fPtMin,fPtMax);  
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

void AliFlowLYZHist2::Fill(Float_t f1, Float_t f2, TComplex C)
{
  //fill the real and imaginary part of fNumer

  fHistProReNumer->Fill(f1, C.Re());  
  fHistProImNumer->Fill(f1, C.Im());
   
  fHistProReNumerPt->Fill(f2, C.Re());  
  fHistProImNumerPt->Fill(f2, C.Im());
  
  fHistProReNumer2D->Fill(f1, f2, C.Re());          
  fHistProImNumer2D->Fill(f1, f2, C.Im());           
}

//-----------------------------------------------------------------------
TComplex AliFlowLYZHist2::GetfNumer(Int_t i)
{
  //get the real and imaginary part of fNumer
  Float_t fReNumer = fHistProReNumer->GetBinContent(i);
  Float_t fImNumer = fHistProImNumer->GetBinContent(i);
  TComplex fNumer(fReNumer,fImNumer);
  //if (fNumer.Rho()==0) {cerr<<"modulus of fNumer is zero in AliFlowLYZHist2::GetfNumer(Int_t i)"<<endl;}
  return fNumer;
}

//----------------------------------------------------------------------- 
TComplex AliFlowLYZHist2::GetfNumerPt(Int_t i)
{
  //get the real and imaginary part of fNumer
  Float_t fReNumer = fHistProReNumerPt->GetBinContent(i);
  Float_t fImNumer = fHistProImNumerPt->GetBinContent(i);
  TComplex fNumer(fReNumer,fImNumer);
  return fNumer;
}

//----------------------------------------------------------------------- 
Int_t AliFlowLYZHist2::GetNprime(Int_t i)
{
  //Get the number of entries in the bin 
  Int_t Nprime = fHistProReNumer->GetBinEntries(i);
  return Nprime;
}

//----------------------------------------------------------------------- 
Int_t AliFlowLYZHist2::GetNprimePt(Int_t i)
{
  //Get the number of entries in the bin 
  Int_t Nprime = fHistProReNumerPt->GetBinEntries(i);
  return Nprime;
}
