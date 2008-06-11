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

//#include "Riostream.h"  //in case one wants to make print statements
#include "TProfile.h"
#include "TString.h" 
#include "TComplex.h"  
#include "AliFlowLYZConstants.h" 
#include "AliFlowLYZHist1.h"


class TH1D; 


// Class to organize the histograms in the first run
// in the Lee Yang Zeros Flow analysis.
// Also contains methods to find the first minimum R0
// of the generating function.
// author: N. van der Kolk (kolk@nikhef.nl)


ClassImp(AliFlowLYZHist1)

  

//-----------------------------------------------------------------------

  AliFlowLYZHist1::AliFlowLYZHist1(Int_t theta):
    fHistGtheta(0),
    fHistProReGtheta(0),
    fHistProImGtheta(0)
{

  //constructor creating histograms 
  Int_t iNbins = AliFlowLYZConstants::kNbins;
  Double_t dMin = AliFlowLYZConstants::fgMin;
  Double_t dMax = AliFlowLYZConstants::fgMax;
  TString title, name;
 
  
  //fHistGtheta
  name = "First_Flow_Gtheta";
  name +=theta;
  name +="_LYZ";
  title = "First_Flow_Gtheta";
  title +=theta;
  title +="_LYZ";
  fHistGtheta = new TH1D(name.Data(),title.Data(),iNbins,dMin,dMax);  
  fHistGtheta->SetXTitle("r");
  fHistGtheta->SetYTitle("|G^{#theta}(ir)|^{2}");
  
  //fHistProReGtheta
  name = "First_FlowPro_ReGtheta";
  name +=theta;
  name +="_LYZ";
  title = "First_FlowPro_ReGtheta";
  title +=theta;
  title +="_LYZ";
  fHistProReGtheta = new TProfile(name.Data(),title.Data(),iNbins,dMin,dMax);
  fHistProReGtheta->SetXTitle("r");
  fHistProReGtheta->SetYTitle("Re G^{#theta}(ir)");
  
  //fHistProImGtheta
  name = "First_FlowPro_ImGtheta";
  name +=theta;
  name +="_LYZ";
  title = "First_FlowPro_ImGtheta";
  title +=theta;
  title +="_LYZ";
  fHistProImGtheta = new TProfile(name.Data(),title.Data(),iNbins,dMin,dMax);
  fHistProImGtheta->SetXTitle("r");
  fHistProImGtheta->SetYTitle("Im G^{#theta}(ir)");
      
}
  
//----------------------------------------------------------------------- 

AliFlowLYZHist1::~AliFlowLYZHist1()
{
  //deletes histograms
  delete fHistGtheta;
  delete fHistProReGtheta;
  delete fHistProImGtheta;
}

//----------------------------------------------------------------------- 

void AliFlowLYZHist1::Fill(Double_t f, TComplex c)
{
  //fill the histograms

  fHistProReGtheta->Fill(f, c.Re());
  fHistProImGtheta->Fill(f, c.Im());
}

//----------------------------------------------------------------------- 

TH1D* AliFlowLYZHist1::FillGtheta()
{
  //fill the fHistGtheta histograms
  Int_t iNbins = fHistGtheta->GetNbinsX();
  for (Int_t bin=1;bin<=iNbins;bin++)
	{
	  //get bincentre of bins in histogram
	  Double_t dBin = fHistGtheta->GetBinCenter(bin); 
	  Double_t dRe = fHistProReGtheta->GetBinContent(bin);
	  Double_t dIm = fHistProImGtheta->GetBinContent(bin);
	  TComplex cGtheta(dRe,dIm);
	  //fill fHistGtheta with the modulus squared of cGtheta
	  fHistGtheta->Fill(dBin,cGtheta.Rho2());
	}

  return fHistGtheta;
}

//----------------------------------------------------------------------- 

Double_t AliFlowLYZHist1::GetR0()
{
  //find the first minimum of the square of the modulus of Gtheta 

  Int_t iNbins = fHistGtheta->GetNbinsX();
  Double_t dR0 = 0.; 

  for (Int_t b=2;b<iNbins;b++)
    {
      Double_t dG0 = fHistGtheta->GetBinContent(b);
      Double_t dGnext = fHistGtheta->GetBinContent(b+1);
      Double_t dGnextnext = fHistGtheta->GetBinContent(b+2);
      
      if (dGnext > dG0 && dGnextnext > dG0)
	{
	  Double_t dGlast = fHistGtheta->GetBinContent(b-1);
	  Double_t dXlast = fHistGtheta->GetBinCenter(b-1);
	  Double_t dX0 = fHistGtheta->GetBinCenter(b);
	  Double_t dXnext = fHistGtheta->GetBinCenter(b+1);

	  dR0 = dX0 - ((dX0-dXlast)*(dX0-dXlast)*(dG0-dGnext) - (dX0-dXnext)*(dX0-dXnext)*(dG0-dGlast))/
	    (2.*((dX0-dXlast)*(dG0-dGnext) - (dX0-dXnext)*(dG0-dGlast))); //interpolated minimum
	  
	  break; //stop loop if minimum is found
	} //if

    }//b

      
  return dR0;
}
   
//----------------------------------------------------------------------- 
Double_t AliFlowLYZHist1::GetBinCenter(Int_t i)
{
  //gets bincenter of histogram
  Double_t dR = fHistGtheta->GetBinCenter(i);
  return dR;
}

//----------------------------------------------------------------------- 

Int_t AliFlowLYZHist1::GetNBins()
{
  //gets iNbins
  Int_t iNbins = fHistGtheta->GetNbinsX();
  return iNbins;
}

