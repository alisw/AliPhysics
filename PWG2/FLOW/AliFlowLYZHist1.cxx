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
#include "AliFlowLYZHist1.h"
#include "AliFlowLYZConstants.h"  //??
#include "TProfile.h"
#include "TString.h" 
#include "TComplex.h" 

class TH1D; 


// Class to organize the histograms in the first run
// in the Lee Yang Zeros Flow analysis.
// Also contains methods to find the first minimum R0
// of the generating function.
// author: N. van der Kolk (kolk@nikhef.nl)


ClassImp(AliFlowLYZHist1)

  

//-----------------------------------------------------------------------

  AliFlowLYZHist1::AliFlowLYZHist1(Int_t theta, Int_t har)
{
  //constructor creating histograms 
  Int_t fNbins = AliFlowLYZConstants::kNbins;
  Float_t fMin = AliFlowLYZConstants::fgMin;
  Float_t fMax = AliFlowLYZConstants::fgMax;
  TString title, name;
 
  
  //fHistGtheta
  name = "FirstHist_FlowLYZ_Gtheta";
  name +=theta;
  name +="_Har";
  name +=har;
  title = "FirstHist_FlowLYZ_Gtheta";
  title +=theta;
  title +="_Har";
  title +=har;
  fHistGtheta = new TH1D(name.Data(),title.Data(),fNbins,fMin,fMax);  
  fHistGtheta->SetXTitle("r");
  fHistGtheta->SetYTitle("|G^{#theta}(ir)|^{2}");
  
  //fHistProReGtheta
  name = "FirstHist_FlowProLYZ_ReGtheta";
  name +=theta;
  name +="_Har";
  name +=har;
  title = "FirstHist_FlowProLYZ_ReGtheta";
  title +=theta;
  title +="_Har";
  title +=har;
  fHistProReGtheta = new TProfile(name.Data(),title.Data(),fNbins,fMin,fMax);
  fHistProReGtheta->SetXTitle("r");
  fHistProReGtheta->SetYTitle("Re G^{#theta}(ir)");
  
  //fHistProImGtheta
  name = "FirstHist_FlowProLYZ_ImGtheta";
  name +=theta;
  name +="_Har";
  name +=har;
  title = "FirstHist_FlowProLYZ_ImGtheta";
  title +=theta;
  title +="_Har";
  title +=har;
  fHistProImGtheta = new TProfile(name.Data(),title.Data(),fNbins,fMin,fMax);
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

void AliFlowLYZHist1::Fill(Float_t f, TComplex C)
{
  //fill the histograms

  fHistProReGtheta->Fill(f, C.Re());
  fHistProImGtheta->Fill(f, C.Im());
  
}


//----------------------------------------------------------------------- 

TH1D* AliFlowLYZHist1::FillGtheta()
{
  //fill the fHistGtheta histograms
  Int_t fNbins = fHistGtheta->GetNbinsX();
  for (Int_t bin=1;bin<=fNbins;bin++)
	{
	  //get bincentre of bins in histogram
	  Float_t fBin = fHistGtheta->GetBinCenter(bin); 
	  Float_t fRe = fHistProReGtheta->GetBinContent(bin);
	  Float_t fIm = fHistProImGtheta->GetBinContent(bin);
	  TComplex fGtheta(fRe,fIm);
	  //fill fHistGtheta with the modulus squared of fGtheta
	  fHistGtheta->Fill(fBin,fGtheta.Rho2());
	}

  return fHistGtheta;
}

//----------------------------------------------------------------------- 

Float_t AliFlowLYZHist1::GetR0()
{
  //find the first minimum of the square of the modulus of Gtheta 

  Int_t fNbins = fHistGtheta->GetNbinsX();
  Float_t fR0 = 0; 

  for (Int_t b=2;b<fNbins;b++)
    {
      Float_t fG0 = fHistGtheta->GetBinContent(b);
      Float_t fGnext = fHistGtheta->GetBinContent(b+1);
      Float_t fGnextnext = fHistGtheta->GetBinContent(b+2);
      //fprintf(fTemp,"%f\t %f\t %f\n",fG0,fGnext,fGnextnext);

      if (fGnext > fG0 && fGnextnext > fG0)
	{
	  Float_t fGlast = fHistGtheta->GetBinContent(b-1);
	  Float_t fXlast = fHistGtheta->GetBinCenter(b-1);
	  Float_t fX0 = fHistGtheta->GetBinCenter(b);
	  Float_t fXnext = fHistGtheta->GetBinCenter(b+1);

	  fR0 = fX0 - ((fX0-fXlast)*(fX0-fXlast)*(fG0-fGnext) - (fX0-fXnext)*(fX0-fXnext)*(fG0-fGlast))/
	    (2.*((fX0-fXlast)*(fG0-fGnext) - (fX0-fXnext)*(fG0-fGlast))); //interpolated minimum
	  
	  break; //stop loop if minimum is found
	} //if

    }//b

      
  return fR0;
}
   
//----------------------------------------------------------------------- 
Float_t AliFlowLYZHist1::GetBinCenter(Int_t i)
{
  //gets bincenter of histogram
  Float_t fR = fHistGtheta->GetBinCenter(i);
  return fR;

}

//----------------------------------------------------------------------- 

Int_t AliFlowLYZHist1::GetNBins()
{
  //gets fNbins
  Int_t fNbins = fHistGtheta->GetNbinsX();
  return fNbins;

}

