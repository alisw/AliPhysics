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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for T0 calibration                       TM-AC-AM_6-02-2006  
// equalize time shift for each time CFD channel
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliT0CalibTimeEq.h"
#include "AliLog.h"
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TProfile.h>
#include <iostream>

ClassImp(AliT0CalibTimeEq)

//________________________________________________________________
  AliT0CalibTimeEq::AliT0CalibTimeEq():TNamed(),
				       fMeanVertex(0),        
				       fRmsVertex(0)      
{
  //

}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const char* name):TNamed(),
				       fMeanVertex(0),        
				       fRmsVertex(0)      
{
  //constructor

  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const AliT0CalibTimeEq& calibda):TNamed(calibda),		
				       fMeanVertex(0),        
				       fRmsVertex(0)      
{
// copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());


}

//________________________________________________________________
AliT0CalibTimeEq &AliT0CalibTimeEq::operator =(const AliT0CalibTimeEq& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
 
  return *this;
}

//________________________________________________________________
AliT0CalibTimeEq::~AliT0CalibTimeEq()
{
  //
  // destrictor
}
//________________________________________________________________
void AliT0CalibTimeEq::Reset()
{
  //reset values

  memset(fCFDvalue,0,120*sizeof(Float_t));
  memset(fTimeEq,1,24*sizeof(Float_t));
}


//________________________________________________________________
void  AliT0CalibTimeEq::Print(Option_t*) const
{
  // print time values

  printf("\n	----	PM Arrays	----\n\n");
  printf(" Time delay CFD \n");
  for (Int_t i=0; i<24; i++) printf(" CFD  %f ",fTimeEq[i]);
  printf("\n Mean Vertex %f \n", fMeanVertex);
} 


//________________________________________________________________
Bool_t AliT0CalibTimeEq::ComputeOnlineParams(const char* filePhys)
{
  // compute online equalized time
  Float_t meandiff, sigmadiff, meanver, meancfdtime, sigmacfdtime;
  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
    // Double_t rms=0, rmscfd=0; 
  Double_t rmsver=0;
 Int_t nent=0;
  Bool_t ok=false;
  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
  }
  else
    {
      //    gFile->ls();
      ok=true;
      for (Int_t i=0; i<24; i++)
	{
	  TH1F *cfd = (TH1F*) gFile->Get(Form("CFD1minCFD%d",i+1));
	  TH1F *cfdtime = (TH1F*) gFile->Get(Form("CFD%d",i+1));
	  if(!cfd) AliWarning(Form("no histograms collected by PHYS DA for channel %i", i));
	  //      printf(" i = %d buf1 = %s\n", i, buf1);
	  if(cfd) {
	    GetMeanAndSigma(cfd, meandiff, sigmadiff);
 	    nent=cfd->GetEntries();
	    if(nent<500 || cfd->GetRMS()>20. ) {
	      ok=false;
	      AliWarning(Form("Data is not good enouph in PMT %i - mean %f rsm %f nentries %i", i,meandiff,sigmadiff , nent));
	      
	    }
	    else 
	      {
		ok=false;
		AliWarning(Form("Data is not good enouph in PMT %i , no clean peak", i));
	      }
	    if(!cfd) AliWarning(Form("no histograms collected by PHYS DA for channel %i", i));
	  }
	    //      printf(" i = %d buf1 = %s\n", i, buf1);
	  if(cfdtime) {
	    GetMeanAndSigma(cfdtime,meancfdtime, sigmacfdtime);
 	    nent=cfdtime->GetEntries();
	    if(nent<500 || sigmacfdtime>30. ) {
	      ok=false;
	      AliWarning(Form("Data is not good enouph in PMT %i CFD data - meancfdtime %f rsm %f nentries %i", i,meancfdtime, sigmacfdtime, nent));
		
	    }
	  }
	  else 
	    {
		ok=false;
		AliWarning(Form("Data is not good enouph in PMT %i , no clean peak", i));
	      }
	  
	  
	
	  SetTimeEq(i,meandiff);
	  SetTimeEqRms(i,sigmadiff);
	  SetCFDvalue(i,0,meancfdtime);
	  SetCFDvalue(i,0,sigmacfdtime);
	  if (cfd) delete cfd;
	  if (cfdtime) delete cfdtime;

	}
      TH1F *ver = (TH1F*) gFile->Get("hVertex");
      if(!ver) AliWarning("no T0 histogram collected by PHYS DA ");
      if(ver) {
	meanver = ver->GetMean();
	rmsver = ver->GetRMS();
      }
      SetMeanVertex(meanver);
      SetRmsVertex(rmsver);
      
      gFile->Close();
      delete gFile;

    }
    return ok; 
}

//________________________________________________________________
Bool_t AliT0CalibTimeEq::ComputeOfflineParams(const char* filePhys)
{
  // compute online equalized time
  Float_t meandiff, sigmadiff, meanver, meancfdtime, sigmacfdtime;
  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
  // Double_t rms=0, rmscfd=0; 
  Double_t rmsver=0;
  Int_t nent=0;
  Bool_t ok=false;
  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
  }
  else
    {
      //    gFile->ls();
      ok=true;
      TObjArray * TzeroObj = (TObjArray*) gFile->Get("fTzeroObject");
      for (Int_t i=0; i<24; i++)
	{
	  TH1F *cfddiff = (TH1F*)TzeroObj->At(i);
	  TH1F *cfdtime = (TH1F*)TzeroObj->At(i+24);
	  if(!cfddiff) AliWarning(Form("no histograms collected by PHYS DA for channel %i", i));
	  //      printf(" i = %d buf1 = %s\n", i, buf1);
	  if(cfddiff) {
 	    GetMeanAndSigma(cfddiff,meandiff, sigmadiff);
 	    nent=cfddiff->GetEntries();
	    if(nent<500 || cfddiff->GetRMS()>20. ) {
	      ok=false;
	      AliWarning(Form("Data is not good enouph in PMT %i - mean %f rsm %f nentries %i", i,meandiff,sigmadiff, nent));
	      	      
	    }
	    else 
	      {
		ok=false;
		AliWarning(Form("Data is not good enouph in PMT %i , no clean peak", i));
	      }
	  if(!cfdtime) AliWarning(Form("no histograms collected by PHYS DA for channel %i", i));
	  //      printf(" i = %d buf1 = %s\n", i, buf1);
	  if(cfdtime) {
	    GetMeanAndSigma(cfdtime,meancfdtime, sigmacfdtime);
 	    nent=cfdtime->GetEntries();
	    if(nent<500 || cfdtime->GetRMS()>30. ) {
	      ok=false;
	      AliWarning(Form("Data is not good enouph in PMT %i CFD data - mean %f rsm %f nentries %i", i,meancfdtime, sigmacfdtime, nent));
		
	    }
	  }
	  else 
	    {
		ok=false;
		AliWarning(Form("Data is not good enouph in PMT %i , no clean peak", i));
	      }
	  }
	  
	
	  SetTimeEq(i,meandiff);
	  SetTimeEqRms(i,sigmadiff);
	  SetCFDvalue(i,0,meancfdtime);
	  SetCFDvalue(i,0,sigmacfdtime);
	  if (cfddiff) delete cfddiff;
	  if (cfdtime) delete cfdtime;

	}
      TH1F *ver = (TH1F*) gFile->Get("hVertex");
      if(!ver) AliWarning("no T0 histogram collected by PHYS DA ");
      if(ver) {
	meanver = ver->GetMean();
	rmsver = ver->GetRMS();
      }
      SetMeanVertex(meanver);
      SetRmsVertex(rmsver);
      
      gFile->Close();
      delete gFile;

    }
    return ok; 
}

//________________________________________________________________________
void AliT0CalibTimeEq::GetMeanAndSigma(TH1F* hist,  Float_t &mean, Float_t &sigma) {

  const double window = 5.;  //fit window 
  double norm  = hist->Integral();  // normalize to one count
  hist->Scale(1./norm); 
 
  double meanEstimate, sigmaEstimate; 
  int maxBin;
  maxBin        =  hist->GetMaximumBin(); //position of maximum
  meanEstimate  =  hist->GetBinCenter( maxBin); // mean of gaussian sitting in maximum
  sigmaEstimate = hist->GetRMS();
  TF1* fit= new TF1("fit","gaus", meanEstimate - window*sigmaEstimate, meanEstimate + window*sigmaEstimate);
  fit->SetParameters(hist->GetBinContent(maxBin), meanEstimate, sigmaEstimate);
  hist->Fit("fit","R");

  mean  = (Float_t) fit->GetParameter(1);
  sigma = (Float_t) fit->GetParameter(2);

  delete fit;
}


