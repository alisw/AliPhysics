
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
  Double_t rmsver=0;
  Int_t nent=0;
  Bool_t ok=false;
  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
  }
  else
    {
      ok=true;
      for (Int_t i=0; i<24; i++)
	{

	  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
	  TH1F *cfd = (TH1F*) gFile->Get(Form("CFD1minCFD%d",i+1));
	  TH1F *cfdtime = (TH1F*) gFile->Get(Form("CFD%d",i+1));
	  if(!cfd) {
	    AliWarning(Form("no Diff histograms collected by PHYS DA for channel %i", i));
	  }
	  if(!cfdtime) {
	    AliWarning(Form("no CFD histograms collected by PHYS DA for channel %i", i));
	  }
	  if(cfd) {
	    nent = Int_t(cfd->GetEntries());
	    if(nent>50)  { 
	      if(cfd->GetRMS()>1.5 )
		GetMeanAndSigma(cfd, meandiff, sigmadiff);
	      if(cfd->GetRMS()<=1.5) 
		{
		  meandiff = cfd->GetMean();
		  sigmadiff=cfd->GetRMS();
		}
	      Int_t   maxBin = cfd->GetMaximumBin(); 
	      Double_t  meanEstimate = cfd->GetBinCenter( maxBin); 
	      if(TMath::Abs(meanEstimate - meandiff) > 20 ) meandiff = meanEstimate; 
	    }
	    else 
	      {
		//	ok=false;
		AliWarning(Form(" Not  enouph data in PMT %i- PMT1:  %i ", i, nent));
	      }
	  }
	  if(cfdtime) {
	    nent = Int_t(cfdtime->GetEntries());
	    if(nent > 50 )  { //!!!!!!!!!!
	      if(cfdtime->GetRMS()>1.5 )
		GetMeanAndSigma(cfdtime,meancfdtime, sigmacfdtime);
	      if(cfdtime->GetRMS()<=1.5) 
		{
		  meancfdtime = cfdtime->GetMean();
		  sigmacfdtime = cfdtime->GetRMS();
		}
	      Int_t   maxBin = cfdtime->GetMaximumBin(); 
	      Double_t  meanEstimate = cfdtime->GetBinCenter( maxBin); 
	      if(TMath::Abs(meanEstimate - meancfdtime) > 20 ) meancfdtime = meanEstimate; 
	    }
	    else 
	      {
		//	ok=false;
		AliWarning(Form(" Not  enouph data in PMT in CFD peak %i - %i ", i, nent));
	      }
	  }
	  SetTimeEq(i,meandiff);
	  SetTimeEqRms(i,sigmadiff);
	  SetCFDvalue(i,0,meancfdtime);
	  if (cfd) delete cfd;
	  if (cfdtime) delete cfdtime;
	  
	}
      TH1F *ver = (TH1F*) gFile->Get("hVertex");
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
  Bool_t AliT0CalibTimeEq::ComputeOfflineParams(const char* filePhys, Float_t *timecdb, Float_t *cfdvalue, Int_t badpmt)
{
  // compute offline equalized time
  Float_t meandiff, sigmadiff, meanver, meancfdtime, sigmacfdtime;
  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
  Int_t nent=0;
  Bool_t ok=false;
  TH1F *cfddiff = NULL; 
  TH1F *cfdtime = NULL;
  TObjArray * TzeroObj = NULL;

  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
    return ok;
  }
  else
    {
      meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
      ok=true;
      TDirectory *dr = (TDirectory*) gFile->Get("T0Calib");
      if (dr)   TzeroObj = (TObjArray*) dr->Get("T0Calib");
      
      for (Int_t i=0; i<24; i++)
	{
	  if (i != badpmt) {
	    if(TzeroObj) {
	      cfddiff = (TH1F*)TzeroObj->FindObject(Form("CFD1minCFD%d",i+1));
	      cfdtime = (TH1F*)TzeroObj->FindObject(Form("CFD%d",i+1));
	    }
	    else
	      {
		cfddiff = (TH1F*)gFile->Get(Form("CFD1minCFD%d",i+1));
		cfdtime = (TH1F*)gFile->Get(Form("CFD%d",i+1));
		
	      }
	    if(!cfddiff ) {
	      AliWarning(Form("no  histograms collected by pass0 for diff channel %i", i));      
	      meandiff = timecdb[i];
	      sigmadiff = 0; 
	    }
	    if(cfddiff) {
	      nent = Int_t(cfddiff->GetEntries());
	      if(nent>100 )  { //!!!!!
		if(cfddiff->GetRMS()>1.5 )
		  GetMeanAndSigma(cfddiff, meandiff, sigmadiff);
		if(cfddiff->GetRMS()<=1.5) 
		  {
		    meandiff = cfddiff->GetMean();
		    sigmadiff = cfddiff->GetRMS();
		  }
		Int_t   maxBin = cfddiff->GetMaximumBin(); 
		Double_t  meanEstimate = cfddiff->GetBinCenter( maxBin); 
		if(TMath::Abs(meanEstimate - meandiff) > 20 ) meandiff = meanEstimate; 	      
	      }
	      else 
		{
		  AliWarning(Form(" Not  enouph data in PMT %i- PMT1:  %i ", i, nent));
		//  ok=false;
		  meandiff = timecdb[i];
		  sigmadiff = 0; 
		  
		}
	    }	    
	    
	    if(!cfdtime ) {
	      AliWarning(Form("no  histograms collected by pass0 for time channel %i", i));
	      meancfdtime = cfdvalue[i];
	      ok = false;
	    }
	    if(cfdtime) {
	      nent = Int_t(cfdtime->GetEntries());
	      if(nent>100 )  { //!!!!!
		if(cfdtime->GetRMS()>1.5 )
		  GetMeanAndSigma(cfdtime,meancfdtime, sigmacfdtime);
		if(cfdtime->GetRMS()<=1.5) 
		  {
		    meancfdtime = cfdtime->GetMean();
		    sigmacfdtime=cfdtime->GetRMS();
		    if(cfdtime->GetRMS() == 0 || cfdtime->GetMean() ==0 ) ok = false;
		}
		Int_t   maxBin = cfdtime->GetMaximumBin(); 
		Double_t  meanEstimate = cfdtime->GetBinCenter( maxBin); 
		if(TMath::Abs(meanEstimate - meancfdtime) > 20 ) meancfdtime = meanEstimate; 
	    }
	    else 
	      {
		AliWarning(Form(" Not  enouph data in PMT in CFD peak %i - %i ", i, nent));
         	ok = false;
	      }
	  }
	  
	  SetTimeEq(i,meandiff);
	  SetTimeEqRms(i,sigmadiff);
	  SetCFDvalue(i,0, meancfdtime );
//	  printf(" pmt %i diff %f sigma %f meancfdtime %f cdbtime %f \n",i, meandiff, sigmadiff, meancfdtime, fCFDvalue[i][0]);
	  if (cfddiff) cfddiff->Reset();
	  if (cfdtime) cfdtime->Reset();
	  } //bad pmt
	}
      
      gFile->Close();
      delete gFile;

    }
    return ok; 
   }

//________________________________________________________________________
void AliT0CalibTimeEq::GetMeanAndSigma(TH1F* hist,  Float_t &mean, Float_t &sigma) {
  
  const double window = 2.;  //fit window 
  
  double meanEstimate, sigmaEstimate; 
  int maxBin;
  maxBin        =  hist->GetMaximumBin(); //position of maximum
  meanEstimate  =  hist->GetBinCenter( maxBin); // mean of gaussian sitting in maximum
  // sigmaEstimate = hist->GetRMS();
  sigmaEstimate = 10;
  TF1* fit= new TF1("fit","gaus", meanEstimate - window*sigmaEstimate, meanEstimate + window*sigmaEstimate);
  fit->SetParameters(hist->GetBinContent(maxBin), meanEstimate, sigmaEstimate);
  hist->Fit("fit","RQ","Q");

  mean  = (Float_t) fit->GetParameter(1);
  sigma = (Float_t) fit->GetParameter(2);

  delete fit;
}


