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

/* $Id: AliT0CalibSeasonTimeShift.cxx 42881 2010-08-16 10:59:14Z alla $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for T0 calibration                       TM-AC-AM_6-02-2006  
// equalize time shift for each time CFD channel
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliT0CalibSeasonTimeShift.h"
#include "AliLog.h"
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TAxis.h>
#include <iostream>

ClassImp(AliT0CalibSeasonTimeShift)

//________________________________________________________________
AliT0CalibSeasonTimeShift::AliT0CalibSeasonTimeShift():TNamed(),
  fT0vsMult(0)
{
  //
  for (Int_t i=0; i<4; i++) 
    fMeanPar[i] = fSigmaPar[i] = 0;
    
  fT0vsMult.SetOwner(kTRUE);
}

//________________________________________________________________
AliT0CalibSeasonTimeShift::AliT0CalibSeasonTimeShift(const char* name):TNamed(),
  fT0vsMult(0)
{
    //constructor
    
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data()); 
  
  for (Int_t i=0; i<4; i++)
    fMeanPar[i] = fSigmaPar[i] = 0; 
 
}

//________________________________________________________________
AliT0CalibSeasonTimeShift::AliT0CalibSeasonTimeShift(const AliT0CalibSeasonTimeShift& calibda):TNamed(calibda)		

{
// copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  ((AliT0CalibSeasonTimeShift &) calibda).Copy(*this);


}

//________________________________________________________________
AliT0CalibSeasonTimeShift &AliT0CalibSeasonTimeShift::operator =(const AliT0CalibSeasonTimeShift& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  if (this != &calibda) ((AliT0CalibSeasonTimeShift &) calibda).Copy(*this);
 
  return *this;
}

//________________________________________________________________
AliT0CalibSeasonTimeShift::~AliT0CalibSeasonTimeShift()
{
  //
  // destrictor
}


//________________________________________________________________
void  AliT0CalibSeasonTimeShift::Print(Option_t*) const
{
  // print time values

  printf("\n	----	T0 results	----\n\n");
  printf(" (T0A+T0C)/2 = %f; T0A = %f; T0C = %f; resolution = %f  \n", fMeanPar[0], fMeanPar[1],fMeanPar[2],fMeanPar[3]);
  printf(" sigma(T0A+T0C)/2 = %f; sigma(T0 = %f; sigma(T0C) = %f; sigma(resolution) = %f  \n" , fSigmaPar[0], fSigmaPar[1], fSigmaPar[2],fSigmaPar[3]);
  
} 

//________________________________________________________________
Bool_t  AliT0CalibSeasonTimeShift::SetT0Par(Float_t par[4],Float_t spar[4])
{
  Bool_t ok=false;
 for (Int_t i=0; i<4; i++)
    {
      fMeanPar[i] = par[i];
      fSigmaPar[i] = spar[i];
      if ( fSigmaPar[i] == 0 ||  fSigmaPar[i] > 500) ok = false;
    }
 return ok;
}

//________________________________________________________________
Int_t AliT0CalibSeasonTimeShift::SetT0Par(const char* filePhys, Float_t *cdbtime)
{
  // compute shifts fo T0A, T0C, T0AC and resolution
  //Status : writeok:
  // 0 OK
  // 2000 - no data file
  // 300 no one histogram or it is empty
  //-100 peak is very narrow

  Float_t mean=0, sigma=0;
  Int_t ok = 0;
  TH1F *hcfd = NULL;
  TH2F *hT0mult = NULL;
  TObjArray * tzeroObj = NULL;

  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
    return 2000;
  }
  else {
    tzeroObj = dynamic_cast<TObjArray*>(gFile->Get("T0Calib"));
    TString histname[4]={"fTzeroORAplusORC", "fTzeroORA", "fTzeroORC",  "fResolution"};
    for (Int_t i=0; i<4; i++)
      {
	if(hcfd) hcfd->Reset();
	if(tzeroObj) 
	  hcfd = (TH1F*)tzeroObj->FindObject( histname[i].Data());
	else
	  hcfd =  (TH1F*)gFile ->Get(histname[i].Data());

	if(!hcfd) {
	  AliError(Form("no histograms collected for %s", histname[i].Data()));
	  return -300;
	}
	if(hcfd) {
	  if( hcfd->GetEntries() == 0) {
	  AliError(Form("%s histogram is empty", histname[i].Data()));
	  return -300;
	  }
	  GetMeanAndSigma(hcfd, mean, sigma);
	  if (sigma == 0 || sigma > 600 || hcfd->GetEntries()<50 ){ //!!!
	    AliError(Form("%s low statsitics or bad histogram, OCDB value is = %f", histname[i].Data(), cdbtime[i]) );
	    return 400;
	  }
	  if ( sigma > 0 && sigma < 600 && hcfd->GetEntries()>=50) //!!!
	    { 
	      fMeanPar[i] =   mean;
	      fSigmaPar[i] = sigma;
	    }
	}
      } 
    /*
    TString histname2D[4]={"hT0AC", "hT0A", "hT0C",  "hResolution"};
    for (Int_t i=0; i<4; i++)
      {
	if(hT0mult) hT0mult->Reset();
	if(tzeroObj) 
	 hT0mult  = (TH2F*)tzeroObj->FindObject( histname2D[i].Data());
	else
	 hT0mult =  (TH2F*)gFile ->Get(histname2D[i].Data());
	if(hT0mult) {
	  Int_t nbins = hT0mult->GetXaxis()->GetNbins();
	  Float_t meanprof[nbins], sigmares[nbins], mult[nbins];
	  Int_t npoints=0;
	  for (int ibin=1; ibin<nbins-2; ibin++) {
	    mult[ibin-1]= hT0mult-> GetXaxis()->GetBinCenter(ibin);
	    TH1D *proj = hT0mult->ProjectionY(Form("prY%i",ibin),ibin, ibin+1);
	      if(proj->GetEntries()>100) {
		GetMeanAndSigma( (TH1F*)proj, mean, sigma);
		sigmares[ibin-1]=sigma;
		meanprof[ibin-1]=mean;
		npoints++;
		printf("@@@ npoints %i \n",npoints);
	      }
	      else 
		{
		sigmares[ibin-1]=fSigmaPar[i]; 
		meanprof[ibin-1]=fMeanPar[i];
		}
	      delete proj; 
	  }
	  TGraph *gr = NULL;
	  if(i<3) gr = new TGraph(npoints-1, mult, meanprof);
	  else
	    gr = new TGraph(npoints-1, mult, sigmares);
	  fT0vsMult.AddAtAndExpand(gr,i);
	  //  gr->Delete();
	}
      }
    */  
  }
  
  
  gFile->Close();
  delete gFile;
  return ok;
}
//________________________________________________________________________
void AliT0CalibSeasonTimeShift::GetMeanAndSigma(TH1F* hist,  Float_t &mean, Float_t &sigma) {

  const double window =3.;  //fit window 
 
  double meanEstimate, sigmaEstimate; 
  int maxBin;
  maxBin        =  hist->GetMaximumBin(); //position of maximum
  meanEstimate  =  hist->GetBinCenter( maxBin); // mean of gaussian sitting in maximum
  sigmaEstimate = hist->GetRMS();
  TF1* fit= new TF1("fit","gaus", meanEstimate - window*sigmaEstimate, meanEstimate + window*sigmaEstimate);
  fit->SetParameters(hist->GetBinContent(maxBin), meanEstimate, sigmaEstimate);
  hist->Fit("fit","RQ","S");

  mean  = (Float_t) fit->GetParameter(1);
  sigma = (Float_t) fit->GetParameter(2);

  delete fit;
}
