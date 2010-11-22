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
#include <TProfile.h>
#include <iostream>

ClassImp(AliT0CalibSeasonTimeShift)

//________________________________________________________________
  AliT0CalibSeasonTimeShift::AliT0CalibSeasonTimeShift():TNamed()
{
  //

}

//________________________________________________________________
AliT0CalibSeasonTimeShift::AliT0CalibSeasonTimeShift(const char* name):TNamed()
{
    //constructor
    
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
}

//________________________________________________________________
AliT0CalibSeasonTimeShift::AliT0CalibSeasonTimeShift(const AliT0CalibSeasonTimeShift& calibda):TNamed(calibda)		

{
// copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());


}

//________________________________________________________________
AliT0CalibSeasonTimeShift &AliT0CalibSeasonTimeShift::operator =(const AliT0CalibSeasonTimeShift& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
 
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
void  AliT0CalibSeasonTimeShift::SetT0Par(Float_t par[4],Float_t spar[4])
{
  for (Int_t i=0; i<4; i++)
    {
      fMeanPar[i] = par[i];
      fSigmaPar[i] = spar[i];
    }

}

//________________________________________________________________
void AliT0CalibSeasonTimeShift::SetT0Par(const char* filePhys)
{
  // compute online equalized time
  Float_t mean[4], sigma[4];

  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
  }
  else
    {
      gFile->ls();
 
      TString histname[4]={"meanAC", "meanA", "meanC", "resolution"};
       for (Int_t i=0; i<4; i++)
	{
	  TH1F *cfd = (TH1F*) gFile->Get(histname[i].Data());
	  if(!cfd) AliWarning(Form("no histograms collected for %s", histname[i].Data()));
	  if(cfd) {
	    TF1 *g = new TF1("g", "gaus",-2,2);
	    cfd->Fit("g"," ","Q",-2,2);
	    Double_t par[3];
	    g->GetParameters(&par[0]);
	    fMeanPar[i] = par[1];
	    fSigmaPar[i]=par[2];

	  }
	} 

	  gFile->Close();
	  delete gFile;
	  
    }
  
}
