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
  AliT0CalibTimeEq::AliT0CalibTimeEq():TNamed()
{
  //

}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const char* name):TNamed()
{
  //constructor

  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const AliT0CalibTimeEq& calibda):TNamed(calibda)		
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
void AliT0CalibTimeEq::ComputeOnlineParams(const char* filePhys)
{
  // compute online equalized time
  Double_t mean=0, meanver=0;
  gFile = TFile::Open(filePhys);
  
  if(!gFile) {
    AliError("No input PHYS data found ");
  }
  else
    {
      Char_t buf1[30];
      for (Int_t i=0; i<24; i++)
	{
	  sprintf(buf1,"CFD1-CFD%d",i+1);
	  TH1F *cfd = (TH1F*) gFile->Get(buf1);
	  if(!cfd) AliWarning(Form("no histograms collected by PHYS DA for channel %i", i));
	  //      printf(" i = %d buf1 = %s\n", i, buf1);
	  if(cfd) mean=cfd->GetMean();
	  SetTimeEq(i,mean);
	  if (cfd) delete cfd;
	}
      TH1F *ver = (TH1F*) gFile->Get("hVertex");
      if(!ver) AliWarning("no Vertex histograms collected by PHYS DA for Zvertex");
      if(ver)  meanver = ver->GetMean();
      SetMeanVertex(meanver);
      
      gFile->Close();
      delete gFile;
    } 
    }


