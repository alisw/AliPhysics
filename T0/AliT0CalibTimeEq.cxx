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
  Double_t mean=0, meanver=0;
  Double_t rms=0, rmsver=0;
  Int_t nent=0;
  Bool_t ok=false;
  Int_t npeaks = 20;
  Int_t sigma=3;
  Bool_t down=false;
  Int_t index[20];
  Int_t nfound=0;
  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
  }
  else
    {
      //    gFile->ls();
      ok=true;
      Char_t buf1[30];
      for (Int_t i=0; i<24; i++)
	{
	  sprintf(buf1,"CFD1minCFD%d",i+1);
	  TH1F *cfd = (TH1F*) gFile->Get(buf1);
	  if(!cfd) AliWarning(Form("no histograms collected by PHYS DA for channel %i", i));
	  //      printf(" i = %d buf1 = %s\n", i, buf1);
	  if(cfd) {
	    TSpectrum *s = new TSpectrum(2*npeaks,1);
	    nfound = s->Search(cfd,sigma," ",0.1);
	    if(nfound!=0){
	      Float_t *xpeak = s->GetPositionX();
	      TMath::Sort(nfound, xpeak, index,down);
	      Float_t xp = xpeak[index[0]];
	      Double_t hmax = xp+3*sigma;
	      Double_t hmin = xp-3*sigma;
	      cfd->GetXaxis()->SetRangeUser(hmin-1,hmax+1);
	      mean=cfd->GetMean();
	      rms=cfd->GetRMS();
	      nent=cfd->GetEntries();
	      if(nent<500 || rms>10. ) {
		ok=false;
		AliWarning(Form("Data is not good enouph in PMT %i - mean %f rsm %f nentries %i", i,mean, rms, nent));
		
	      }
	    }
	    else 
	      {
		ok=false;
		AliWarning(Form("Data is not good enouph in PMT %i , no clean peak", i));
	      }
	  }
	  
	
	   SetTimeEq(i,mean);
	   SetTimeEqRms(i,rms);
	  if (cfd) delete cfd;
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


