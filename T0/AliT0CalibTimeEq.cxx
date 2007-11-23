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
// class for T0 calibration                       TM-AC-AM_6-02-2006         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliT0CalibTimeEq.h"
#include "AliT0LookUpValue.h"
#include "AliLog.h"
#include "AliRun.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TSystem.h>
#include <Riostream.h>
#include <TSpectrum.h>
#include <TVirtualFitter.h>
#include <TProfile.h>

#include <Riostream.h>
#include <string>

ClassImp(AliT0CalibTimeEq)

//________________________________________________________________
  AliT0CalibTimeEq::AliT0CalibTimeEq():   TNamed(),
                                      fTimeDelayTVD(0),
                                      fMeanT0(0) 
{
  //
}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const char* name):TNamed(),
                                      fTimeDelayTVD(0),
                                      fMeanT0(0)
{
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const AliT0CalibTimeEq& calibda) :
  TNamed(calibda),		
  fTimeDelayTVD(0),
  fMeanT0(0)
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
}
//________________________________________________________________
void AliT0CalibTimeEq::Reset()
{
  memset(fCFDvalue,0,120*sizeof(Float_t));
  memset(fTimeEq,1,24*sizeof(Float_t));
}


//________________________________________________________________
void  AliT0CalibTimeEq::Print(Option_t*) const
{

  printf("\n	----	PM Arrays	----\n\n");
  printf(" Time delay CFD \n");
  for (Int_t i=0; i<24; i++) printf(" CFD  %f ",fTimeEq[i]);
} 


//________________________________________________________________
void AliT0CalibTimeEq::ComputeOnlineParams(char* name1, char* name2, char* canv, Int_t npeaks, Double_t sigma)
{
  TFile *gFile = TFile::Open("daq1Phys.root");
  gSystem->Load("libSpectrum");
  npeaks = 20;
  sigma=3.;
  Bool_t down=false;
  Int_t index[20];
  Char_t buf1[15];
  Char_t temp[10];
  Float_t p[12][3]={0.,0.,0.};
  for (Int_t i=12; i<24; i++)
    {
      sprintf(buf1,name1);
      sprintf(temp,"%i",i+1);
      strcat (buf1,temp);
      strcat (buf1,name2);
      TH1F *cfd = (TH1F*) gFile->Get(buf1);
      TSpectrum *s = new TSpectrum(2*npeaks,1.);
      Int_t nfound = s->Search(cfd,sigma," ",0.2);
      if(nfound!=0){
        Float_t *xpeak = s->GetPositionX();
        TMath::Sort(nfound, xpeak, index,down);
        Float_t xp = xpeak[index[0]];
        Float_t hmax = xp+3*sigma;
        Float_t hmin = xp-3*sigma;
        cfd->GetXaxis()->SetRange((Int_t)hmin-20,(Int_t)hmax+20);
        TF1 *g1 = new TF1("g1", "gaus", hmin, hmax);
              cfd->Fit("g1","IR");
	      
	      for(Int_t j =0; j<3; j++){
                p[i][j] = g1->GetParameter(j);
                SetCFDvalue(i, j, p[i][j]);
	      }
	      SetCFDvalue(i, 3, hmin);
	      SetCFDvalue(i, 4, hmax);
	      SetTimeEq(i,p[i][2]);
      }
    }
 gFile->Close();
  delete gFile;
}

 
