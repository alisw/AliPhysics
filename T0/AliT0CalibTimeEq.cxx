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
} 


//________________________________________________________________
void AliT0CalibTimeEq::ComputeOnlineParams(const char* filePhys)
{
  // compute online equalized time
  // Int_t npeaks = 20;
  // Double_t sigma = 4.;
  //  Bool_t down=false;
  // Int_t index[20];

  gFile = TFile::Open(filePhys);
  gFile->ls();
  Char_t buf1[30];
  for (Int_t i=0; i<24; i++)
  {
    //    if(i<12) sprintf(buf1,"CFD1minCFD%d",i+1);
    // if(i>11) sprintf(buf1,"CFD13minCFD%d",i+1);
    if(i<12) sprintf(buf1,"T0_C_%i_CFD",i+1);
    if(i>11) sprintf(buf1,"T0_A_%i_CFD",i+1-12);
    printf(" i = %d buf1 = %s\n", i, buf1);
    TH1F *cfd = (TH1F*) gFile->Get(buf1);
    
    //    printf(" i = %d buf1 = %s\n", i, buf1);
    Double_t mean=cfd->GetMean();
    printf(" ipmt %i mean = %f \n", i,mean) ;

    /*    
    TSpectrum *s = new TSpectrum(2*npeaks,1.);
    // printf(" buf1 = %s cfd = %x\n", buf1, cfd);
    Int_t nfound = s->Search(cfd,sigma,"goff",0.2);
    printf(" nfound = %d\n", nfound);
    if(nfound!=0)
    {
      Float_t *xpeak = s->GetPositionX();
      TMath::Sort(nfound, xpeak, index,down);
      Float_t xp = xpeak[index[0]];
      Float_t hmax = xp+3*sigma;
      Float_t hmin = xp-3*sigma;
      cfd->GetXaxis()->SetRangeUser(hmin-10,hmax+10);
      TF1 *g1 = new TF1("g1", "gaus", hmin, hmax);
      cfd->Fit("g1","IRQN");
     printf(" ipmt %i fit mean = %f \n", i, g1->GetParameter(1)) ;
     
    }
    */
      SetTimeEq(i,Int_t(mean));
  }
  // delete cfd;
  
   gFile->Close();
   delete gFile;
   printf("\n\n");
   for(int j=0;j<24;j++)
   {
                printf("fTimeEq[%d]=%f\n",j,fTimeEq[j]);
   }

}


