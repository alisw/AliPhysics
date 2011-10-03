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

/* $Id: AliZDCTDCCalib.cxx 46092 2010-12-16 13:18:21Z coppedis $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ZDC calibration  -> values for TDC offset calibration           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCTDCCalib.h"

ClassImp(AliZDCTDCCalib)

//________________________________________________________________
AliZDCTDCCalib::AliZDCTDCCalib():
TNamed()
{
  Reset();
}

//________________________________________________________________
AliZDCTDCCalib::AliZDCTDCCalib(const char* name):
TNamed()
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  for(Int_t i=0; i<6; i++){
    fMeanTDC[i] = 0.;
    fWidthTDC[i] = 0.;
  }
  
  
}

//________________________________________________________________
AliZDCTDCCalib::AliZDCTDCCalib(const AliZDCTDCCalib& calibda) :
  TNamed(calibda)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int t=0; t<6; t++){
     fMeanTDC[t] = calibda.GetMeanTDC(t);
     fWidthTDC[t] = calibda.GetWidthTDC(t);
  }
}

//________________________________________________________________
AliZDCTDCCalib &AliZDCTDCCalib::operator =(const AliZDCTDCCalib& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int t=0; t<6; t++){
     fMeanTDC[t] = calibda.GetMeanTDC(t);
     fWidthTDC[t] = calibda.GetWidthTDC(t);
  }

  return *this;
}

//________________________________________________________________
AliZDCTDCCalib::~AliZDCTDCCalib()
{
}

//________________________________________________________________
void AliZDCTDCCalib::Reset()
{
  // Reset
  memset(fMeanTDC,0,6*sizeof(Float_t));
  memset(fWidthTDC,0,6*sizeof(Float_t));
}                                                                                       


//________________________________________________________________
void  AliZDCTDCCalib::Print(Option_t *) const
{
   // Printing of calibration object
   printf("\n ####### TDC calibration values ####### \n");
   for(int t=0; t<6; t++) 
      printf("\t ch.%d (%f, %f)\n",t,fMeanTDC[t],fWidthTDC[t]);
 
} 

//________________________________________________________________
void AliZDCTDCCalib::SetMeanTDC(Float_t* mean)
{
  if(mean) for(int t=0; t<6; t++) fMeanTDC[t] = mean[t];
  else for(int t=0; t<6; t++) fMeanTDC[t] = 0.;
}
//________________________________________________________________
void AliZDCTDCCalib::SetWidthTDC(Float_t* width)
{
  if(width) for(int t=0; t<6; t++) fWidthTDC[t] = width[t];
  else for(int t=0; t<6; t++) fWidthTDC[t] = 0.;
}
