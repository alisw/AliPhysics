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
// class for ZDC calibration      -> values for pedestal subtraction         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCPedestals.h"

ClassImp(AliZDCPedestals)

//________________________________________________________________
AliZDCPedestals::AliZDCPedestals():
TNamed()
{
  Reset();
  SetPedModeBit(kFALSE);
}

//________________________________________________________________
AliZDCPedestals::AliZDCPedestals(const char* name):
TNamed(),
fPedSubModefromOCDB(0)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  for(Int_t i=0; i<48; i++){
    fMeanPedestal[i] = 0.;
    fMeanPedWidth[i] = 0.;
    fOOTPedestal[i] = 0.;
    fOOTPedWidth[i] = 0.;
    for(Int_t j=0; j<2; j++) fPedCorrCoeff[j][i] = 0.;
  }
  for(Int_t i=0; i<24; i++) fUseCorrFit[i] = kFALSE;
  SetPedModeBit(kFALSE);
  
}

//________________________________________________________________
AliZDCPedestals::AliZDCPedestals(const AliZDCPedestals& calibda) :
  TNamed(calibda)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int t=0; t<48; t++){
     fMeanPedestal[t] = calibda.GetMeanPed(t);
     fMeanPedWidth[t] = calibda.GetMeanPedWidth(t);
     fOOTPedestal[t]  = calibda.GetOOTPed(t);
     fOOTPedWidth[t]  = calibda.GetOOTPedWidth(t);
     fPedCorrCoeff[0][t] = calibda.GetPedCorrCoeff0(t);
     fPedCorrCoeff[1][t] = calibda.GetPedCorrCoeff1(t);
  }
  for(Int_t i=0; i<24; i++) fUseCorrFit[i] = calibda.GetUseCorrFit(i);
}

//________________________________________________________________
AliZDCPedestals &AliZDCPedestals::operator =(const AliZDCPedestals& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int t=0; t<48; t++){
     fMeanPedestal[t] = calibda.GetMeanPed(t);
     fMeanPedWidth[t] = calibda.GetMeanPedWidth(t);
     fOOTPedestal[t]  = calibda.GetOOTPed(t);
     fOOTPedWidth[t]  = calibda.GetOOTPedWidth(t);
     fPedCorrCoeff[0][t] = calibda.GetPedCorrCoeff0(t);
     fPedCorrCoeff[1][t] = calibda.GetPedCorrCoeff1(t);
  }
  for(Int_t i=0; i<24; i++) fUseCorrFit[i] = calibda.GetUseCorrFit(i);

  return *this;
}

//________________________________________________________________
AliZDCPedestals::~AliZDCPedestals()
{
}

//________________________________________________________________
void AliZDCPedestals::Reset()
{
  // Reset
  memset(fMeanPedestal,0,48*sizeof(Float_t));
  memset(fMeanPedWidth,0,48*sizeof(Float_t));
  memset(fOOTPedestal,0,48*sizeof(Float_t));
  memset(fOOTPedWidth,0,48*sizeof(Float_t));
}                                                                                       


//________________________________________________________________
void  AliZDCPedestals::Print(Option_t *) const
{
   // Printing of calibration object
   printf(" \n  Pedestal subtraction mode bit %d\n",AliZDCPedestals::TestPedModeBit());
   printf("\n ####### In-time pedestal values (mean value, sigma) ####### \n");
   for(int t=0; t<48; t++){ 
      if(t<24) printf("\t ADC%d (%.1f, %.1f) ped.sub.mode %d\n",t,fMeanPedestal[t],fMeanPedWidth[t],fUseCorrFit[t]);
      else printf("\t ADC%d (%.1f, %.1f) \n",t,fMeanPedestal[t],fMeanPedWidth[t]);
   }
   //
   /*if(AliZDCPedestals::TestPedModeBit()){
     printf("\n\n ####### Out-of-time pedestal values (mean value, sigma) ####### \n");
     for(int t=0; t<48; t++)
       printf("\t ADC-OoT%d (%.1f, %.1f)\n",t,fOOTPedestal[t],fOOTPedWidth[t]);
   }*/
} 

//________________________________________________________________
void AliZDCPedestals::SetMeanPed(Float_t* MeanPed)
{
  if(MeanPed) for(int t=0; t<48; t++) fMeanPedestal[t] = MeanPed[t];
  else for(int t=0; t<48; t++) fMeanPedestal[t] = 0.;
}
//________________________________________________________________
void AliZDCPedestals::SetMeanPedWidth(Float_t* MeanPedWidth)
{
  if(MeanPedWidth) for(int t=0; t<48; t++) fMeanPedWidth[t] = MeanPedWidth[t];
  else for(int t=0; t<48; t++) fMeanPedWidth[t] = 0.;
}

//________________________________________________________________
void AliZDCPedestals::SetOOTPed(Float_t* OOTPed)
{
  if(OOTPed) for(int t=0; t<48; t++) fOOTPedestal[t] = OOTPed[t];
  else for(int t=0; t<48; t++) fOOTPedestal[t] = 0.;
}

//________________________________________________________________
void AliZDCPedestals::SetOOTPedWidth(Float_t* OOTPedWidth)
{
  if(OOTPedWidth) for(int t=0; t<48; t++) fOOTPedWidth[t] = OOTPedWidth[t];
  else for(int t=0; t<48; t++) fOOTPedWidth[t] = 0.;
}

//________________________________________________________________
void AliZDCPedestals:: SetPedCorrCoeff(Float_t* PedCorrCoeff)
{
  // Set coefficients for pedestal correlations
  if(PedCorrCoeff){
    for(Int_t j=0; j<2; j++){
     for(int t=0; t<48; t++)
       fPedCorrCoeff[j][t] = PedCorrCoeff[t];
    }
  }
  else{
    for(Int_t j=0; j<2; j++){
     for(int t=0; t<48; t++)
       fPedCorrCoeff[j][t] = 0.;
    }
  }
 
}

//________________________________________________________________
void AliZDCPedestals:: SetPedCorrCoeff(Float_t* PedCorrCoeff0, Float_t* PedCorrCoeff1)
{
  // Set coefficients for pedestal correlations
  if(PedCorrCoeff0 && PedCorrCoeff1){
    for(int t=0; t<48; t++){
       fPedCorrCoeff[0][t] = PedCorrCoeff0[t];
       fPedCorrCoeff[1][t] = PedCorrCoeff1[t];
    }
  }
  else{
     for(int t=0; t<48; t++){
       fPedCorrCoeff[0][t] = 0.;
       fPedCorrCoeff[1][t] = 0.;
    }
  }
 
}
