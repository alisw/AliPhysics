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
// class for ZDC calibration                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCRecParam.h"

ClassImp(AliZDCRecParam)

//________________________________________________________________
AliZDCRecParam::AliZDCRecParam():
TNamed()
{
  Reset();
}

//________________________________________________________________
AliZDCRecParam::AliZDCRecParam(const char* name):
TNamed(),
fZEMEndValue(0),
fZEMCutFraction(0),
fDZEMSup(0),
fDZEMInf(0),
fEZN1MaxValue(0),
fEZP1MaxValue(0),
fEZDC1MaxValue(0),
fEZN2MaxValue(0),
fEZP2MaxValue(0),
fEZDC2MaxValue(0)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliZDCRecParam::AliZDCRecParam(const AliZDCRecParam& calibda) :
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
  for(int t=0; t<6; t++)  fEnCalibration[t] = calibda.GetEnCalib(t);
  //
  fZEMEndValue    = calibda.GetZEMEndValue();   
  fZEMCutFraction = calibda.GetZEMCutFraction();
  fDZEMSup	  = calibda.GetDZEMSup();
  fDZEMInf	  = calibda.GetDZEMInf();
}

//________________________________________________________________
AliZDCRecParam &AliZDCRecParam::operator =(const AliZDCRecParam& calibda)
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
  for(int t=0; t<6; t++) fEnCalibration[t] = calibda.GetEnCalib(t);
  fZEMEndValue    = calibda.GetZEMEndValue();
  fZEMCutFraction = calibda.GetZEMCutFraction();

  return *this;
}

//________________________________________________________________
AliZDCRecParam::~AliZDCRecParam()
{
}

//________________________________________________________________
void AliZDCRecParam::Reset()
{
  // Reset
  memset(fMeanPedestal,0,48*sizeof(Float_t));
  memset(fMeanPedWidth,0,48*sizeof(Float_t));
  memset(fOOTPedestal,0,48*sizeof(Float_t));
  memset(fOOTPedWidth,0,48*sizeof(Float_t));
  memset(fEnCalibration,0,6*sizeof(Float_t));
  memset(fZN1EqualCoeff,0,5*sizeof(Float_t));
  memset(fZP1EqualCoeff,0,5*sizeof(Float_t));
  memset(fZN2EqualCoeff,0,5*sizeof(Float_t));
  memset(fZP2EqualCoeff,0,5*sizeof(Float_t));
}                                                                                       


//________________________________________________________________
void  AliZDCRecParam::Print(Option_t *) const
{
   // Printing calibration object
   printf("\n\n ####### Parameters from EZDC vs. ZEM correlation #######	\n");
   printf("  ZEMEndPoint = %1.2f, ZEMCutFraction = %1.2f \n"
     "  DZEMInf = %1.2f, DZEMSup = %1.2f\n",
     fZEMEndValue, fZEMCutFraction, fDZEMInf, fDZEMSup);
 
   printf("\n\n ####### Parameters from EZDC vs. Nspec correlation #######	\n");
   printf("  EZN1MaxValue = %1.2f, EZP1MaxValue = %1.2f, EZDC1MaxValue = %1.2f \n"
     "  EZN2MaxValue = %1.2f, EZP2MaxValue = %1.2f, EZDC2MaxValue = %1.2f \n\n",
     fEZN1MaxValue, fEZP1MaxValue, fEZDC1MaxValue,
     fEZN2MaxValue, fEZP2MaxValue, fEZDC2MaxValue);

} 

//________________________________________________________________
void AliZDCRecParam::SetMeanPed(Float_t* MeanPed)
{
  if(MeanPed) for(int t=0; t<48; t++) fMeanPedestal[t] = MeanPed[t];
  else for(int t=0; t<48; t++) fMeanPedestal[t] = 0.;
}
//________________________________________________________________
void AliZDCRecParam::SetMeanPedWidth(Float_t* MeanPedWidth)
{
  if(MeanPedWidth) for(int t=0; t<48; t++) fMeanPedWidth[t] = MeanPedWidth[t];
  else for(int t=0; t<48; t++) fMeanPedWidth[t] = 0.;
}

//________________________________________________________________
void AliZDCRecParam::SetOOTPed(Float_t* OOTPed)
{
  if(OOTPed) for(int t=0; t<48; t++) fOOTPedestal[t] = OOTPed[t];
  else for(int t=0; t<48; t++) fOOTPedestal[t] = 0.;
}

//________________________________________________________________
void AliZDCRecParam::SetOOTPedWidth(Float_t* OOTPedWidth)
{
  if(OOTPedWidth) for(int t=0; t<48; t++) fOOTPedWidth[t] = OOTPedWidth[t];
  else for(int t=0; t<48; t++) fOOTPedWidth[t] = 0.;
}

//________________________________________________________________
void AliZDCRecParam:: SetPedCorrCoeff(Float_t* PedCorrCoeff)
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
void AliZDCRecParam:: SetPedCorrCoeff(Float_t* PedCorrCoeff0, Float_t* PedCorrCoeff1)
{
  // Set coefficients for pedestal correlations
  if(PedCorrCoeff0 && PedCorrCoeff1){
    for(int t=0; t<48; t++){
       fPedCorrCoeff[0][t] = PedCorrCoeff0[t];
       fPedCorrCoeff[0][t] = PedCorrCoeff1[t];
    }
  }
  else{
     for(int t=0; t<48; t++){
       fPedCorrCoeff[0][t] = 0.;
       fPedCorrCoeff[1][t] = 0.;
    }
  }
 
}

//________________________________________________________________
void AliZDCRecParam::SetEnCalib(Float_t* EnCalib) 
{
  // Set energy calibration coefficients
  if(EnCalib) for(int t=0; t<6; t++) fEnCalibration[t] = EnCalib[t];
  else for(int t=0; t<6; t++) fEnCalibration[t] = 0.;
}

//________________________________________________________________
void AliZDCRecParam::SetZN1EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZN1 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZN1EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZN1EqualCoeff[t] = 1.;
}
 
//________________________________________________________________
void AliZDCRecParam::SetZP1EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZP1 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZP1EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZP1EqualCoeff[t] = 1.;
}
//________________________________________________________________
void AliZDCRecParam::SetZN2EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZN2 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZN2EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZN2EqualCoeff[t] = 1.;
}
 
//________________________________________________________________
void AliZDCRecParam::SetZP2EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZN1 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZP2EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZP2EqualCoeff[t] = 1.;
}
 
