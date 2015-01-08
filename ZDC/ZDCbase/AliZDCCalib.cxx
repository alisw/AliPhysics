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
// class for ZDC calibration      					     //
// -> values for energy calibration and relative sector calibration          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCCalib.h"

ClassImp(AliZDCCalib)

//________________________________________________________________
AliZDCCalib::AliZDCCalib():
TNamed()
{
  Reset();
}

//________________________________________________________________
AliZDCCalib::AliZDCCalib(const char* name):
TNamed()
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  for(Int_t i=0; i<6; i++){
    fEnCalibration[i] = 0.;
    if(i<5){
      fZN1EqualCoeff[i] = 0.;
      fZP1EqualCoeff[i] = 0.;
      fZN2EqualCoeff[i] = 0.;
      fZP2EqualCoeff[i] = 0.;
    }
  }
}

//________________________________________________________________
AliZDCCalib::AliZDCCalib(const AliZDCCalib& calibda) :
  TNamed(calibda)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int i=0; i<6; i++){
    fEnCalibration[i] = calibda.GetEnCalib(i);
    if(i<5){
      fZN1EqualCoeff[i] =  calibda.GetZN1EqualCoeff(i);
      fZP1EqualCoeff[i] =  calibda.GetZP1EqualCoeff(i);
      fZN2EqualCoeff[i] =  calibda.GetZN2EqualCoeff(i);
      fZP2EqualCoeff[i] =  calibda.GetZP2EqualCoeff(i);
    }
  }
}

//________________________________________________________________
AliZDCCalib &AliZDCCalib::operator =(const AliZDCCalib& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int i=0; i<6; i++){
     fEnCalibration[i] = calibda.GetEnCalib(i);
     if(i<5){
      fZN1EqualCoeff[i] =  calibda.GetZN1EqualCoeff(i);
      fZP1EqualCoeff[i] =  calibda.GetZP1EqualCoeff(i);
      fZN2EqualCoeff[i] =  calibda.GetZN2EqualCoeff(i);
      fZP2EqualCoeff[i] =  calibda.GetZP2EqualCoeff(i);
     }
  }
  
  return *this;
}

//________________________________________________________________
AliZDCCalib::~AliZDCCalib()
{
}

//________________________________________________________________
void AliZDCCalib::Reset()
{
  // Reset
  memset(fEnCalibration,0,6*sizeof(Float_t));
  memset(fZN1EqualCoeff,0,5*sizeof(Float_t));
  memset(fZP1EqualCoeff,0,5*sizeof(Float_t));
  memset(fZN2EqualCoeff,0,5*sizeof(Float_t));
  memset(fZP2EqualCoeff,0,5*sizeof(Float_t));
}                                                                                       


//________________________________________________________________
void  AliZDCCalib::Print(Option_t *) const
{
   // Printing of calibration object
   printf("\n\n ####### Energy calibration coefficients #######	\n");
   printf("  ZN1 = %.4f (E[TeV]/ADCch.) \n",fEnCalibration[0]);
   printf("  ZP1 = %.4f (E[TeV]/ADCch.) \n",fEnCalibration[1]);
   printf("  ZN2 = %.4f (E[TeV]/ADCch.) \n",fEnCalibration[2]);
   printf("  ZP2 = %.4f (E[TeV]/ADCch.) \n",fEnCalibration[3]);
   printf("  ZEM1 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[4]);
   printf("  ZEM2 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[5]);
 
   printf("\n\n ####### Equalization coefficients ####### \n");
   printf("  ZN1 -> %1.2f %1.2f %1.2f %1.2f %1.2f  \n",
    fZN1EqualCoeff[0],fZN1EqualCoeff[1],fZN1EqualCoeff[2],fZN1EqualCoeff[3],fZN1EqualCoeff[4]);
   printf("  ZP1 -> %1.2f %1.2f %1.2f %1.2f %1.2f  \n",
    fZP1EqualCoeff[0],fZP1EqualCoeff[1],fZP1EqualCoeff[2],fZP1EqualCoeff[3],fZP1EqualCoeff[4]);
   printf("  ZN2 -> %1.2f %1.2f %1.2f %1.2f %1.2f  \n",
    fZN2EqualCoeff[0],fZN2EqualCoeff[1],fZN2EqualCoeff[2],fZN2EqualCoeff[3],fZN2EqualCoeff[4]);
   printf("  ZP2 -> %1.2f %1.2f %1.2f %1.2f %1.2f  \n",
    fZP2EqualCoeff[0],fZP2EqualCoeff[1],fZP2EqualCoeff[2],fZP2EqualCoeff[3],fZP2EqualCoeff[4]);

} 

//________________________________________________________________
void AliZDCCalib::SetEnCalib(Float_t* EnCalib) 
{
  // Set energy calibration coefficients
  if(EnCalib) for(int t=0; t<6; t++) fEnCalibration[t] = EnCalib[t];
  else for(int t=0; t<6; t++) fEnCalibration[t] = 0.;
}

//________________________________________________________________
void AliZDCCalib::SetZN1EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZN1 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZN1EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZN1EqualCoeff[t] = 1.;
}
 
//________________________________________________________________
void AliZDCCalib::SetZP1EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZP1 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZP1EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZP1EqualCoeff[t] = 1.;
}
//________________________________________________________________
void AliZDCCalib::SetZN2EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZN2 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZN2EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZN2EqualCoeff[t] = 1.;
}
 
//________________________________________________________________
void AliZDCCalib::SetZP2EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZN1 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZP2EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZP2EqualCoeff[t] = 1.;
}
 
