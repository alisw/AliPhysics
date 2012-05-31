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

/* $Id: AliZDCTowerCalib.cxx 22045 2007-11-08 13:31:24Z coppedis $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ZDC calibration      					     //
// -> values for energy calibration and relative sector calibration          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCTowerCalib.h"

ClassImp(AliZDCTowerCalib)

//________________________________________________________________
AliZDCTowerCalib::AliZDCTowerCalib():
TNamed()
{
  Reset();
}

//________________________________________________________________
AliZDCTowerCalib::AliZDCTowerCalib(const char* name):
TNamed()
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  for(Int_t i=0; i<5; i++){
    fZN1EqualCoeff[i] = 0.;
    fZP1EqualCoeff[i] = 0.;
    fZN2EqualCoeff[i] = 0.;
    fZP2EqualCoeff[i] = 0.;
  }
}

//________________________________________________________________
AliZDCTowerCalib::AliZDCTowerCalib(const AliZDCTowerCalib& calibda) :
  TNamed(calibda)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int i=0; i<5; i++){
     fZN1EqualCoeff[i] =  calibda.GetZN1EqualCoeff(i);
     fZP1EqualCoeff[i] =  calibda.GetZP1EqualCoeff(i);
     fZN2EqualCoeff[i] =  calibda.GetZN2EqualCoeff(i);
     fZP2EqualCoeff[i] =  calibda.GetZP2EqualCoeff(i);
  }
}

//________________________________________________________________
AliZDCTowerCalib &AliZDCTowerCalib::operator =(const AliZDCTowerCalib& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int i=0; i<5; i++){
     fZN1EqualCoeff[i] =  calibda.GetZN1EqualCoeff(i);
     fZP1EqualCoeff[i] =  calibda.GetZP1EqualCoeff(i);
     fZN2EqualCoeff[i] =  calibda.GetZN2EqualCoeff(i);
     fZP2EqualCoeff[i] =  calibda.GetZP2EqualCoeff(i);
  }
  
  return *this;
}

//________________________________________________________________
AliZDCTowerCalib::~AliZDCTowerCalib()
{
}

//________________________________________________________________
void AliZDCTowerCalib::Reset()
{
  // Reset
  memset(fZN1EqualCoeff,0,5*sizeof(Float_t));
  memset(fZP1EqualCoeff,0,5*sizeof(Float_t));
  memset(fZN2EqualCoeff,0,5*sizeof(Float_t));
  memset(fZP2EqualCoeff,0,5*sizeof(Float_t));
}                                                                                       


//________________________________________________________________
void  AliZDCTowerCalib::Print(Option_t *) const
{
   // Printing of calibration object
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
void AliZDCTowerCalib::SetZN1EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZN1 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZN1EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZN1EqualCoeff[t] = 1.;
}
 
//________________________________________________________________
void AliZDCTowerCalib::SetZP1EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZP1 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZP1EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZP1EqualCoeff[t] = 1.;
}
//________________________________________________________________
void AliZDCTowerCalib::SetZN2EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZN2 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZN2EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZN2EqualCoeff[t] = 1.;
}
 
//________________________________________________________________
void AliZDCTowerCalib::SetZP2EqualCoeff(Float_t* EqualCoeff)
{
  // Set ZN1 equalization coefficients
  if(EqualCoeff) for(int t=0; t<5; t++) fZP2EqualCoeff[t] = EqualCoeff[t];
  else for(int t=0; t<5; t++) fZP2EqualCoeff[t] = 1.;
}
 
