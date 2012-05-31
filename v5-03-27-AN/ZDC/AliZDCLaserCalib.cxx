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

///////////////////////////////////////////////////////////////
//                                                           //
// Class for ZDC calibration -> LASER signal monitoring      //
// author: Chiara Oppedisano			             //
//                                                           //
///////////////////////////////////////////////////////////////

#include "AliZDCLaserCalib.h"

ClassImp(AliZDCLaserCalib)

//________________________________________________________________
AliZDCLaserCalib::AliZDCLaserCalib():
TNamed()
{
  Reset();
}

//________________________________________________________________
AliZDCLaserCalib::AliZDCLaserCalib(const char* name):
TNamed()
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  for(Int_t i=0; i<24; i++){
    fDetector[i] = 0;
    fSector[i] = 0;
    fPMValue[i] = 0.;
    fPMWidth[i] = 0.;
  }
  
  
}

//________________________________________________________________
AliZDCLaserCalib::AliZDCLaserCalib(const AliZDCLaserCalib& calibda) :
  TNamed(calibda)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int t=0; t<22; t++){
     fDetector[t] = calibda.GetDetector(t);
     fSector[t]   = calibda.GetSector(t);
     fPMValue[t]  = calibda.GetPMValue(t);
     fPMWidth[t]  = calibda.GetPMWidth(t);
  }
}

//________________________________________________________________
AliZDCLaserCalib &AliZDCLaserCalib::operator =(const AliZDCLaserCalib& calibda) 
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int t=0; t<22; t++){
     fDetector[t] = calibda.GetDetector(t);
     fSector[t]   = calibda.GetSector(t);
     fPMValue[t]  = calibda.GetPMValue(t);
     fPMWidth[t]  = calibda.GetPMWidth(t);
  }

  return *this;
}

//________________________________________________________________
AliZDCLaserCalib::~AliZDCLaserCalib()
{
}

//________________________________________________________________
void AliZDCLaserCalib::Reset()
{
  // Reset
  memset(fDetector,0,24*sizeof(Float_t));
  memset(fSector,0,24*sizeof(Float_t));
  memset(fPMValue,0,24*sizeof(Float_t));
  memset(fPMWidth,0,24*sizeof(Float_t));
}                                                                                       


//________________________________________________________________
void  AliZDCLaserCalib::Print(Option_t *) const
{
   // Printing of calibration object
   printf("\n\n ******************* AliZDCLaserCalib object *******************\n");
   for(Int_t i=0; i<24; i++){
     printf("  Detector %d, sector %d, PMValue = %1.1f +/- %1.1f\n",
     	fDetector[i], fSector[i],fPMValue[i],fPMWidth[i]);
   }
   printf(" ***************************************************************\n\n");
} 
