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
  for(Int_t i=0; i<2; i++){
    fReferenceValue[i] = 0.;
    fSector[i] = 0;
    fGain[i] = 0;
    fPMRefValue[i] = 0.;
    fPMRefWidth[i] = 0.;
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
  for(int t=0; t<2; t++){
     fReferenceValue[t] = calibda.GetReferenceValue(t);
     fSector[t] = calibda.GetSector(t);
     fGain[t] = calibda.GetGain(t);
     fPMRefValue[t] = calibda.GetPMRefValue(t);
     fPMRefWidth[t] = calibda.GetPMRefWidth(t);
  }
}

//________________________________________________________________
AliZDCLaserCalib &AliZDCLaserCalib::operator =(const AliZDCLaserCalib& calibda) 
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int t=0; t<2; t++){
     fReferenceValue[t] = calibda.GetReferenceValue(t);
     fSector[t]     = calibda.GetSector(t);
     fGain[t] = calibda.GetGain(t);
     fPMRefValue[t] = calibda.GetPMRefValue(t);
     fPMRefWidth[t] = calibda.GetPMRefWidth(t);
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
  memset(fReferenceValue,0,2*sizeof(Float_t));
  memset(fSector,0,2*sizeof(Float_t));
  memset(fGain,0,2*sizeof(Float_t));
  memset(fPMRefValue,0,2*sizeof(Float_t));
  memset(fPMRefWidth,0,2*sizeof(Float_t));
}                                                                                       


//________________________________________________________________
void  AliZDCLaserCalib::Print(Option_t *) const
{
   // Printing of calibration object
   printf("\n\n ******************* AliZDCLaserCalib object *******************\n");
   for(Int_t i=0; i<2; i++){
     if(fSector[i]==1.) printf("  side C ->");
     else if(fSector[i]==4.) printf("  side A->");
     printf("  Ref.Values %1.2f  Gain %1.0f",fReferenceValue[i], fGain[i]);
     printf(" PMRefValue = %1.0f  PMRefWidth = %1.0f\n",fPMRefValue[i],fPMRefWidth[i]);
   }
   printf(" ***************************************************************\n\n");
} 
