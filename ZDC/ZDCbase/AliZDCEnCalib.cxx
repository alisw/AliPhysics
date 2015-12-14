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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// 		class for ZDC ENERGY calibration      			     //
// 		-> values for energy calibration 		             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCEnCalib.h"

ClassImp(AliZDCEnCalib)

//________________________________________________________________
AliZDCEnCalib::AliZDCEnCalib():
TNamed()
{
  Reset();
  for(Int_t i=0; i<6; i++) fEnCalibration[i] = 0.;
}

//________________________________________________________________
AliZDCEnCalib::AliZDCEnCalib(const char* name):
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
  }
}

//________________________________________________________________
AliZDCEnCalib::AliZDCEnCalib(const AliZDCEnCalib& calibda) :
  TNamed(calibda)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int i=0; i<6; i++){
    fEnCalibration[i] = calibda.GetEnCalib(i);
  }
}

//________________________________________________________________
AliZDCEnCalib &AliZDCEnCalib::operator =(const AliZDCEnCalib& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int i=0; i<6; i++){
     fEnCalibration[i] = calibda.GetEnCalib(i);
  }
  
  return *this;
}

//________________________________________________________________
AliZDCEnCalib::~AliZDCEnCalib()
{
}

//________________________________________________________________
void AliZDCEnCalib::Reset()
{
  // Reset
}                                                                                       


//________________________________________________________________
void  AliZDCEnCalib::Print(Option_t *) const
{
   // Printing of calibration object
   printf("\n\n ####### Energy calibration coefficients #######	\n");
   printf("  ZNC = %.4f (E[TeV]/ADCch.) \n",fEnCalibration[0]);
   printf("  ZPC = %.4f (E[TeV]/ADCch.) \n",fEnCalibration[1]);
   printf("  ZNA = %.4f (E[TeV]/ADCch.) \n",fEnCalibration[2]);
   printf("  ZPA = %.4f (E[TeV]/ADCch.) \n",fEnCalibration[3]);
   printf("  ZEM1 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[4]);
   printf("  ZEM2 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[5]);

} 

//________________________________________________________________
void AliZDCEnCalib::SetEnCalib(Float_t* EnCalib) 
{
  // Set energy calibration coefficients
  if(EnCalib) for(int t=0; t<6; t++) fEnCalibration[t] = EnCalib[t];
  else for(int t=0; t<6; t++) fEnCalibration[t] = 0.;
}
