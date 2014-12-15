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
// 		class for ZDC  calibration      			     //
// 		-> values for calibration 		                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCSaturationCalib.h"

ClassImp(AliZDCSaturationCalib)

//________________________________________________________________
AliZDCSaturationCalib::AliZDCSaturationCalib():
TNamed()
{
  Reset();
  for(Int_t i=0; i<4; i++){
     fZNASatCalibration[i] = 0.;
     fZNCSatCalibration[i] = 0.;
  }
}

//________________________________________________________________
AliZDCSaturationCalib::AliZDCSaturationCalib(const char* name):
TNamed()
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  for(Int_t i=0; i<4; i++){
     fZNASatCalibration[i] = 0.;
     fZNCSatCalibration[i] = 0.;
  }
}

//________________________________________________________________
AliZDCSaturationCalib::AliZDCSaturationCalib(const AliZDCSaturationCalib& calibda) :
  TNamed(calibda)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int i=0; i<4; i++){
    fZNASatCalibration[i] = calibda.GetZNASatCalib(i);
    fZNCSatCalibration[i] = calibda.GetZNCSatCalib(i);
  }
}

//________________________________________________________________
AliZDCSaturationCalib &AliZDCSaturationCalib::operator =(const AliZDCSaturationCalib& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int i=0; i<4; i++){
     fZNASatCalibration[i] = calibda.GetZNASatCalib(i);
     fZNCSatCalibration[i] = calibda.GetZNCSatCalib(i);
  }
  
  return *this;
}

//________________________________________________________________
AliZDCSaturationCalib::~AliZDCSaturationCalib()
{
}

//________________________________________________________________
void AliZDCSaturationCalib::Reset()
{
  // Reset
}                                                                                       


//________________________________________________________________
void  AliZDCSaturationCalib::Print(Option_t *) const
{
   // Printing of calibration object
   printf("\n\n ####### Parameters for ZN knee saturation coeffieicients #######	\n");
   printf("  ZNA => %e %e %e %e \n",fZNASatCalibration[0],fZNASatCalibration[1],fZNASatCalibration[2],fZNASatCalibration[3]);
   printf("  ZNC => %e %e %e %e \n",fZNCSatCalibration[0],fZNCSatCalibration[1],fZNCSatCalibration[2],fZNCSatCalibration[3]);

} 

//________________________________________________________________
void AliZDCSaturationCalib::SetZNASatCalib(Float_t* satCalib) 
{
  // Set energy calibration coefficients
  if(satCalib) for(int t=0; t<4; t++) fZNASatCalibration[t] = satCalib[t];
  else for(int t=0; t<4; t++) fZNASatCalibration[t] = 0.;
}

//________________________________________________________________
void AliZDCSaturationCalib::SetZNCSatCalib(Float_t* satCalib) 
{
  // Set energy calibration coefficients
  if(satCalib) for(int t=0; t<4; t++) fZNCSatCalibration[t] = satCalib[t];
  else for(int t=0; t<4; t++) fZNCSatCalibration[t] = 0.;
}
