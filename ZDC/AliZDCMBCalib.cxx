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

#include "AliZDCMBCalib.h"

ClassImp(AliZDCMBCalib)

//________________________________________________________________
AliZDCMBCalib::AliZDCMBCalib():
  TNamed(),
  fhZDCvsZEM(0x0),
  fhZDCCvsZEM(0x0),
  fhZDCAvsZEM(0x0)
{
  Reset();
}

//________________________________________________________________
AliZDCMBCalib::AliZDCMBCalib(const char* name):
  TNamed(),
  fhZDCvsZEM(0x0),
  fhZDCCvsZEM(0x0),
  fhZDCAvsZEM(0x0)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliZDCMBCalib::AliZDCMBCalib(const char* name, 
  	      TH2F *hzdcvszem, TH2F *hzdccvszem, TH2F *hzdcavszem):
  TNamed(),
  fhZDCvsZEM(hzdcvszem),
  fhZDCCvsZEM(hzdccvszem),
  fhZDCAvsZEM(hzdcavszem)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}
  
//________________________________________________________________
AliZDCMBCalib::AliZDCMBCalib(const AliZDCMBCalib& calibda) :
  TNamed(calibda),
  fhZDCvsZEM(0), 
  fhZDCCvsZEM(0),
  fhZDCAvsZEM(0)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  fhZDCvsZEM = calibda.GethZDCvsZEM(); 
  fhZDCCvsZEM = calibda.GethZDCCvsZEM();
  fhZDCAvsZEM = calibda.GethZDCAvsZEM();


}

//________________________________________________________________
AliZDCMBCalib &AliZDCMBCalib::operator =(const AliZDCMBCalib& calibda)
{
// assignment operator
  if(&calibda == this) return *this;
  
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  fhZDCvsZEM = calibda.GethZDCvsZEM(); 
  fhZDCCvsZEM = calibda.GethZDCCvsZEM();
  fhZDCAvsZEM = calibda.GethZDCAvsZEM();
  
  return *this;
}

//________________________________________________________________
AliZDCMBCalib::~AliZDCMBCalib()
{
}

//________________________________________________________________
void AliZDCMBCalib::Reset()
{
  // Reset
  
}                                                                                       

