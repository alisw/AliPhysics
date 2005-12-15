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

/* $Id: AliVZEROCalibData.cxx,                                            */

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// class for VZERO calibration                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "AliVZEROCalibData.h"

ClassImp(AliVZEROCalibData)

//________________________________________________________________
AliVZEROCalibData::AliVZEROCalibData()
{
  
}

//________________________________________________________________
AliVZEROCalibData::AliVZEROCalibData(const char* name)
{
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());

}

//________________________________________________________________
AliVZEROCalibData::AliVZEROCalibData(const AliVZEROCalibData& calibda) :
  TNamed(calibda)
{
// copy constructor

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  for(int t=0; t<80; t++) { 
      fPedestal[t] = calibda.GetPedestal(t);
      fGain[t]     = calibda.GetGain(t); }
  
}

//________________________________________________________________
AliVZEROCalibData &AliVZEROCalibData::operator =(const AliVZEROCalibData& calibda)
{
// assignment operator

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  for(int t=0; t<80; t++) {
      fPedestal[t] = calibda.GetPedestal(t);
      fGain[t]     = calibda.GetGain(t); }
  
  return *this;
  
}

//________________________________________________________________
AliVZEROCalibData::~AliVZEROCalibData()
{
  
}

                                                                                   

//________________________________________________________________
void AliVZEROCalibData::SetPedestal(Float_t* Pedestal)
{
  if(Pedestal) for(int t=0; t<80; t++) fPedestal[t] = Pedestal[t];
  else for(int t=0; t<80; t++) fPedestal[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetGain(Float_t* Gain) 
{
  if(Gain) for(int t=0; t<80; t++) fGain[t] = Gain[t];
  else for(int t=0; t<80; t++) fGain[t] = 0.0;
}

