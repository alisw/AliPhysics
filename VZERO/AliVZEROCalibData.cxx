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
  // 
}

//________________________________________________________________
void AliVZEROCalibData::Reset()
{
  // reset 
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
  
  for(int t=0; t<128; t++) { 
      fPedestal[t] = calibda.GetPedestal(t);
      fSigma[t]    = calibda.GetSigma(t);
      fGain[t]     = calibda.GetGain(t); }
      
  for(int t=0; t<64; t++) { 
      fTimeOffset[t]   = calibda.GetTimeOffset(t);
      fTimeGain[t]     = calibda.GetTimeGain(t); 
      fMeanHV[t]       = calibda.GetMeanHV(t);
      fWidthHV[t]      = calibda.GetWidthHV(t); }   
  
}

//________________________________________________________________
AliVZEROCalibData &AliVZEROCalibData::operator =(const AliVZEROCalibData& calibda)
{
// assignment operator

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  for(int t=0; t<128; t++) {
      fPedestal[t] = calibda.GetPedestal(t);
      fSigma[t]    = calibda.GetSigma(t);
      fGain[t]     = calibda.GetGain(t); }
      
  for(int t=0; t<64; t++) { 
      fTimeOffset[t]   = calibda.GetTimeOffset(t);
      fTimeGain[t]     = calibda.GetTimeGain(t); 
      fMeanHV[t]       = calibda.GetMeanHV(t);
      fWidthHV[t]      = calibda.GetWidthHV(t); }     
  
  return *this;
  
}

//________________________________________________________________
AliVZEROCalibData::~AliVZEROCalibData()
{
  // destructor
}

                                                                                   

//________________________________________________________________
void AliVZEROCalibData::SetPedestal(Float_t* Pedestal)
{
  if(Pedestal) for(int t=0; t<128; t++) fPedestal[t] = Pedestal[t];
  else for(int t=0; t<128; t++) fPedestal[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetSigma(Float_t* Sigma)
{
  if(Sigma) for(int t=0; t<128; t++) fSigma[t] = Sigma[t];
  else for(int t=0; t<128; t++) fSigma[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetGain(Float_t* Gain) 
{
  if(Gain) for(int t=0; t<128; t++) fGain[t] = Gain[t];
  else for(int t=0; t<128; t++) fGain[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetTimeOffset(Float_t* TimeOffset) 
{
  if(TimeOffset) for(int t=0; t<64; t++) fTimeOffset[t] = TimeOffset[t];
  else for(int t=0; t<64; t++) fTimeOffset[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetTimeGain(Float_t* TimeGain) 
{
  if(TimeGain) for(int t=0; t<64; t++) fTimeGain[t] = TimeGain[t];
  else for(int t=0; t<64; t++) fTimeGain[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetMeanHV(Float_t* MeanHVs) 
{
  if(MeanHVs) for(int t=0; t<64; t++) fMeanHV[t] = MeanHVs[t];
  else for(int t=0; t<64; t++) fMeanHV[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetWidthHV(Float_t* WidthHVs) 
{
  if(WidthHVs) for(int t=0; t<64; t++) fWidthHV[t] = WidthHVs[t];
  else for(int t=0; t<64; t++) fWidthHV[t] = 0.0;
}
