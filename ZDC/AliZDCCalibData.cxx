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

#include "AliZDCCalibData.h"

ClassImp(AliZDCCalibData)

//________________________________________________________________
AliZDCCalibData::AliZDCCalibData()
{
  fHistMeanPed=0;
  Reset();
}

//________________________________________________________________
AliZDCCalibData::AliZDCCalibData(const char* name)
{
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  fHistMeanPed=0;
  Reset();
}

//________________________________________________________________
AliZDCCalibData::~AliZDCCalibData()
{
  CleanHistos();
}

//________________________________________________________________
void AliZDCCalibData::Reset()
{
  memset(fMeanPedestal,0,47*sizeof(Float_t));
  memset(fEnCalibration,0,4*sizeof(Float_t));
}                                                                                       

//________________________________________________________________
void AliZDCCalibData::CleanHistos()
{
  if (fHistMeanPed) delete fHistMeanPed; fHistMeanPed = 0;
}

//________________________________________________________________
void AliZDCCalibData::PrepHistos()
{
  CleanHistos();
  Int_t   kNChannels = 100;
  Float_t kMaxPedVal = 400.;
  TString hname = GetName();  hname += "_Pedestals";
  fHistMeanPed   = new TH1F(hname.Data(),hname.Data(),kNChannels,0.,kMaxPedVal);
  for (int i=0; i<64; i++) {
    fHistMeanPed->SetBinContent(i+1,GetMeanPed(i));
  }
}

//________________________________________________________________
void  AliZDCCalibData::Print() const
{
   printf("\n	----	Mean pedestal values	----\n\n");
   for(int t=0; t<47; t++){
      printf("	MeanPed Value[ADC%d] = %f\n",t,fMeanPedestal[t]);
   }
 
   printf("\n	----	Energy calibration coefficients 	----\n\n");
   for(int t=0; t<4; t++){
      printf("	En Calib Coeff. [ZDC%d] = %f\n",t,fEnCalibration[t]);
   }
} 

//________________________________________________________________
void AliZDCCalibData::SetMeanPed(Float_t* MeanPed)
{
  if(MeanPed) for(int t=0; t<47; t++) fMeanPedestal[t] = MeanPed[t];
  else for(int t=0; t<43; t++) fMeanPedestal[t] = 0.;
}

//________________________________________________________________
void AliZDCCalibData::SetEnCalib(Float_t* EnCalib) 
{
  if(EnCalib) for(int t=0; t<4; t++) fEnCalibration[t] = EnCalib[t];
  else for(int t=0; t<4; t++) fEnCalibration[t] = 0.;
}

