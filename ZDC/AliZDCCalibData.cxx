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
//  fHistMeanPed=0;
  Reset();
}

//________________________________________________________________
AliZDCCalibData::AliZDCCalibData(const char* name)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
//  fHistMeanPed=0;
  Reset();
}

//________________________________________________________________
AliZDCCalibData::AliZDCCalibData(const AliZDCCalibData& calibda) :
  TNamed(calibda)
{
// copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int t=0; t<47; t++){
     fMeanPedestal[t] = calibda.GetMeanPed(t);
     fMeanPedWidth[t] = calibda.GetMeanPedWidth(t);
  }
  for(int t=0; t<44; t++){
     fOOTPedestal[t]  = calibda.GetOOTPed(t);
     fOOTPedWidth[t]  = calibda.GetOOTPedWidth(t);
     fPedCorrCoeff[0][t] = calibda.GetPedCorrCoeff0(t);
     fPedCorrCoeff[1][t] = calibda.GetPedCorrCoeff1(t);
  }
  for(int t=0; t<6; t++)  fEnCalibration[t] = calibda.GetEnCalib(t);
//  PrepHistos();
}

//________________________________________________________________
AliZDCCalibData &AliZDCCalibData::operator =(const AliZDCCalibData& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(int t=0; t<47; t++){
     fMeanPedestal[t] = calibda.GetMeanPed(t);
     fMeanPedWidth[t] = calibda.GetMeanPedWidth(t);
  }
  for(int t=0; t<44; t++){
     fOOTPedestal[t]  = calibda.GetOOTPed(t);
     fOOTPedWidth[t]  = calibda.GetOOTPedWidth(t);
     fPedCorrCoeff[0][t] = calibda.GetPedCorrCoeff0(t);
     fPedCorrCoeff[1][t] = calibda.GetPedCorrCoeff1(t);
  }
  for(int t=0; t<6; t++) fEnCalibration[t] = calibda.GetEnCalib(t);
//  PrepHistos();
  return *this;
}

//________________________________________________________________
AliZDCCalibData::~AliZDCCalibData()
{
//  CleanHistos();
}

//________________________________________________________________
void AliZDCCalibData::Reset()
{
  // Reset
  memset(fMeanPedestal,0,47*sizeof(Float_t));
  memset(fMeanPedWidth,0,47*sizeof(Float_t));
  memset(fOOTPedestal,0,44*sizeof(Float_t));
  memset(fOOTPedWidth,0,44*sizeof(Float_t));
  memset(fEnCalibration,0,6*sizeof(Float_t));
}                                                                                       

/*
//________________________________________________________________
void AliZDCCalibData::CleanHistos()
{
  if (fHistMeanPed) delete fHistMeanPed; fHistMeanPed = 0;
}

//________________________________________________________________
void AliZDCCalibData::PrepHistos()
{
//  CleanHistos(); // this gives a segm.viol!
  Int_t   kNChannels = 47;
  Float_t kMaxPedVal = 47.;
  TString hname = GetName();  hname += "_Pedestals";
  fHistMeanPed = new TH1F(hname.Data(),hname.Data(),kNChannels,0.,kMaxPedVal);
  for(int i=0; i<47; i++)  fHistMeanPed->SetBinContent(i+1,GetMeanPed(i));
}
*/

//________________________________________________________________
void  AliZDCCalibData::Print(Option_t *) const
{
   // Printing of calibration object
   printf("\n	#######	Mean pedestal values	####### \n");
   for(int t=0; t<47; t++){
     if(t==0 || t==24) printf("\n\t -------- ZN HighRes -------- \n");
     else if(t==5 || t==29) printf("\n\t -------- ZN LowRes -------- \n");
     else if(t==10 || t==34) printf("\n\t -------- ZP HighRes -------- \n");
     else if(t==15 || t==39) printf("\n\t -------- ZP LowRes -------- \n");
     else if(t==20) printf("\n\t -------- ZEM1 HighRes -------- ");
     else if(t==21) printf("\n\t -------- ZEM1 LowRes -------- ");
     else if(t==22) printf("\n\t -------- ZEM2 HighRes -------- ");
     else if(t==23) printf("\n\t -------- ZEM2 LowRes -------- ");
     printf("   MeanPed[ADC%d] = %.1f   ",t,fMeanPedestal[t]);
   }
   
   printf("\n\n\n	#######	Out Of Time pedestal values	####### \n");
   for(int t=0; t<44; t++){
     if(t==0 || t==24) printf("\n\t -------- ZN HighRes -------- \n");
     else if(t==5 || t==29) printf("\n\t -------- ZN LowRes -------- \n");
     else if(t==10 || t==34) printf("\n\t -------- ZP HighRes -------- \n");
     else if(t==15 || t==39) printf("\n\t -------- ZP LowRes -------- \n");
     else if(t==20) printf("\n\t -------- ZEM1 HighRes -------- ");
     else if(t==21) printf("\n\t -------- ZEM1 LowRes -------- ");
     else if(t==22) printf("\n\t -------- ZEM2 HighRes -------- ");
     else if(t==23) printf("\n\t -------- ZEM2 LowRes -------- ");
     printf("   OOTPed[ADC%d] = %.1f   ",t,fOOTPedestal[t]);
   }
 
   printf("\n\n\n	#######	Energy calibration coefficients #######	\n");
   printf("	ZN1 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[0]);
   printf("	ZP1 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[1]);
   printf("	ZN2 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[2]);
   printf("	ZP2 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[3]);
   printf("	ZEM1 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[4]);
   printf("	ZEM2 = %.2f (E[TeV]/ADCch.) \n",fEnCalibration[5]);
} 

//________________________________________________________________
void AliZDCCalibData::SetMeanPed(Float_t* MeanPed)
{
  if(MeanPed) for(int t=0; t<47; t++) fMeanPedestal[t] = MeanPed[t];
  else for(int t=0; t<47; t++) fMeanPedestal[t] = 0.;
}
//________________________________________________________________
void AliZDCCalibData::SetMeanPedWidth(Float_t* MeanPedWidth)
{
  if(MeanPedWidth) for(int t=0; t<47; t++) fMeanPedWidth[t] = MeanPedWidth[t];
  else for(int t=0; t<47; t++) fMeanPedWidth[t] = 0.;
}

//________________________________________________________________
void AliZDCCalibData::SetOOTPed(Float_t* OOTPed)
{
  if(OOTPed) for(int t=0; t<44; t++) fOOTPedestal[t] = OOTPed[t];
  else for(int t=0; t<44; t++) fOOTPedestal[t] = 0.;
}

//________________________________________________________________
void AliZDCCalibData::SetOOTPedWidth(Float_t* OOTPedWidth)
{
  if(OOTPedWidth) for(int t=0; t<44; t++) fOOTPedWidth[t] = OOTPedWidth[t];
  else for(int t=0; t<44; t++) fOOTPedWidth[t] = 0.;
}

//________________________________________________________________
void AliZDCCalibData:: SetPedCorrCoeff(Float_t* PedCorrCoeff0, 
	Float_t* PedCorrCoeff1)
{
  // Set coefficients for pedestal correlations
  if(PedCorrCoeff0 && PedCorrCoeff1) for(int t=0; t<44; t++){
    fPedCorrCoeff[0][t] = PedCorrCoeff0[t];
    fPedCorrCoeff[1][t] = PedCorrCoeff1[t];
  }
  else for(int t=0; t<44; t++){
    fPedCorrCoeff[0][t] = 0.;
    fPedCorrCoeff[1][t] = 0.;
  }
}

//________________________________________________________________
void AliZDCCalibData::SetEnCalib(Float_t* EnCalib) 
{
  if(EnCalib) for(int t=0; t<6; t++) fEnCalibration[t] = EnCalib[t];
  else for(int t=0; t<6; t++) fEnCalibration[t] = 0.;
}

