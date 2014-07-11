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

/* $Id: AliADCalibData.cxx,                                            */

#include <TMath.h>
#include <TObjString.h>
#include <TMap.h>
#include <TH1F.h>
#include <TH2F.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliADCalibData.h"
#include "AliADConst.h"
#include "AliLog.h"

ClassImp(AliADCalibData)


//________________________________________________________________
AliADCalibData::AliADCalibData():
  fLightYields(NULL),
  fPMGainsA(NULL),
  fPMGainsB(NULL)
{
  // default constructor
  
    for(int t=0; t<16; t++) {
        fMeanHV[t]      = 100.0;
        fWidthHV[t]     = 0.0; 
	fTimeOffset[t]  = 5.0;
        fTimeGain[t]    = 1.0;
	fDeadChannel[t]= kFALSE;
	fDiscriThr[t]  = 2.5;
    }
    for(int t=0; t<32; t++) {
        fPedestal[t]    = 0.0;     
        fSigma[t]       = 0.0;        
        fADCmean[t]     = 0.0;      
        fADCsigma[t]    = 0.0;
    }
    for(int i=0; i<kNCIUBoards ;i++) {
	fTimeResolution[i]  = 25./256.;     // Default time resolution
	fWidthResolution[i] = 25./64.;     // Default time width resolution
	fMatchWindow[i] = 4;
	fSearchWindow[i] = 16;
	fTriggerCountOffset[i] = 3247;
	fRollOver[i] = 3563;
    }

}
//________________________________________________________________
void AliADCalibData::Reset()
{
  
}

//________________________________________________________________
AliADCalibData::AliADCalibData(const char* name) :
  fLightYields(NULL),
  fPMGainsA(NULL),
  fPMGainsB(NULL)
{
  // Constructor
   TString namst = "Calib_";
   namst += name;
   SetName(namst.Data());
   SetTitle(namst.Data());
   for(int t=0; t<16; t++) {
       fMeanHV[t]      = 100.0;
       fWidthHV[t]     = 0.0; 
       fTimeOffset[t]  = 5.0;
       fTimeGain[t]    = 1.0;
       fDeadChannel[t]= kFALSE;
       fDiscriThr[t]  = 2.5;
    }
   for(int t=0; t<32; t++) {
       fPedestal[t]    = 0.0;     
       fSigma[t]       = 0.0;        
       fADCmean[t]     = 0.0;      
       fADCsigma[t]    = 0.0;
   }
   for(int i=0; i<kNCIUBoards ;i++) {
       fTimeResolution[i]  = 25./256.;    // Default time resolution in ns / channel
       fWidthResolution[i] = 25./64.;     // Default time width resolution in ns / channel
       fMatchWindow[i] = 4;
       fSearchWindow[i] = 16;
       fTriggerCountOffset[i] = 3247;
       fRollOver[i] = 3563;
   }
}

//________________________________________________________________
AliADCalibData::AliADCalibData(const AliADCalibData& calibda) :
  TNamed(calibda),
  fLightYields(NULL),
  fPMGainsA(NULL),
  fPMGainsB(NULL)
{
// copy constructor

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  for(int t=0; t<32; t++) { 
      fPedestal[t] = calibda.GetPedestal(t);
      fSigma[t]    = calibda.GetSigma(t);
      fADCmean[t]  = calibda.GetADCmean(t);
      fADCsigma[t] = calibda.GetADCsigma(t); }
      
  for(int t=0; t<16; t++) { 
      fMeanHV[t]       = calibda.GetMeanHV(t);
      fWidthHV[t]      = calibda.GetWidthHV(t);        
      fTimeOffset[t]   = calibda.GetTimeOffset(t);
      fTimeGain[t]     = calibda.GetTimeGain(t); 
      fDeadChannel[t]  = calibda.IsChannelDead(t);
      fDiscriThr[t]    = calibda.GetDiscriThr(t);
  }  
  
  for(int i=0; i<kNCIUBoards ;i++) {
      fTimeResolution[i]  = calibda.GetTimeResolution(i);
      fWidthResolution[i] = calibda.GetWidthResolution(i);	  
      fMatchWindow[i] = calibda.GetMatchWindow(i);
      fSearchWindow[i] = calibda.GetSearchWindow(i);
      fTriggerCountOffset[i] = calibda.GetTriggerCountOffset(i);
      fRollOver[i] = calibda.GetRollOver(i);
  }
  
}

//________________________________________________________________
AliADCalibData &AliADCalibData::operator =(const AliADCalibData& calibda)
{
// assignment operator

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  for(int t=0; t<32; t++) {
      fPedestal[t] = calibda.GetPedestal(t);
      fSigma[t]    = calibda.GetSigma(t);
      fADCmean[t]  = calibda.GetADCmean(t);
      fADCsigma[t] = calibda.GetADCsigma(t); }
      
  for(int t=0; t<16; t++) {
      fMeanHV[t]       = calibda.GetMeanHV(t);
      fWidthHV[t]      = calibda.GetWidthHV(t);        
      fTimeOffset[t]   = calibda.GetTimeOffset(t);
      fTimeGain[t]     = calibda.GetTimeGain(t); 
      fDeadChannel[t]  = calibda.IsChannelDead(t);
      fDiscriThr[t]    = calibda.GetDiscriThr(t);
  }   
  for(int i=0; i<kNCIUBoards ;i++) {
      fTimeResolution[i]  = calibda.GetTimeResolution(i);
      fWidthResolution[i] = calibda.GetWidthResolution(i);	  
      fMatchWindow[i] = calibda.GetMatchWindow(i);
      fSearchWindow[i] = calibda.GetSearchWindow(i);
      fTriggerCountOffset[i] = calibda.GetTriggerCountOffset(i);
      fRollOver[i] = calibda.GetRollOver(i);
  }
   
  return *this;
  
}

//________________________________________________________________
AliADCalibData::~AliADCalibData()
{
  // destructor
  if (fLightYields)
    delete [] fLightYields;
  if (fPMGainsA)
    delete [] fPMGainsA;
  if (fPMGainsB)
    delete [] fPMGainsB;
}

//________________________________________________________________
Int_t AliADCalibData::GetBoardNumber(Int_t channel)
{
  // Get FEE board number
  // from offline channel index
  if (channel >= 0 && channel < 16) return (channel);
  //if (channel >=8 && channel < 16) return (channel / 2);

  AliErrorClass(Form("Wrong channel index: %d",channel));
  return -1;
}

//________________________________________________________________
Float_t AliADCalibData::GetLightYields(Int_t channel)
{
  // Get the light yield efficiency
  // for a given channel
  if (!fLightYields) InitLightYields();

  if (channel >= 0 && channel < 64) {
    return fLightYields[channel];
  }

  AliError(Form("Wrong channel index: %d",channel));
  return 0;
}

//________________________________________________________________
void  AliADCalibData::InitLightYields()
{
  // Initialize the light yield factors
  // Read from a separate OCDB entry
  if (fLightYields) return;

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("VZERO/Calib/LightYields");
  if (!entry) AliFatal("AD light yields are not found in OCDB !");
  TH1F *yields = (TH1F*)entry->GetObject();

  fLightYields = new Float_t[16];
  for(Int_t i = 0 ; i < 16; ++i) {
    fLightYields[i] = yields->GetBinContent(i+1);
  }
}

//________________________________________________________________
void  AliADCalibData::InitPMGains()
{
  // Initialize the PM gain factors
  // Read from a separate OCDB entry
  if (fPMGainsA) return;

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("VZERO/Calib/PMGains");
  if (!entry) AliFatal("VZERO PM gains are not found in OCDB !");
  TH2F *gains = (TH2F*)entry->GetObject();

  fPMGainsA = new Float_t[16];
  fPMGainsB = new Float_t[16];
  for(Int_t i = 0 ; i < 16; ++i) {
    fPMGainsA[i] = gains->GetBinContent(i+1,1);
    fPMGainsB[i] = gains->GetBinContent(i+1,2);
  }
}

//________________________________________________________________
Float_t AliADCalibData::GetGain(Int_t channel)
{
  // Computes the PM gains
  // Argument passed is the PM number (aliroot numbering)
  if (!fPMGainsA) InitPMGains();

  // High Voltage retrieval from Calibration Data Base:  
  Float_t hv = fMeanHV[channel];
  Float_t gain = 0;
  if (hv>0)
    gain = TMath::Exp(fPMGainsA[channel]+fPMGainsB[channel]*TMath::Log(hv));
  return gain;
}

//________________________________________________________________
Float_t AliADCalibData::GetCalibDiscriThr(Int_t channel, Bool_t scaled)
{
  // The method returns actual TDC discri threshold
  // extracted from the data.
  //
  // In case scaled flag is set the threshold is scaled
  // so that to get 4.0 for a FEE threshold of 4.0.
  // In this way we avoid a change in the slewing correction
  // for the entire 2010 p-p data.
  //
  // The method is to be moved to OCDB object.

  Float_t thr = GetDiscriThr(channel);

  Float_t calThr = 0;
  if (thr <= 1.) 
    calThr = 3.1;
  else if (thr >= 2.)
    calThr = (3.1+1.15*thr-1.7);
  else
    calThr = (3.1-0.3*thr+0.3*thr*thr);

  if (scaled) calThr *= 4./(3.1+1.15*4.-1.7);

  return calThr;
}
