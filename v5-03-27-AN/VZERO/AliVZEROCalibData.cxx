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

#include <TMath.h>
#include <TObjString.h>
#include <TMap.h>
#include <TH1F.h>
#include <TH2F.h>

#include "AliDCSValue.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliVZEROCalibData.h"
#include "AliVZERODataDCS.h"
#include "AliVZEROConst.h"
#include "AliLog.h"

ClassImp(AliVZEROCalibData)

//________________________________________________________________
AliVZEROCalibData::AliVZEROCalibData():
  fLightYields(NULL),
  fPMGainsA(NULL),
  fPMGainsB(NULL)
{
  // default constructor
  
    for(int t=0; t<64; t++) {
        fMeanHV[t]      = 100.0;
        fWidthHV[t]     = 0.0; 
	fTimeOffset[t]  = 5.0;
        fTimeGain[t]    = 1.0;
	fDeadChannel[t]= kFALSE;
	fDiscriThr[t]  = 2.5;
    }
    for(int t=0; t<128; t++) {
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
void AliVZEROCalibData::Reset()
{
  // reset 
}

//________________________________________________________________
AliVZEROCalibData::AliVZEROCalibData(const char* name):
  fLightYields(NULL),
  fPMGainsA(NULL),
  fPMGainsB(NULL)
{
  // Constructor
   TString namst = "Calib_";
   namst += name;
   SetName(namst.Data());
   SetTitle(namst.Data());
   for(int t=0; t<64; t++) {
       fMeanHV[t]      = 100.0;
       fWidthHV[t]     = 0.0; 
       fTimeOffset[t]  = 5.0;
       fTimeGain[t]    = 1.0;
       fDeadChannel[t]= kFALSE;
       fDiscriThr[t]  = 2.5;
    }
   for(int t=0; t<128; t++) {
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
AliVZEROCalibData::AliVZEROCalibData(const AliVZEROCalibData& calibda) :
  TNamed(calibda),
  fLightYields(NULL),
  fPMGainsA(NULL),
  fPMGainsB(NULL)
{
// copy constructor

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  for(int t=0; t<128; t++) { 
      fPedestal[t] = calibda.GetPedestal(t);
      fSigma[t]    = calibda.GetSigma(t);
      fADCmean[t]  = calibda.GetADCmean(t);
      fADCsigma[t] = calibda.GetADCsigma(t); }
      
  for(int t=0; t<64; t++) { 
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
AliVZEROCalibData &AliVZEROCalibData::operator =(const AliVZEROCalibData& calibda)
{
// assignment operator

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  for(int t=0; t<128; t++) {
      fPedestal[t] = calibda.GetPedestal(t);
      fSigma[t]    = calibda.GetSigma(t);
      fADCmean[t]  = calibda.GetADCmean(t);
      fADCsigma[t] = calibda.GetADCsigma(t); }
      
  for(int t=0; t<64; t++) {
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
AliVZEROCalibData::~AliVZEROCalibData()
{
  // destructor
  if (fLightYields)
    delete [] fLightYields;
  if (fPMGainsA)
    delete [] fPMGainsA;
  if (fPMGainsB)
    delete [] fPMGainsB;
}
//_____________________________________________________________________________
void AliVZEROCalibData::FillDCSData(AliVZERODataDCS * data){
	// Set all parameters from the data get by the shuttle
	TMap * params = data->GetFEEParameters();
	TIter iter(params);	
	TObjString* aliasName;
	
	while ((  aliasName = (TObjString*) iter.Next() ))  {
		AliDCSValue* aValue = (AliDCSValue*) params->GetValue(aliasName);
		Int_t val;
		if(aValue) {
			val = aValue->GetUInt();
			AliInfo(Form("%s : %d",aliasName->String().Data(), val));
			SetParameter(aliasName->String(),val);
		}
	}	
	
	SetMeanHV(data->GetMeanHV());
	SetWidthHV(data->GetWidthHV());
  	SetDeadMap(data->GetDeadMap());

}
//_____________________________________________________________________________
void AliVZEROCalibData::SetParameter(TString name, Int_t val){
	// Set given parameter
	
	Int_t iBoard = -1;
	Int_t iChannel = -1;

	TSeqCollection* nameSplit = name.Tokenize("/");
	TObjString * boardName = (TObjString *)nameSplit->At(2);
	sscanf(boardName->String().Data(),"CIU%d",&iBoard);

	TString paramName = ((TObjString *)nameSplit->At(3))->String();
	Char_t channel[2] ; channel[1] = '\0';
	channel[0] = paramName[paramName.Sizeof()-2];
	sscanf(channel,"%d",&iChannel);
		
	if(name.Contains("TimeResolution")) SetTimeResolution((UShort_t) val,iBoard);
	else if(name.Contains("WidthResolution")) SetWidthResolution((UShort_t) val,iBoard);
	else if(name.Contains("MatchWindow")) SetMatchWindow((UInt_t) val,iBoard);
	else if(name.Contains("SearchWindow")) SetSearchWindow((UInt_t) val,iBoard);
	else if(name.Contains("TriggerCountOffset")) SetTriggerCountOffset((UInt_t) val,iBoard);
	else if(name.Contains("RollOver")) SetRollOver((UInt_t) val,iBoard);
	else if(name.Contains("DelayHit")) SetTimeOffset(0.01*(Float_t)val,iBoard,(iChannel-1));
	else if(name.Contains("DiscriThr")) SetDiscriThr(((Float_t)val-2040.)/112.,iBoard,(iChannel-1));
	else AliError(Form("No Setter found for FEE parameter : %s",name.Data()));
}

//________________________________________________________________
void AliVZEROCalibData::SetPedestal(const Float_t* Pedestal)
{
  if(Pedestal) for(int t=0; t<128; t++) fPedestal[t] = Pedestal[t];
  else for(int t=0; t<128; t++) fPedestal[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetSigma(const Float_t* Sigma)
{
  if(Sigma) for(int t=0; t<128; t++) fSigma[t] = Sigma[t];
  else for(int t=0; t<128; t++) fSigma[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetADCmean(const Float_t* ADCmean) 
{
  if(ADCmean) for(int t=0; t<128; t++) fADCmean[t] = ADCmean[t];
  else for(int t=0; t<128; t++) fADCmean[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetADCsigma(const Float_t* ADCsigma) 
{
  if(ADCsigma) for(int t=0; t<128; t++) fADCsigma[t] = ADCsigma[t];
  else for(int t=0; t<128; t++) fADCsigma[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetMeanHV(const Float_t* MeanHV) 
{
  if(MeanHV) for(int t=0; t<64; t++) fMeanHV[t] = MeanHV[t];
  else for(int t=0; t<64; t++) fMeanHV[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetWidthHV(const Float_t* WidthHV) 
{
  if(WidthHV) for(int t=0; t<64; t++) fWidthHV[t] = WidthHV[t];
  else for(int t=0; t<64; t++) fWidthHV[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetDeadMap(const Bool_t* deadMap) 
{
  if(deadMap) for(int t=0; t<64; t++) fDeadChannel[t] = deadMap[t];
  else for(int t=0; t<64; t++) fDeadChannel[t] = kFALSE;
}

//________________________________________________________________
void AliVZEROCalibData::SetTimeOffset(Float_t val, Int_t board, Int_t channel)
{
  Int_t ch = AliVZEROCalibData::GetOfflineChannelNumber(board,channel);
  if(ch >= 0){
    fTimeOffset[ch]=val;
    AliInfo(Form("Time offset for channel %d set to %f",ch,fTimeOffset[ch]));
  }
  else
    AliError("Board/Channel numbers are not valid");
}

//________________________________________________________________
void AliVZEROCalibData::SetTimeOffset(const Float_t* TimeOffset) 
{
  if(TimeOffset) for(int t=0; t<64; t++) fTimeOffset[t] = TimeOffset[t];
  else for(int t=0; t<64; t++) fTimeOffset[t] = 5.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetTimeGain(const Float_t* TimeGain) 
{
  if(TimeGain) for(int t=0; t<64; t++) fTimeGain[t] = TimeGain[t];
  else for(int t=0; t<64; t++) fTimeGain[t] = 0.0;
}

//_____________________________________________________________________________
Float_t AliVZEROCalibData::GetMIPperADC(Int_t channel) {
	
	// Computes the MIP conversion factor - MIP per ADC channel - 
	// Argument passed is the PM number (aliroot numbering)
	
  Float_t nPhPerMIP = (channel < 32) ? 6950 : 33690; 
  return 1./(nPhPerMIP*GetLightYields(channel)*0.18*TMath::Qe()*GetGain(channel)/kChargePerADC);
	
}

//_____________________________________________________________________________
Float_t AliVZEROCalibData::GetHV(Int_t channel, Float_t adcPerMip) {
	
  // Computes the HV value for certain ADC per MIP value
  // Arguments passed is the PM number (aliroot numbering) and
  // required value of ADC per MIP
  if (!fPMGainsA) InitPMGains();

  if (adcPerMip <= 0) return 0;
  Float_t nPhPerMIP = (channel < 32) ? 6950 : 33690;
  Float_t gain = adcPerMip/(nPhPerMIP*GetLightYields(channel)*0.18*TMath::Qe())*kChargePerADC;
  return TMath::Exp((TMath::Log(gain)-fPMGainsA[channel])/fPMGainsB[channel]);
}

//________________________________________________________________
Float_t AliVZEROCalibData::GetGain(Int_t channel)
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
void AliVZEROCalibData::SetTimeResolution(UShort_t *resols){
	// Set Time Resolution of the TDC
	if(resols)  for(int t=0; t<kNCIUBoards; t++) SetTimeResolution(resols[t],t);
	else AliError("Time Resolution not defined.");
	
}
//________________________________________________________________
void AliVZEROCalibData::SetTimeResolution(UShort_t resol, Int_t board)
{
	// Set Time Resolution of the TDC
	if((board>=0) && (board<kNCIUBoards)) {
		switch(resol){
			case 0:
				fTimeResolution[board] = 25./256.;
				break;
			case 1:
				fTimeResolution[board] = 25./128.;
				break;
			case 2:
				fTimeResolution[board] = 25./64.;
				break;
			case 3:
				fTimeResolution[board] = 25./32.;
				break;
			case 4:
				fTimeResolution[board] = 25./16.;
				break;
			case 5:
				fTimeResolution[board] = 25./8.;
				break;
			case 6:
				fTimeResolution[board] = 6.25;
				break;
			case 7:
				fTimeResolution[board] = 12.5;
				break;
		}
		AliInfo(Form("Time Resolution of board %d set to %f",board,fTimeResolution[board]));
	} else AliError(Form("Board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROCalibData::SetWidthResolution(UShort_t *resols){
	// Set Time Width Resolution of the TDC
	if(resols)  for(int t=0; t<kNCIUBoards; t++) SetWidthResolution(resols[t],t);
	else AliError("Width Resolution not defined.");
	
}
//________________________________________________________________
void AliVZEROCalibData::SetWidthResolution(UShort_t resol, Int_t board)
{
	// Set Time Width Resolution of the TDC
	if((board>=0) && (board<kNCIUBoards)){
		switch(resol){
			case 0:
				fWidthResolution[board] = 25./256.;
				break;
			case 1:
				fWidthResolution[board] = 25./128.;
				break;
			case 2:
				fWidthResolution[board] = 25./64.;
				break;
			case 3:
				fWidthResolution[board] = 25./32.;
				break;
			case 4:
				fWidthResolution[board] = 25./16.;
				break;
			case 5:
				fWidthResolution[board] = 25./8.;
				break;
			case 6:
				fWidthResolution[board] = 6.25;
				break;
			case 7:
				fWidthResolution[board] = 12.5;
				break;
			case 8:
				fWidthResolution[board] = 25.;
				break;
			case 9:
				fWidthResolution[board] = 50.;
				break;
			case 10:
				fWidthResolution[board] = 100.;
				break;
			case 11:
				fWidthResolution[board] = 200.;
				break;
			case 12:
				fWidthResolution[board] = 400.;
				break;
			case 13:
				fWidthResolution[board] = 800.;
				break;
				
		}
		AliInfo(Form("Width Resolution of board %d set to %f",board,fWidthResolution[board]));
	}else AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliVZEROCalibData::SetMatchWindow(UInt_t *windows)
{
  // Set Match window of the HPTDC
  // The units are 25ns
  if(windows)  for(Int_t b=0; b<kNCIUBoards; b++) SetMatchWindow(windows[b],b);
  else AliError("Match windows not defined.");
}

//________________________________________________________________
void AliVZEROCalibData::SetMatchWindow(UInt_t window, Int_t board)
{
  // Set Match window of the HPTDC
  // The units are 25ns
  if((board>=0) && (board<kNCIUBoards)){
    fMatchWindow[board] = window;
    AliInfo(Form("Match window of board %d set to %d",board,fMatchWindow[board]));
  }
  else
    AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliVZEROCalibData::SetSearchWindow(UInt_t *windows)
{
  // Set Search window of the HPTDC
  // The units are 25ns
  if(windows)  for(Int_t b=0; b<kNCIUBoards; b++) SetSearchWindow(windows[b],b);
  else AliError("Search windows not defined.");
}

//________________________________________________________________
void  AliVZEROCalibData::SetSearchWindow(UInt_t window, Int_t board)
{
  // Set Search window of the HPTDC
  // The units are 25ns
  if((board>=0) && (board<kNCIUBoards)){
    fSearchWindow[board] = window;
    AliInfo(Form("Search window of board %d set to %d",board,fSearchWindow[board]));
  }
  else
    AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliVZEROCalibData::SetTriggerCountOffset(UInt_t *offsets)
{
  // Set trigger-count offset of the HPTDC
  // The units are 25ns
  if(offsets)  for(Int_t b=0; b<kNCIUBoards; b++) SetTriggerCountOffset(offsets[b],b);
  else AliError("Trigger count offsets not defined.");
}

//________________________________________________________________
void AliVZEROCalibData::SetTriggerCountOffset(UInt_t offset, Int_t board)
{
  // Set trigger-count offsets of the HPTDC
  // The units are 25ns
  if((board>=0) && (board<kNCIUBoards)){
    fTriggerCountOffset[board] = offset;
    AliInfo(Form("Trigger-count offset of board %d set to %d",board,fTriggerCountOffset[board]));
  }
  else
    AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliVZEROCalibData::SetRollOver(UInt_t *offsets)
{
  // Set Roll-over of the HPTDC
  // The units are 25ns
  if(offsets)  for(Int_t b=0; b<kNCIUBoards; b++) SetRollOver(offsets[b],b);
  else AliError("Roll-over offsets not defined.");
}

//________________________________________________________________
void AliVZEROCalibData::SetRollOver(UInt_t offset, Int_t board)
{
  // Set Roll-over of the HPTDC
  // The units are 25ns
  if((board>=0) && (board<kNCIUBoards)){
    fRollOver[board] = offset;
    AliInfo(Form("Roll-over offset of board %d set to %d",board,fRollOver[board]));
  }
  else
    AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliVZEROCalibData::SetDiscriThr(Float_t thr, Int_t board, Int_t channel)
{
  // Set the TDC discriminator
  // threshold values expressed in units of ADC
  Int_t ch = AliVZEROCalibData::GetOfflineChannelNumber(board,channel);
  if(ch >= 0){
    if (thr > 0) {
      fDiscriThr[ch]=thr;
      AliInfo(Form("Discriminator threshold for channel %d set to %f",ch,fDiscriThr[ch]));
    }
    else {
      AliWarning(Form("Ignore wrong threshold value (%f) for channel %d !",thr,ch));
    }
  }
  else
    AliError("Board/Channel numbers are not valid");
}

//________________________________________________________________
void AliVZEROCalibData::SetDiscriThr(const Float_t* thresholds) 
{
  // Set the TDC discriminator
  // threshold values expressed in units of ADC
  if(thresholds) for(int t=0; t<64; t++) fDiscriThr[t] = thresholds[t];
  else for(int t=0; t<64; t++) fDiscriThr[t] = 2.5;
}

Int_t AliVZEROCalibData::GetOfflineChannelNumber(Int_t board, Int_t channel)
{
  // Get the offline channel number from
  // the FEE board and channel indexes

  if (board < 0 || board >= 8) {
    AliErrorClass(Form("Wrong FEE board number: %d",board));
    return -1;
  }
  if (channel < 0 || channel >= 8) {
    AliErrorClass(Form("Wrong FEE channel number: %d",channel));
    return -1;
  }

  Int_t offCh = (board < 4) ? (8 * board + 32) : (8 * board -32);
  offCh += (7 - channel);

  return offCh;
}

Int_t AliVZEROCalibData::GetBoardNumber(Int_t channel)
{
  // Get FEE board number
  // from offline channel index
  if (channel >= 0 && channel < 32) return (channel / 8 + 4);
  if (channel >=32 && channel < 64) return (channel / 8 - 4);

  AliErrorClass(Form("Wrong channel index: %d",channel));
  return -1;
}

Int_t AliVZEROCalibData::GetFEEChannelNumber(Int_t channel)
{
  // Get FEE channel number
  // from offline channel index
  if (channel >= 0 && channel < 64) return (7 - (channel % 8));

  AliErrorClass(Form("Wrong channel index: %d",channel));
  return -1;
}

Float_t AliVZEROCalibData::GetLightYields(Int_t channel)
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

void  AliVZEROCalibData::InitLightYields()
{
  // Initialize the light yield factors
  // Read from a separate OCDB entry
  if (fLightYields) return;

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("VZERO/Calib/LightYields");
  if (!entry) AliFatal("VZERO light yields are not found in OCDB !");
  TH1F *yields = (TH1F*)entry->GetObject();

  fLightYields = new Float_t[64];
  for(Int_t i = 0 ; i < 64; ++i) {
    fLightYields[i] = yields->GetBinContent(i+1);
  }
}

void  AliVZEROCalibData::InitPMGains()
{
  // Initialize the PM gain factors
  // Read from a separate OCDB entry
  if (fPMGainsA) return;

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("VZERO/Calib/PMGains");
  if (!entry) AliFatal("VZERO PM gains are not found in OCDB !");
  TH2F *gains = (TH2F*)entry->GetObject();

  fPMGainsA = new Float_t[64];
  fPMGainsB = new Float_t[64];
  for(Int_t i = 0 ; i < 64; ++i) {
    fPMGainsA[i] = gains->GetBinContent(i+1,1);
    fPMGainsB[i] = gains->GetBinContent(i+1,2);
  }
}

Float_t AliVZEROCalibData::GetCalibDiscriThr(Int_t channel, Bool_t scaled)
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
