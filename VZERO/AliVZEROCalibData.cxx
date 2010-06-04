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

#include "AliDCSValue.h"
#include "AliVZEROCalibData.h"
#include "AliVZERODataDCS.h"
#include "AliLog.h"

ClassImp(AliVZEROCalibData)

//________________________________________________________________
AliVZEROCalibData::AliVZEROCalibData()
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
AliVZEROCalibData::AliVZEROCalibData(const char* name)
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
  TNamed(calibda)
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
Float_t AliVZEROCalibData::GetMIPperADC(Int_t channel) const {
	
	// Computes the MIP conversion factor - MIP per ADC channel - 
	// Argument passed is the PM number (aliroot numbering)
	
	Float_t p0[64] = {
		7.094891, 7.124938, 7.089708, 7.098169, 7.094482, 7.147250, 7.170978, 7.183392, 
		7.145760, 7.148096, 7.153840, 7.143544, 7.186069, 7.194580, 7.203516, 7.195176, 
		7.188333, 7.198607, 7.209412, 7.226565, 7.221695, 7.205132, 7.191238, 7.227724, 
		7.232810, 7.252655, 7.230309, 7.140891, 7.273518, 7.242969, 7.252859, 7.252655, 
		7.026802, 7.079913, 7.134147, 7.092387, 7.079561, 7.072848, 7.123192, 7.003141, 
		7.024667, 7.124784, 7.123442, 7.129744, 7.110671, 7.143031, 7.139439, 7.178109, 
		7.247803, 7.139396, 7.293809, 7.094454, 6.992198, 7.206448, 7.244765, 7.056197, 
		7.263595, 7.138569, 7.089582, 7.215683, 7.266183, 7.165123, 7.243276, 7.235135 };
	Float_t p1[64] = {
		0.135569, 0.146405, 0.142425, 0.144278, 0.142307, 0.141648, 0.128477, 0.138239, 
		0.144173, 0.143419, 0.143572, 0.144482, 0.138024, 0.136542, 0.135955, 0.138537, 
		0.148521, 0.141999, 0.139627, 0.130014, 0.134970, 0.135635, 0.139094, 0.140634, 
		0.137971, 0.142080, 0.142793, 0.136054, 0.142778, 0.146045, 0.139133, 0.142080, 
		0.144121, 0.142311, 0.136564, 0.142686, 0.138792, 0.166285, 0.136387, 0.155391, 
		0.176082, 0.140408, 0.164738, 0.144270, 0.142766, 0.147486, 0.141951, 0.138012, 
		0.132394, 0.142849, 0.140477, 0.144592, 0.141558, 0.157646, 0.143758, 0.173385, 
		0.146489, 0.143279, 0.145230, 0.147203, 0.147333, 0.144979, 0.148597, 0.138985 };
	
	// High Voltage retrieval from Calibration Data Base:  
	Float_t  hv = fMeanHV[channel];  
	Float_t mip = -1;
	if (hv>0)
	  mip = 0.6/TMath::Exp((TMath::Log(hv) - p0[channel] )/p1[channel]);
	return mip; 
	
}

//________________________________________________________________
Float_t AliVZEROCalibData::GetGain(Int_t channel) const
{
  // Computes the PM gain factors
  // Argument passed is the PM number (aliroot numbering)
  Float_t a[64] = {-39.68,-35.83,-36.92,-36.42,-37.02,-37.50,-43.05,-39.39,
		   -36.62,-36.93,-37.30,-36.46,-39.51,-40.32,-39.92,-39.20,
		   -35.39,-37.95,-38.85,-42.76,-40.68,-40.32,-39.00,-37.36,
		   -39.64,-38.86,-37.59,-39.59,-37.97,-36.32,-38.88,-41.35,
		   -36.01,-36.82,-39.48,-36.86,-38.22,-32.55,-39.44,-35.08,
		   -29.91,-37.88,-33.25,-36.49,-37.25,-35.89,-40.31,-39.15,
		   -41.71,-37.07,-38.94,-36.04,-36.62,-32.96,-36.99,-30.71,
		   -36.66,-37.23,-35.98,-36.56,-35.64,-36.97,-35.88,-38.78};
  Float_t b[64] = {  7.40,  6.83,  7.02,  6.94,  7.03,  7.04,  7.79,  7.27,
		     6.92,  6.96,  7.01,  6.90,  7.28,  7.38,  7.33,  7.23,
		     6.71,  7.05,  7.17,  7.69,  7.41,  7.38,  7.21,  7.11,
		     7.26,  7.12,  6.98,  7.35,  6.99,  6.79,  7.13,  7.58,
		     6.95,  7.01,  7.33,  7.01,  7.21,  6.01,  7.34,  6.44,
		     5.68,  7.12,  6.07,  6.92,  7.04,  6.82,  7.04,  7.24,
		     7.53,  6.99,  7.10,  6.89,  7.07,  6.35,  6.88,  5.77,
		     6.81,  7.01,  6.89,  6.84,  6.68,  6.95,  6.73,  7.14};

  // High Voltage retrieval from Calibration Data Base:  
  Float_t hv = fMeanHV[channel];
  Float_t gain = 0;
  if (hv>0)
    gain = TMath::Exp(a[channel]+b[channel]*TMath::Log(hv));
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
