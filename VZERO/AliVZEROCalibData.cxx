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

ClassImp(AliVZEROCalibData)

//________________________________________________________________
AliVZEROCalibData::AliVZEROCalibData()
{
  // default constructor
  
    for(int t=0; t<64; t++) {
        fMeanHV[t]      = 100.0;
        fWidthHV[t]     = 0.0; 
	fTimeOffset[t]  = 0.0;
        fTimeGain[t]    = 1.0;
	fDeadChannel[t]= kFALSE;
    }
    for(int t=0; t<128; t++) {
        fPedestal[t]    = 0.0;     
        fSigma[t]       = 0.0;        
        fADCmean[t]     = 0.0;      
        fADCsigma[t]    = 0.0;
        fGain[t]        = 1.0;
    }
    for(int i=0; i<kNCIUBoards ;i++) {
	fTimeResolution[i]  = 25./256.;     // Default time resolution
	fWidthResolution[i] = 25./64.;     // Default time width resolution
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
   TString namst = "Calib_";
   namst += name;
   SetName(namst.Data());
   SetTitle(namst.Data());
   for(int t=0; t<64; t++) {
       fMeanHV[t]      = 100.0;
       fWidthHV[t]     = 0.0; 
       fTimeOffset[t]  = 0.0;
       fTimeGain[t]    = 1.0;
       fDeadChannel[t]= kFALSE;
    }
   for(int t=0; t<128; t++) {
       fPedestal[t]    = 0.0;     
       fSigma[t]       = 0.0;        
       fADCmean[t]     = 0.0;      
       fADCsigma[t]    = 0.0;
       fGain[t]        = 1.0;
   }
   for(int i=0; i<kNCIUBoards ;i++) {
       fTimeResolution[i]  = 25./256.;    // Default time resolution in ns / channel
       fWidthResolution[i] = 25./64.;     // Default time width resolution in ns / channel
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
      fADCsigma[t] = calibda.GetADCsigma(t);
      fGain[t]     = calibda.GetGain(t); }
      
  for(int t=0; t<64; t++) { 
      fMeanHV[t]       = calibda.GetMeanHV(t);
      fWidthHV[t]      = calibda.GetWidthHV(t);        
      fTimeOffset[t]   = calibda.GetTimeOffset(t);
      fTimeGain[t]     = calibda.GetTimeGain(t); 
      fDeadChannel[t]  = calibda.IsChannelDead(t);
  }  
  
  for(int i=0; i<kNCIUBoards ;i++) {
      fTimeResolution[i]  = calibda.GetTimeResolution(i);
      fWidthResolution[i] = calibda.GetWidthResolution(i);	  
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
      fADCsigma[t] = calibda.GetADCsigma(t);
      fGain[t]     = calibda.GetGain(t); }
      
  for(int t=0; t<64; t++) {
      fMeanHV[t]       = calibda.GetMeanHV(t);
      fWidthHV[t]      = calibda.GetWidthHV(t);        
      fTimeOffset[t]   = calibda.GetTimeOffset(t);
      fTimeGain[t]     = calibda.GetTimeGain(t); 
      fDeadChannel[t]  = calibda.IsChannelDead(t);
  }   
  for(int i=0; i<kNCIUBoards ;i++) {
      fTimeResolution[i]  = calibda.GetTimeResolution(i);
      fWidthResolution[i] = calibda.GetWidthResolution(i);	  
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
		Float_t val;
		if(aValue) {
			val = aValue->GetFloat();
			//AliInfo(Form("%s : %f",aliasName->String().Data(), val));
			SetParameter(aliasName->String(),val);
		}
	}	
	
	SetMeanHV(data->GetMeanHV());
	SetWidthHV(data->GetWidthHV());
  	SetDeadMap(data->GetDeadMap());

}
//_____________________________________________________________________________
void AliVZEROCalibData::SetParameter(TString name, Float_t val){
	// Set given parameter
	
	Int_t iBoard = -1;

	TSeqCollection* nameSplit = name.Tokenize("/");
	TObjString * boardName = (TObjString *)nameSplit->At(2);
	sscanf(boardName->String().Data(),"CIU%d",&iBoard);
		
	if(name.Contains("TimeResolution")) SetTimeResolution((UShort_t) val,iBoard-1);
	else if(name.Contains("WidthResolution")) SetWidthResolution((UShort_t) val,iBoard-1);
	else AliError(Form("No Setter found for FEE parameter : %s",name.Data()));
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
void AliVZEROCalibData::SetADCmean(Float_t* ADCmean) 
{
  if(ADCmean) for(int t=0; t<128; t++) fADCmean[t] = ADCmean[t];
  else for(int t=0; t<128; t++) fADCmean[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetADCsigma(Float_t* ADCsigma) 
{
  if(ADCsigma) for(int t=0; t<128; t++) fADCsigma[t] = ADCsigma[t];
  else for(int t=0; t<128; t++) fADCsigma[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetMeanHV(Float_t* MeanHV) 
{
  if(MeanHV) for(int t=0; t<64; t++) fMeanHV[t] = MeanHV[t];
  else for(int t=0; t<64; t++) fMeanHV[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetWidthHV(Float_t* WidthHV) 
{
  if(WidthHV) for(int t=0; t<64; t++) fWidthHV[t] = WidthHV[t];
  else for(int t=0; t<64; t++) fWidthHV[t] = 0.0;
}

//________________________________________________________________
void AliVZEROCalibData::SetDeadMap(Bool_t* deadMap) 
{
  if(deadMap) for(int t=0; t<64; t++) fDeadChannel[t] = deadMap[t];
  else for(int t=0; t<64; t++) fDeadChannel[t] = kFALSE;
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

//_____________________________________________________________________________
Float_t AliVZEROCalibData::GetMIPperADC(Int_t channel) const {
	
	// Computes the MIP conversion factor - MIP per ADC channel - 
	// Argument passed is the PM number (aliroot numbering)
	
	Float_t P0[64] = {
		7.094891, 7.124938, 7.089708, 7.098169, 7.094482, 7.147250, 7.170978, 7.183392, 
		7.145760, 7.148096, 7.153840, 7.143544, 7.186069, 7.194580, 7.203516, 7.195176, 
		7.188333, 7.198607, 7.209412, 7.226565, 7.221695, 7.205132, 7.191238, 7.227724, 
		7.232810, 7.252655, 7.230309, 7.273518, 7.273518, 7.242969, 7.252859, 7.252655, 
		7.026802, 7.079913, 7.134147, 7.092387, 7.079561, 7.072848, 7.123192, 7.003141, 
		7.024667, 7.124784, 7.123442, 7.129744, 7.110671, 7.143031, 7.139439, 7.178109, 
		7.247803, 7.139396, 7.293809, 7.094454, 6.992198, 7.206448, 7.244765, 7.056197, 
		7.263595, 7.138569, 7.089582, 7.215683, 7.266183, 7.165123, 7.243276, 7.235135 };
	Float_t P1[64] = {
		0.135569, 0.146405, 0.142425, 0.144278, 0.142307, 0.141648, 0.128477, 0.138239, 
		0.144173, 0.143419, 0.143572, 0.144482, 0.138024, 0.136542, 0.135955, 0.138537, 
		0.148521, 0.141999, 0.139627, 0.130014, 0.134970, 0.135635, 0.139094, 0.140634, 
		0.137971, 0.142080, 0.142793, 0.142778, 0.142778, 0.146045, 0.139133, 0.142080, 
		0.144121, 0.142311, 0.136564, 0.142686, 0.138792, 0.166285, 0.136387, 0.155391, 
		0.176082, 0.140408, 0.164738, 0.144270, 0.142766, 0.147486, 0.141951, 0.138012, 
		0.132394, 0.142849, 0.140477, 0.144592, 0.141558, 0.157646, 0.143758, 0.173385, 
		0.146489, 0.143279, 0.145230, 0.147203, 0.147333, 0.144979, 0.148597, 0.138985 };
	
	// High Voltage retrieval from Calibration Data Base:  
	Float_t  HV = fMeanHV[channel];  
	Float_t MIP = -1;
	if (HV>0)
	  MIP = 0.6/TMath::Exp((TMath::Log(HV) - P0[channel] )/P1[channel]);
	return MIP; 
	
}
//________________________________________________________________
void AliVZEROCalibData::SetTimeResolution(UShort_t *resols){
	// Set Time Resolution of the TDC
	if(resols)  for(int t=0; t<kNCIUBoards; t++) SetTimeResolution(resols[t],t);
	else AliFatal("Time Resolution not defined.");
	
}
//________________________________________________________________
void AliVZEROCalibData::SetTimeResolution(UShort_t resol, Int_t board)
{
	// Set Time Resolution of the TDC
	if((board<kNCIUBoards)) {
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
	} else AliError(Form("Board %d is not valid",board));
}
//________________________________________________________________
void AliVZEROCalibData::SetWidthResolution(UShort_t *resols){
	// Set Time Width Resolution of the TDC
	if(resols)  for(int t=0; t<kNCIUBoards; t++) SetWidthResolution(resols[t],t);
	else AliFatal("Width Resolution not defined.");
	
}
//________________________________________________________________
void AliVZEROCalibData::SetWidthResolution(UShort_t resol, Int_t board)
{
	// Set Time Width Resolution of the TDC
	if((board<kNCIUBoards)){
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
	}else AliError(Form("Board %d is not valid",board));
}

