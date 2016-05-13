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
#include <Riostream.h>
#include <bitset>

#include <TMath.h>
#include <TObjString.h>
#include <TMap.h>
#include <TH1F.h>
#include <TH2F.h>

#include "AliDCSValue.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliADDataDCS.h"
#include "AliADCalibData.h"
#include "AliADConst.h"
#include "AliADLogicalSignal.h"
#include "AliLog.h"

ClassImp(AliADCalibData);

//________________________________________________________________
AliADCalibData::AliADCalibData():
  fLightYields(NULL),
  fPMGainsA(NULL),
  fPMGainsB(NULL),
  fThrCalibA(NULL),
  fThrCalibB(NULL),
  fBBAThreshold(0),
  fBBCThreshold(0) ,  
  fBGAThreshold(0) ,  
  fBGCThreshold(0) ,  
  fBBAForBGThreshold(0) ,  
  fBBCForBGThreshold(0) ,   
  fMultADAThrLow(0) ,  
  fMultADAThrHigh(0) , 
  fMultADCThrLow(0) ,  
  fMultADCThrHigh(0)
{
  // default constructor
  
  for(int t=0; t<16; t++) {
    fMeanHV[t]      = 1400.0;
    fWidthHV[t]     = 0.0; 
    fTimeOffset[t]  = 250.0;
    fTimeGain[t]    = 1.0;
    fDeadChannel[t]= kFALSE;
    fDiscriThr[t]  = 2.5;
  }
  for(int t=0; t<32; t++) {
    fPedestal[t]    = 3.0;     
    fSigma[t]       = 1.0;        
    fADCmean[t]     = 0.0;      
    fADCsigma[t]    = 0.0;
  }
  for(int i=0; i<kNCIUBoards ;i++) {
    fTimeResolution[i]  = 25./256.;     // Default time resolution
    fWidthResolution[i] = 25./64.;     // Default time width resolution
    fMatchWindow[i] = 16;
    fSearchWindow[i] = 16;
    fTriggerCountOffset[i] = 3553;
    fRollOver[i] = 3563;
  }
  for(int i=0; i<kNCIUBoards ;i++) {
    fClk1Win1[i] = fClk1Win2[i] = 0;
    fDelayClk1Win1[i] = fDelayClk1Win2[i] = 0;
    fClk2Win1[i] = fClk2Win2[i] = 0;
    fDelayClk2Win1[i] = fDelayClk2Win2[i] = 0;
    fLatchWin1[i] = fLatchWin2[i] = 0;
    fResetWin1[i] = fResetWin2[i] = 0;
    fPedestalSubtraction[i] = kFALSE;
  }
  for(Int_t j = 0; j < 16; ++j) {
    fEnableCharge[j] = kFALSE;
    fEnableTiming[j] = kTRUE;
    fPedestalOdd[j] = fPedestalEven[j] = 0;
    fPedestalCutOdd[j] = fPedestalCutEven[j] = 0;
  }
  for(Int_t i = 0; i < 5; ++i) fTriggerSelected[i] = 0;

}
//________________________________________________________________
void AliADCalibData::Reset()
{
  
}

//________________________________________________________________
AliADCalibData::AliADCalibData(const char* name) :
  fLightYields(NULL),
  fPMGainsA(NULL),
  fPMGainsB(NULL),
  fThrCalibA(NULL),
  fThrCalibB(NULL),
  fBBAThreshold(0),
  fBBCThreshold(0) ,  
  fBGAThreshold(0) ,  
  fBGCThreshold(0) ,  
  fBBAForBGThreshold(0) ,  
  fBBCForBGThreshold(0) ,   
  fMultADAThrLow(0) ,  
  fMultADAThrHigh(0) , 
  fMultADCThrLow(0) ,  
  fMultADCThrHigh(0)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  for(int t=0; t<16; t++) {
    fMeanHV[t]      = 1500.0;
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
  for(int i=0; i<kNCIUBoards ;i++) {
    fClk1Win1[i] = fClk1Win2[i] = 0;
    fDelayClk1Win1[i] = fDelayClk1Win2[i] = 0;
    fClk2Win1[i] = fClk2Win2[i] = 0;
    fDelayClk2Win1[i] = fDelayClk2Win2[i] = 0;
    fLatchWin1[i] = fLatchWin2[i] = 0;
    fResetWin1[i] = fResetWin2[i] = 0;
    fPedestalSubtraction[i] = kFALSE;
  }
  for(Int_t j = 0; j < 16; ++j) {
    fEnableCharge[j] = kFALSE;
    fEnableTiming[j] = kTRUE;
    fPedestalOdd[j] = fPedestalEven[j] = 0;
    fPedestalCutOdd[j] = fPedestalCutEven[j] = 0;
  }
  for(Int_t i = 0; i < 5; ++i) fTriggerSelected[i] = 0;
}

//________________________________________________________________
AliADCalibData::AliADCalibData(const AliADCalibData& calibda) :
  TNamed(calibda),
  fLightYields(NULL),
  fPMGainsA(NULL),
  fPMGainsB(NULL),
  fThrCalibA(NULL),
  fThrCalibB(NULL)
{
  // copy constructor

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  fBBAThreshold = calibda.GetBBAThreshold();
  fBBCThreshold = calibda.GetBBCThreshold();  
  fBGAThreshold = calibda.GetBGAThreshold();  
  fBGCThreshold = calibda.GetBGCThreshold();  
  fBBAForBGThreshold = calibda.GetBBAForBGThreshold();  
  fBBCForBGThreshold = calibda.GetBBCForBGThreshold();   
  fMultADAThrLow = calibda.GetMultADAThrLow();  
  fMultADAThrHigh = calibda.GetMultADAThrHigh();
  fMultADCThrLow = calibda.GetMultADCThrLow();  
  fMultADCThrHigh = calibda.GetMultADCThrHigh();
  
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
  
  for(int i=0; i<kNCIUBoards ;i++) {
    fClk1Win1[i] = calibda.GetClk1Win1(i);
    fClk1Win2[i] = calibda.GetClk1Win2(i);
    fDelayClk1Win1[i] = calibda.GetDelayClk1Win1(i);
    fDelayClk1Win2[i] = calibda.GetDelayClk1Win2(i);
    fClk2Win1[i] = calibda.GetClk2Win1(i);
    fClk2Win2[i] = calibda.GetClk2Win2(i);
    fDelayClk2Win1[i] = calibda.GetDelayClk2Win1(i);
    fDelayClk2Win2[i] = calibda.GetDelayClk2Win2(i);
    fLatchWin1[i] = calibda.GetLatchWin1(i);
    fLatchWin2[i] = calibda.GetLatchWin2(i);
    fResetWin1[i] = calibda.GetResetWin1(i);
    fResetWin2[i] = calibda.GetResetWin2(i);
    fPedestalSubtraction[i] = GetPedestalSubtraction(i);
  }
  for(Int_t j = 0; j < 16; ++j) {
    fEnableCharge[j] = calibda.GetEnableCharge(j);
    fEnableTiming[j] = calibda.GetEnableTiming(j);
    fPedestalOdd[j] = calibda.GetOnlinePedestal(0,j);
    fPedestalEven[j] = calibda.GetOnlinePedestal(1,j);;
    fPedestalCutOdd[j] = calibda.GetOnlinePedestalCut(0,j);
    fPedestalCutEven[j] = calibda.GetOnlinePedestal(1,j);
  }
  for(Int_t i = 0; i < 5; ++i) fTriggerSelected[i] = calibda.GetTriggerSelected(i);
  
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
  
  for(int i=0; i<kNCIUBoards ;i++) {
    fClk1Win1[i] = calibda.GetClk1Win1(i);
    fClk1Win2[i] = calibda.GetClk1Win2(i);
    fDelayClk1Win1[i] = calibda.GetDelayClk1Win1(i);
    fDelayClk1Win2[i] = calibda.GetDelayClk1Win2(i);
    fClk2Win1[i] = calibda.GetClk2Win1(i);
    fClk2Win2[i] = calibda.GetClk2Win2(i);
    fDelayClk2Win1[i] = calibda.GetDelayClk2Win1(i);
    fDelayClk2Win2[i] = calibda.GetDelayClk2Win2(i);
    fLatchWin1[i] = calibda.GetLatchWin1(i);
    fLatchWin2[i] = calibda.GetLatchWin2(i);
    fResetWin1[i] = calibda.GetResetWin1(i);
    fResetWin2[i] = calibda.GetResetWin2(i);
    fPedestalSubtraction[i] = GetPedestalSubtraction(i);
  }
  for(Int_t j = 0; j < 16; ++j) {
    fEnableCharge[j] = calibda.GetEnableCharge(j);
    fEnableTiming[j] = calibda.GetEnableTiming(j);
    fPedestalOdd[j] = calibda.GetOnlinePedestal(0,j);
    fPedestalEven[j] = calibda.GetOnlinePedestal(1,j);;
    fPedestalCutOdd[j] = calibda.GetOnlinePedestalCut(0,j);
    fPedestalCutEven[j] = calibda.GetOnlinePedestal(1,j);
  }
  for(Int_t i = 0; i < 5; ++i) fTriggerSelected[i] = calibda.GetTriggerSelected(i); 
   
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
  if (fThrCalibA)
    delete [] fThrCalibA;
  if (fThrCalibB)
    delete [] fThrCalibB;
}

//________________________________________________________________
Float_t AliADCalibData::GetLightYields(Int_t channel)
{
  // Get the light yield efficiency
  // for a given channel
  if (!fLightYields) InitLightYields();

  if (channel >= 0 && channel < 16) {
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

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("AD/Calib/LightYields");
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

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("AD/Calib/PMGains");
  if (!entry) AliFatal("AD PM gains are not found in OCDB !");
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
    //gain = TMath::Exp(fPMGainsA[channel]+fPMGainsB[channel]*TMath::Log(hv));
    gain = TMath::Power(hv/fPMGainsA[channel],fPMGainsB[channel])*kADChargePerADC/kNPhotonsPerMIP;
  return gain;
}
//________________________________________________________________
Float_t AliADCalibData::GetADCperMIP(Int_t channel)
{
  // Computes the MIP position wrt gain curve and HV
  // Argument passed is the PM number (aliroot numbering)
  if (!fPMGainsA) InitPMGains();

  // High Voltage retrieval from Calibration Data Base:  
  Float_t hv = fMeanHV[channel];
  Float_t ADCperMIP = 0;
  if (hv>0)
    ADCperMIP = TMath::Power(hv/fPMGainsA[channel],fPMGainsB[channel]);
  return ADCperMIP;
}
//________________________________________________________________
void  AliADCalibData::InitCalibThresholds()
{
  // Initialize the PM gain factors
  // Read from a separate OCDB entry
  if (fThrCalibA) return;

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("AD/Calib/Thresholds");
  if (!entry) AliFatal("AD Thresholds are not found in OCDB !");
  TH2F *thrCalib = (TH2F*)entry->GetObject();

  fThrCalibA = new Float_t[16];
  fThrCalibB = new Float_t[16];
  for(Int_t i = 0 ; i < 16; ++i) {
    fThrCalibA[i] = thrCalib->GetBinContent(i+1,1);
    fThrCalibB[i] = thrCalib->GetBinContent(i+1,2);
  }
}

//________________________________________________________________
Float_t AliADCalibData::GetCalibDiscriThr(Int_t channel)
{
  // The method returns actual TDC discri threshold
  // extracted from the data.
  if (!fThrCalibA) InitCalibThresholds();
  
  Float_t thr = GetDiscriThr(channel);

  Float_t calThr = fThrCalibA[channel]*thr + fThrCalibB[channel];

  return calThr;
}
//_____________________________________________________________________________
void AliADCalibData::FillDCSData(AliADDataDCS * data){
  // Set all parameters from the data get by the shuttle
  TMap * params = data->GetFEEParameters();
  TIter iter(params);	
  TObjString* aliasName;
	
  while ((  aliasName = (TObjString*) iter.Next() ))  {
    AliDCSValue* aValue = (AliDCSValue*) params->GetValue(aliasName);
    UInt_t val;
    if(aValue) {
      val = aValue->GetUInt();
      SetParameter(aliasName->String(),val);
    }
  }	
	
  SetMeanHV(data->GetMeanHV());
  SetWidthHV(data->GetWidthHV());
  SetDeadMap(data->GetDeadMap());

}
//_____________________________________________________________________________
void AliADCalibData::SetADCperMIP(Int_t nADCperMIP){
  //Sets HV in a way to have uniform gains on PM according  
  //to number of ADC per MIP 
  if (!fPMGainsA) InitPMGains();
  Double_t hv = 0;
  for(Int_t channel = 0; channel<16; channel++){
    hv = TMath::Power(nADCperMIP,1/fPMGainsB[channel])*fPMGainsA[channel];
    SetMeanHV(hv,channel);
    AliInfo(Form("HV on channel %d set to %f V",channel,hv));
  }
}

//_____________________________________________________________________________
void AliADCalibData::SetParameter(TString name, Int_t val){
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
		
  else if(name.Contains("DelayClk1Win1")) SetDelayClk1Win1((UShort_t) val,iBoard);
  else if(name.Contains("Clk1Win1")) SetClk1Win1((UShort_t) val,iBoard);
  else if(name.Contains("DelayClk1Win2")) SetDelayClk1Win2((UShort_t) val,iBoard);
  else if(name.Contains("Clk1Win2")) SetClk1Win2((UShort_t) val,iBoard);
  else if(name.Contains("DelayClk2Win1")) SetDelayClk2Win1((UShort_t) val,iBoard);
  else if(name.Contains("Clk2Win1")) SetClk2Win1((UShort_t) val,iBoard);
  else if(name.Contains("DelayClk2Win2")) SetDelayClk2Win2((UShort_t) val,iBoard);
  else if(name.Contains("Clk2Win2")) SetClk2Win2((UShort_t) val,iBoard);
  else if(name.Contains("LatchWin1")) SetLatchWin1((UShort_t) val,iBoard);
  else if(name.Contains("LatchWin2")) SetLatchWin2((UShort_t) val,iBoard);
  else if(name.Contains("ResetWin1")) SetResetWin1((UShort_t) val,iBoard);
  else if(name.Contains("ResetWin2")) SetResetWin2((UShort_t) val,iBoard);
  else if(name.Contains("PedestalSubtraction")) SetPedestalSubtraction((Bool_t) val,iBoard);
  else if(name.Contains("BBAThreshold")) SetBBAThreshold((UShort_t) val);
  else if(name.Contains("BBCThreshold")) SetBBCThreshold((UShort_t) val);
  else if(name.Contains("BGAThreshold")) SetBGAThreshold((UShort_t) val);
  else if(name.Contains("BGCThreshold")) SetBGCThreshold((UShort_t) val);
  else if(name.Contains("BBAForBGThreshold")) SetBBAForBGThreshold((UShort_t) val);
  else if(name.Contains("BBCForBGThreshold")) SetBBCForBGThreshold((UShort_t) val);
  else if(name.Contains("MultADAThrLow")) SetMultADAThrLow((UShort_t) val);
  else if(name.Contains("MultADAThrHigh")) SetMultADAThrHigh((UShort_t) val);
  else if(name.Contains("MultADCThrLow")) SetMultADCThrLow((UShort_t) val);
  else if(name.Contains("MultADCThrHigh")) SetMultADCThrHigh((UShort_t) val);
  else if(name.Contains("TriggerSelect")) SetTriggerSelected((UShort_t) val, iChannel-1 );
  else if(name.Contains("EnableCharge")) SetEnableCharge((Bool_t) val, iBoard , iChannel-1);
  else if(name.Contains("EnableTiming")) SetEnableTiming((Bool_t) val, iBoard , iChannel-1);
  else if(name.Contains("PedOdd")) SetOnlinePedestal((UShort_t) val, 1, iBoard, iChannel-1);
  else if(name.Contains("PedEven")) SetOnlinePedestal((UShort_t) val, 0, iBoard, iChannel-1);
  else if(name.Contains("PedCutOdd")) SetOnlinePedestalCut((UShort_t) val, 1, iBoard, iChannel-1);
  else if(name.Contains("PedCutEven")) SetOnlinePedestalCut((UShort_t) val, 0, iBoard, iChannel-1);
	
  else AliError(Form("No Setter found for FEE parameter : %s",name.Data()));
  //
  delete nameSplit;
}

//________________________________________________________________
void AliADCalibData::SetPedestal(const Float_t* Pedestal)
{
  if(Pedestal) for(int t=0; t<32; t++) fPedestal[t] = Pedestal[t];
  else for(int t=0; t<32; t++) fPedestal[t] = 0.0;
}

//________________________________________________________________
void AliADCalibData::SetSigma(const Float_t* Sigma)
{
  if(Sigma) for(int t=0; t<32; t++) fSigma[t] = Sigma[t];
  else for(int t=0; t<32; t++) fSigma[t] = 0.0;
}

//________________________________________________________________
void AliADCalibData::SetADCmean(const Float_t* ADCmean) 
{
  if(ADCmean) for(int t=0; t<32; t++) fADCmean[t] = ADCmean[t];
  else for(int t=0; t<32; t++) fADCmean[t] = 0.0;
}

//________________________________________________________________
void AliADCalibData::SetADCsigma(const Float_t* ADCsigma) 
{
  if(ADCsigma) for(int t=0; t<32; t++) fADCsigma[t] = ADCsigma[t];
  else for(int t=0; t<32; t++) fADCsigma[t] = 0.0;
}

//________________________________________________________________
void AliADCalibData::SetMeanHV(const Float_t* MeanHV) 
{
  if(MeanHV) for(int t=0; t<16; t++) fMeanHV[t] = MeanHV[t];
  else for(int t=0; t<16; t++) fMeanHV[t] = 0.0;
}

//________________________________________________________________
void AliADCalibData::SetWidthHV(const Float_t* WidthHV) 
{
  if(WidthHV) for(int t=0; t<16; t++) fWidthHV[t] = WidthHV[t];
  else for(int t=0; t<16; t++) fWidthHV[t] = 0.0;
}

//________________________________________________________________
void AliADCalibData::SetDeadMap(const Bool_t* deadMap) 
{
  if(deadMap) for(int t=0; t<16; t++) fDeadChannel[t] = deadMap[t];
  else for(int t=0; t<16; t++) fDeadChannel[t] = kFALSE;
}

//________________________________________________________________
void AliADCalibData::SetTimeOffset(Float_t val, Int_t board, Int_t channel)
{
  Int_t ch = AliADCalibData::GetOfflineChannelNumber(board,channel);
  if(ch >= 0) fTimeOffset[ch]=val;
  else
    AliError("Board/Channel numbers are not valid");
}

//________________________________________________________________
void AliADCalibData::SetTimeOffset(const Float_t* TimeOffset) 
{
  if(TimeOffset) for(int t=0; t<16; t++) fTimeOffset[t] = TimeOffset[t];
  else for(int t=0; t<16; t++) fTimeOffset[t] = 5.0;
}
//________________________________________________________________
void AliADCalibData::SetTimeGain(const Float_t* TimeGain) 
{
  if(TimeGain) for(int t=0; t<16; t++) fTimeGain[t] = TimeGain[t];
  else for(int t=0; t<16; t++) fTimeGain[t] = 0.0;
}
//________________________________________________________________
void AliADCalibData::SetTimeResolution(UShort_t *resols){
  // Set Time Resolution of the TDC
  if(resols)  for(int t=0; t<kNCIUBoards; t++) SetTimeResolution(resols[t],t);
  else AliError("Time Resolution not defined.");
	
}
//________________________________________________________________
void AliADCalibData::SetTimeResolution(UShort_t resol, Int_t board)
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
  } else AliError(Form("Board %d is not valid",board));
}
//________________________________________________________________
void AliADCalibData::SetWidthResolution(UShort_t *resols){
  // Set Time Width Resolution of the TDC
  if(resols)  for(int t=0; t<kNCIUBoards; t++) SetWidthResolution(resols[t],t);
  else AliError("Width Resolution not defined.");
	
}
//________________________________________________________________
void AliADCalibData::SetWidthResolution(UShort_t resol, Int_t board)
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
  }else AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliADCalibData::SetMatchWindow(UInt_t *windows)
{
  // Set Match window of the HPTDC
  // The units are 25ns
  if(windows)  for(Int_t b=0; b<kNCIUBoards; b++) SetMatchWindow(windows[b],b);
  else AliError("Match windows not defined.");
}

//________________________________________________________________
void AliADCalibData::SetMatchWindow(UInt_t window, Int_t board)
{
  // Set Match window of the HPTDC
  // The units are 25ns
  if((board>=0) && (board<kNCIUBoards)) fMatchWindow[board] = window;
  else
    AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliADCalibData::SetSearchWindow(UInt_t *windows)
{
  // Set Search window of the HPTDC
  // The units are 25ns
  if(windows)  for(Int_t b=0; b<kNCIUBoards; b++) SetSearchWindow(windows[b],b);
  else AliError("Search windows not defined.");
}

//________________________________________________________________
void  AliADCalibData::SetSearchWindow(UInt_t window, Int_t board)
{
  // Set Search window of the HPTDC
  // The units are 25ns
  if((board>=0) && (board<kNCIUBoards)) fSearchWindow[board] = window;
  else
    AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliADCalibData::SetTriggerCountOffset(UInt_t *offsets)
{
  // Set trigger-count offset of the HPTDC
  // The units are 25ns
  if(offsets)  for(Int_t b=0; b<kNCIUBoards; b++) SetTriggerCountOffset(offsets[b],b);
  else AliError("Trigger count offsets not defined.");
}

//________________________________________________________________
void AliADCalibData::SetTriggerCountOffset(UInt_t offset, Int_t board)
{
  // Set trigger-count offsets of the HPTDC
  // The units are 25ns
  if((board>=0) && (board<kNCIUBoards)) fTriggerCountOffset[board] = offset;
  else
    AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliADCalibData::SetRollOver(UInt_t *offsets)
{
  // Set Roll-over of the HPTDC
  // The units are 25ns
  if(offsets)  for(Int_t b=0; b<kNCIUBoards; b++) SetRollOver(offsets[b],b);
  else AliError("Roll-over offsets not defined.");
}

//________________________________________________________________
void AliADCalibData::SetRollOver(UInt_t offset, Int_t board)
{
  // Set Roll-over of the HPTDC
  // The units are 25ns
  if((board>=0) && (board<kNCIUBoards)) fRollOver[board] = offset;
  else
    AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
void AliADCalibData::SetDiscriThr(Float_t thr, Int_t board, Int_t channel)
{
  // Set the TDC discriminator
  // threshold values expressed in units of ADC
  Int_t ch = AliADCalibData::GetOfflineChannelNumber(board,channel);
  if(ch >= 0) fDiscriThr[ch]=thr;
  else
    AliError("Board/Channel numbers are not valid");
}

//________________________________________________________________
void AliADCalibData::SetDiscriThr(const Float_t* thresholds) 
{
  // Set the TDC discriminator
  // threshold values expressed in units of ADC
  if(thresholds) for(int t=0; t<16; t++) fDiscriThr[t] = thresholds[t];
  else for(int t=0; t<16; t++) fDiscriThr[t] = 2.5;
}


//________________________________________________________________
void AliADCalibData::SetOnlinePedestalCut(UShort_t val,Int_t integrator, Int_t channel)
{
  // Set Pedestal Cut of individual channel 
  if(channel>=0 && channel<16) {
    if(integrator) fPedestalCutOdd[channel] = val;
    else fPedestalCutEven[channel] = val;
  } else AliError(Form("Impossible to write at : Channel %d",channel));
}
//________________________________________________________________
void AliADCalibData::SetOnlinePedestalCut(UShort_t val,Int_t integrator, Int_t board, Int_t channel)
{
  Int_t ch = AliADCalibData::GetOfflineChannelNumber(board,channel);
  if(ch>=0 && ch<16) {
    if(integrator) fPedestalCutOdd[ch] = val;
    else fPedestalCutEven[ch] = val;
  }
  else
    AliError("Board/Channel numbers are not valid");
}
//________________________________________________________________
UShort_t AliADCalibData::GetOnlinePedestalCut(Int_t integrator, Int_t channel) const
{
  // Get Pedestal Cut of individual channel 
  if(channel>=0 && channel<16) {
    if(integrator) return(fPedestalCutOdd[channel]);
    else return(fPedestalCutEven[channel]);
  }else AliError(Form("Impossible to read at : Channel %d",channel));
  return 0;
}
//________________________________________________________________
void AliADCalibData::SetOnlinePedestal(UShort_t val, Int_t integrator, Int_t channel)
{
  // Set Pedestal of individual channel 
  if(channel>=0 && channel<16) {
    if(integrator) fPedestalOdd[channel] = val;
    else fPedestalEven[channel] = val;
  } else AliError(Form("Impossible to write at : Channel %d ; Integrator %d ",channel,integrator));
}
//________________________________________________________________
void AliADCalibData::SetOnlinePedestal(UShort_t val,Int_t integrator, Int_t board, Int_t channel)
{
  Int_t ch = AliADCalibData::GetOfflineChannelNumber(board,channel);
  if(ch>=0 && ch<16) {
    if(integrator) fPedestalOdd[ch] = val;
    else fPedestalEven[ch] = val;
  }
  else
    AliError("Board/Channel numbers are not valid");
}
//________________________________________________________________
UShort_t AliADCalibData::GetOnlinePedestal(Int_t integrator, Int_t channel) const
{
  // Get Pedestal of individual channel 
  if(channel>=0 && channel<16) {
    if(integrator) return(fPedestalOdd[channel]);
    else return(fPedestalEven[channel]);
  } else AliError(Form("Impossible to read at : Channel %d",channel));
  return 0;
}
//________________________________________________________________
void AliADCalibData::SetEnableCharge(Bool_t val, Int_t channel)
{
  // Set the channels enabled for Charge triggers
  if(channel>=0 && channel<16) fEnableCharge[channel] = val;
  else AliError(Form("Impossible to write at : Channel %d",channel));
}
//________________________________________________________________
Bool_t AliADCalibData::GetEnableCharge(Int_t channel) const
{
  // Get the channels enabled for Charge triggers
  if(channel>=0 && channel<16) return(fEnableCharge[channel]);
  else AliError(Form("Impossible to read at : Channel %d",channel));
  return kFALSE;
}
//________________________________________________________________
void AliADCalibData::SetEnableTiming(Bool_t val, Int_t channel)
{
  // Set the channels enabled for Timing triggers
  if(channel>=0 && channel<16) fEnableTiming[channel] = val;
  else AliError(Form("Impossible to write at : Channel %d",channel));
}
//________________________________________________________________
Bool_t AliADCalibData::GetEnableTiming(Int_t channel) const
{
  // Get the channels enabled for Timing triggers
  if(channel>=0 && channel<16) return(fEnableTiming[channel]);
  else AliError(Form("Impossible to read at : Channel %d",channel));
  return kFALSE;
}
//________________________________________________________________
void AliADCalibData::SetEnableCharge(Bool_t val,Int_t board, Int_t channel)
{
  Int_t ch = AliADCalibData::GetOfflineChannelNumber(board,channel);
  // Set the channels enabled for Charge triggers
  if(ch>=0) fEnableCharge[ch] = val;
  else AliError(Form("Impossible to write at : Board %d ; Channel %d",board,channel));
}
//________________________________________________________________
void AliADCalibData::SetEnableTiming(Bool_t val,Int_t board, Int_t channel)
{
  Int_t ch = AliADCalibData::GetOfflineChannelNumber(board,channel);
  // Set the channels enabled for Timing triggers
  if(ch>=0) fEnableTiming[ch] = val;
  else AliError(Form("Impossible to write at : Board %d ; Channel %d",board,channel));
}
//________________________________________________________________
void AliADCalibData::SetTriggerSelected(UShort_t trigger, Int_t output)
{
  // Set the trigger selected on the outputs to CTP
  if(output>=0 && output<5) fTriggerSelected[output] = trigger;
  else AliError(Form("Trigger output number %d not valid",output));
}

//________________________________________________________________
void AliADCalibData::SetClk1Win1(UShort_t* clks)
{
  // Set Win clock of BB
  if(clks) for(int t=0; t<kNCIUBoards; t++) SetClk1Win1(clks[t],t);
  else AliError("Profil Clock1 Win1 Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetClk2Win1(UShort_t* clks)
{
  // Set Win clock of BB
  if(clks) for(int t=0; t<kNCIUBoards; t++) SetClk2Win1(clks[t],t);
  else AliError("Profil Clock2 Win1 Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetClk1Win1(UShort_t clk, Int_t board)
{
  // Set Win clock of BB
  if((board>=0) && (board<kNCIUBoards)) {
    fClk1Win1[board] = clk;
    if(!IsClkValid(clk)) AliWarning(Form("Profil Clock1 Win1 of board %d is not valid : %d",board,clk));
  }else {
    AliError(Form("Impossible to Write at Board %d",board));
  }
}
//________________________________________________________________
void AliADCalibData::SetClk2Win1(UShort_t clk, Int_t board)
{
  // Set Win clock of BB
  if((board>=0) && (board<kNCIUBoards)) {
    fClk2Win1[board] = clk;
    if(!IsClkValid(clk)) AliWarning(Form("Profil Clock2 Win1 of board %d is not valid : %d",board,clk));
  }else {
    AliError(Form("Impossible to Write at Board %d",board));
  }
}
//________________________________________________________________
void AliADCalibData::SetClk1Win2(UShort_t* clks)
{
  // Set Win clock of BG
  if(clks) for(int t=0; t<kNCIUBoards; t++) SetClk1Win2(clks[t],t);
  else AliError("Profil Clock1 Win2 Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetClk2Win2(UShort_t* clks)
{
  // Set Win clock of BG
  if(clks) for(int t=0; t<kNCIUBoards; t++) SetClk2Win2(clks[t],t);
  else AliError("Profil Clock2 Win2 Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetClk1Win2(UShort_t clk, Int_t board)
{
  // Set Win clock of BG
  if((board>=0) && (board<kNCIUBoards)) {
    fClk1Win2[board] = clk;
    if(!IsClkValid(clk)) AliWarning(Form("Profil Clock1 Win2 of board %d is not valid : %d",board,clk));
  }else {
    AliError(Form("Impossible to Write at Board %d",board));
  }
}
//________________________________________________________________
void AliADCalibData::SetClk2Win2(UShort_t clk, Int_t board)
{
  // Set Win clock of BG
  if((board>=0) && (board<kNCIUBoards)) {
    fClk2Win2[board] = clk;
    if(!IsClkValid(clk)) AliWarning(Form("Profil Clock2 Win2 of board %d is not valid : %d",board,clk));
  }else {
    AliError(Form("Impossible to Write at Board %d",board));
  }
}
//________________________________________________________________
void AliADCalibData::SetDelayClk1Win1(UShort_t* delays)
{
  // Set Delay for Win clock of BB
  if(delays) for(int t=0; t<kNCIUBoards; t++) SetDelayClk1Win1(delays[t],t);
  else AliError("Profil Clock1 Win1 Delays Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetDelayClk1Win1(UShort_t delay, Int_t board)
{
  // Set Delay for Win clock of BB
  if(delay>1023){
    AliWarning(Form("Profil Clock1 Win1 Delay of board %d should be less 1023 is currently %d. Truncated to the first 10 bits",board, delay));
    delay = delay & 0x3FF;
  }
  if((board>=0) && (board<kNCIUBoards))	fDelayClk1Win1[board] = delay;
  else AliError(Form("Trying to write out of the array Board = %d",board));
}
//________________________________________________________________
void AliADCalibData::SetDelayClk2Win1(UShort_t* delays)
{
  // Set Delay for Win clock of BB
  if(delays) for(int t=0; t<kNCIUBoards; t++) SetDelayClk2Win1(delays[t],t);
  else AliError("Profil Clock2 Win1 Delays Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetDelayClk2Win1(UShort_t delay, Int_t board)
{
  // Set Delay for Win clock of BB
  if(delay>1023){
    AliWarning(Form("Profil Clock2 Win1 Delay of board %d should be less 1023 is currently %d. Truncated to the first 10 bits",board, delay));
    delay = delay & 0x3FF;
  }
  if((board>=0) && (board<kNCIUBoards))	fDelayClk2Win1[board] = delay;
  else AliError(Form("Trying to write out of the array Board = %d",board));
}
//________________________________________________________________
void AliADCalibData::SetDelayClk1Win2(UShort_t* delays)
{
  // Set Delay for Win clock of BG
  if(delays) for(int t=0; t<kNCIUBoards; t++) SetDelayClk1Win2(delays[t],t);
  else AliError("Profil Clock1 Win2 Delays Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetDelayClk1Win2(UShort_t delay, Int_t board)
{
  // Set Delay for Win clock of BG
  if(delay>1023){
    AliWarning(Form("Profil Clock1 Win2 Delay of board %d should be less 1023 is currently %d. Truncated to the first 10 bits",board, delay));
    delay = delay & 0x3FF;
  }
  if((board>=0) && (board<kNCIUBoards))	fDelayClk1Win2[board] = delay;
  else AliError(Form("Trying to write out of the array Board = %d",board));
}
//________________________________________________________________
void AliADCalibData::SetDelayClk2Win2(UShort_t* delays)
{
  // Set Delay for Win clock of BG
  if(delays) for(int t=0; t<kNCIUBoards; t++) SetDelayClk2Win2(delays[t],t);
  else AliError("Profil Clock2 Win2 Delays Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetDelayClk2Win2(UShort_t delay, Int_t board)
{
  // Set Delay for Win clock of BG
  if(delay>1023){
    AliWarning(Form("Profil Clock2 Win2 Delay of board %d should be less 1023 is currently %d. Truncated to the first 10 bits",board, delay));
    delay = delay & 0x3FF;
  }
  if((board>=0) && (board<kNCIUBoards))	fDelayClk2Win2[board] = delay;
  else AliError(Form("Trying to write out of the array Board = %d",board));
}
//________________________________________________________________
void AliADCalibData::SetLatchWin1(UShort_t *latchs){
  // Set Latch Win clock for BB
  if(latchs) for(int t=0; t<kNCIUBoards; t++) SetLatchWin1(latchs[t],t);
  else AliError("Latch Win1 profil Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetLatchWin1(UShort_t latch, Int_t board)
{
  // Set Latch Win clock for BB
  if((board>=0) && (board<kNCIUBoards)) {
    fLatchWin1[board] = latch;
    if(!IsClkValid(latch)) AliWarning(Form("Latch Win1 of board %d is not valid : %d",board,latch));
  }else {
    AliError(Form("Impossible to Write at Board %d",board));
  }
}
//________________________________________________________________
void AliADCalibData::SetLatchWin2(UShort_t *latchs){
  // Set Latch Win clock for BG
  if(latchs) for(int t=0; t<kNCIUBoards; t++) SetLatchWin2(latchs[t],t);
  else AliError("Latch Win2 profil Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetLatchWin2(UShort_t latch, Int_t board)
{
  // Set Latch Win clock for BG
  if((board>=0) && (board<kNCIUBoards)) {
    fLatchWin2[board] = latch;
    if(!IsClkValid(latch)) AliWarning(Form("Latch Win2 of board %d is not valid : %d",board,latch));
  }else {
    AliError(Form("Impossible to Write at Board %d",board));
  }
}
//________________________________________________________________
void AliADCalibData::SetResetWin1(UShort_t *resets){
  // Set Reset Win clock for BB
  if(resets) for(int t=0; t<kNCIUBoards; t++) SetResetWin1(resets[t],t);
  else AliError("Reset Win1 profil Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetResetWin1(UShort_t reset, Int_t board)
{
  // Set Reset Win clock for BB
  if((board>=0) && (board<kNCIUBoards)) {
    fResetWin1[board] = reset;
    if(!IsClkValid(reset)) AliWarning(Form("Reset Win1 of board %d is not valid : %d",board,reset));
  }else {
    AliError(Form("Impossible to Write at Board %d",board));
  }
}
//________________________________________________________________
void AliADCalibData::SetResetWin2(UShort_t *resets){
  // Set Reset Win clock for BG
  if(resets)  for(int t=0; t<kNCIUBoards; t++) SetResetWin2(resets[t],t);
  else AliError("Reset Win2 profil Not defined.");
}
//________________________________________________________________
void AliADCalibData::SetResetWin2(UShort_t reset, Int_t board)
{
  // Set Reset Win clock for BG
  if((board>=0) && (board<kNCIUBoards)) {
    fResetWin2[board] = reset;
    if(!IsClkValid(reset)) AliWarning(Form("Reset Win2 of board %d is not valid : %d",board,reset));
  }else {
    AliError(Form("Impossible to Write at Board %d",board));
  }
}
//________________________________________________________________
void AliADCalibData::SetPedestalSubtraction(Bool_t *peds){
  // Set Pedestal Subtraction Parameter
  if(peds)  for(int t=0; t<kNCIUBoards; t++) SetPedestalSubtraction(peds[t],t);
  else AliError("Pedestal Subtraction Not defined.");
	
}
//________________________________________________________________
void AliADCalibData::SetPedestalSubtraction(Bool_t ped, Int_t board)
{
  // Set Pedestal Subtraction Parameter
  if((board>=0) && (board<kNCIUBoards)) fPedestalSubtraction[board] = ped;
  else AliError(Form("Board %d is not valid",board));
}

//________________________________________________________________
Bool_t	AliADCalibData::IsClkValid(UShort_t clock) const {
  // Check if the given clock has a valid profil.
  Bool_t word[5];
  Bool_t isValid = kTRUE;
  Short_t risingEdge = 0;
  Short_t fallingEdge = 0;
  for(int i=0 ; i<5 ; i++) word[i] = (clock >> i) & 0x1;
	
  if(word[0] != word[4]){
    if(word[4]) fallingEdge++;
    else risingEdge++;
  }	
  for(int i=1 ; i<5 ; i++){
    if(word[i] != word[i-1]) {
      if(word[i-1]) fallingEdge++;
      else risingEdge++;
    }
  }
  if((fallingEdge>1)||(risingEdge>1)) isValid = kFALSE;
  if(((risingEdge==0)&&(fallingEdge==0)) &&(!word[0]))  isValid = kFALSE;
  return isValid;
}

//________________________________________________________________
Int_t AliADCalibData::GetOfflineChannelNumber(Int_t board, Int_t channel)
{
  // Get the offline channel number from
  // the FEE board and channel indexes

  if (board < 0 || board >= 2) {
    AliErrorClass(Form("Wrong FEE board number: %d",board));
    return -1;
  }
  if (channel < 0 || channel >= 8) {
    AliErrorClass(Form("Wrong FEE channel number: %d",channel));
    return -1;
  }

  Int_t offCh = kOfflineChannel[channel+(8*board)];

  return offCh;
}
//________________________________________________________________
Int_t AliADCalibData::GetFEEChannelNumber(Int_t channel)
{
  // Get FEE channel number
  // from offline channel index
  if (channel >= 0 && channel < 16) return ((kOnlineChannel[channel] % 8));

  AliErrorClass(Form("Wrong channel index: %d",channel));
  return -1;
}
//________________________________________________________________
Int_t AliADCalibData::GetBoardNumber(Int_t channel)
{
  // Get FEE board number
  // from offline channel index
  Int_t OnChannel = kOnlineChannel[channel];
  if (OnChannel >= 0 && OnChannel < 8) return (0);
  if (OnChannel >= 8 && OnChannel < 16) return (1);

  AliErrorClass(Form("Wrong channel index: %d",channel));
  return -1;
}
//________________________________________________________________
void AliADCalibData::PrintConfig()
{

  printf("\n");  
  printf("======================================================\n");
  printf("=======================CCIU config====================\n");
  printf("======================================================\n");
  printf("\n"); 
 
  printf("Selected triggers = ");
  for(Int_t i = 0; i<5 ; i++) printf("%d ",GetTriggerSelected(i));
  printf("\n");

  printf("BBA thr = %d, BBC thr = %d\n",GetBBAThreshold(),GetBBCThreshold());
  printf("BGA thr = %d, BGC thr = %d\n",GetBGAThreshold(),GetBGCThreshold());
  printf("BBAforBG thr = %d, BBCforBG thr = %d\n",GetBBAForBGThreshold(),GetBBCForBGThreshold());

  
  printf("\n"); 
  printf("======================================================\n");
  printf("=======================CIUs config====================\n");
  printf("======================================================\n");
  printf("\n");
 
  for(Int_t CIUNumber = 0; CIUNumber < 2; ++CIUNumber) {
    printf("CIU= %d, Time Res= %f, Width Res= %f, RollOver= %d, TrigCountOffset= %d, SearchWin= %d, MatchWin= %d\n",
	   CIUNumber,
	   GetTimeResolution(CIUNumber),
	   GetWidthResolution(CIUNumber),
	   GetRollOver(CIUNumber),
	   GetTriggerCountOffset(CIUNumber),
	   GetSearchWindow(CIUNumber),
	   GetMatchWindow(CIUNumber));
	   
	  
  }

  printf("\n");
  printf("=======================Trigger windows================\n"); 
  printf("\n");
  Bool_t gatesOpen = kTRUE;
  for (Int_t CIUNumber = 0; CIUNumber < 2; ++CIUNumber) {
    if (GetDelayClk1Win1(CIUNumber)!=0 || GetDelayClk2Win1(CIUNumber)!=0 || GetDelayClk1Win2(CIUNumber)!=0 || GetDelayClk2Win2(CIUNumber)!=0) gatesOpen = kFALSE;
  }
  if(!gatesOpen){
    for(Int_t CIUNumber = 0; CIUNumber < 2; ++CIUNumber) {
      std::cout <<"CIU = "<<CIUNumber<<std::endl;
      std::cout << "P1 BB " << std::bitset<5>(GetClk1Win1(CIUNumber))<<  "   P1 BG " << std::bitset<5>(GetClk1Win2(CIUNumber))<< std::endl;
      std::cout << "P2 BB " << std::bitset<5>(GetClk2Win1(CIUNumber))<<  "   P2 BG " << std::bitset<5>(GetClk2Win2(CIUNumber))<< std::endl;
      std::cout << "  L   " << std::bitset<5>(GetLatchWin1(CIUNumber))<< "     L   " << std::bitset<5>(GetLatchWin2(CIUNumber))<< std::endl;
      std::cout << "  R   " << std::bitset<5>(GetResetWin1(CIUNumber))<< "     R   " << std::bitset<5>(GetResetWin2(CIUNumber))<< std::endl;
      printf("\n");
      std::cout << "BB Delay1 = " <<GetDelayClk1Win1(CIUNumber)<< "  BG Delay1 = " <<GetDelayClk1Win2(CIUNumber)<<std::endl;
      std::cout << "BB Delay2 = " <<GetDelayClk2Win1(CIUNumber)<< "  BG Delay2 = " <<GetDelayClk2Win2(CIUNumber)<<std::endl;
      printf("\n");
  
      AliADLogicalSignal clk1BB(GetClk1Win1(CIUNumber),GetDelayClk1Win1(CIUNumber),GetLatchWin1(CIUNumber),GetResetWin1(CIUNumber));
      AliADLogicalSignal clk2BB(GetClk2Win1(CIUNumber),GetDelayClk2Win1(CIUNumber),GetLatchWin1(CIUNumber),GetResetWin1(CIUNumber));
      AliADLogicalSignal *fBBGate = new AliADLogicalSignal(clk1BB & clk2BB);
	
      std::cout<<"BB window = ["<<fBBGate->GetStartTime()<<" , "<<fBBGate->GetStopTime()<<"]"<<std::endl;
	
      AliADLogicalSignal clk1BG(GetClk1Win2(CIUNumber),GetDelayClk1Win2(CIUNumber),GetLatchWin2(CIUNumber),GetResetWin2(CIUNumber));
      clk1BG.SetStartTime(clk1BG.GetStartTime()+2);
      clk1BG.SetStopTime(clk1BG.GetStopTime()+2);
      AliADLogicalSignal clk2BG(GetClk2Win2(CIUNumber),GetDelayClk2Win2(CIUNumber),GetLatchWin2(CIUNumber),GetResetWin2(CIUNumber));
      clk2BG.SetStartTime(clk2BG.GetStartTime()-2);
      clk2BG.SetStopTime(clk2BG.GetStopTime()-2);
      AliADLogicalSignal *fBGGate = new AliADLogicalSignal(clk1BG & clk2BG);
	
      std::cout<<"BG window = ["<<fBGGate->GetStartTime()<<" , "<<fBGGate->GetStopTime()<<"]"<<std::endl;
      printf("\n");	  
    }
  }
  else printf("Test window 25ns\n");

  printf("\n");  
  printf("======================================================\n");
  printf("====================Channels config===================\n");
  printf("======================================================\n");
  printf("\n"); 
 
  for(Int_t pmNumber = 0; pmNumber < 16; ++pmNumber) {
    printf("ChOff = %d, ChOn = %d, HV = %.1f, MIP = %.1f ADC, Dead = %s, DelayHit = %.2f, Thr_DCS = %.1f, Thr_Calib = %.1f\n",
	   pmNumber,
	   kOnlineChannel[pmNumber],
	   GetMeanHV(pmNumber),
	   GetADCperMIP(pmNumber),
	   IsChannelDead(pmNumber)?"yes":"no",
	   GetTimeOffset(pmNumber),
	   GetDiscriThr(pmNumber),
	   GetCalibDiscriThr(pmNumber)
	   );
  }
  
  printf("\n");
  printf("======================================================\n");
  printf("======================= Pedestal =====================\n");
  printf("======================================================\n");
  printf("\n");

  for(Int_t pmNumber = 0; pmNumber < 16; ++pmNumber) {
    for(Int_t integrator = 0; integrator < 2; ++integrator){
      if(integrator == 0)printf("ChOff = %d, ChOn = %d, Int = %d, Pedestal = %.3f, Width = %3f,", pmNumber, kOnlineChannel[pmNumber],integrator, GetPedestal(pmNumber+16*integrator),GetSigma(pmNumber+16*integrator));
      else printf(" Int = %d, Pedestal = %.3f, Width = %3f\n", integrator, GetPedestal(pmNumber+16*integrator),GetSigma(pmNumber+16*integrator));	
    }
  }



}

//________________________________________________________________
void AliADCalibData::PrintConfigShuttle()
{

  printf("\n");  
  printf("======================================================\n");
  printf("=======================CCIU config====================\n");
  printf("======================================================\n");
  printf("\n"); 
 
  printf("Selected triggers = ");
  for(Int_t i = 0; i<5 ; i++) printf("%d ",GetTriggerSelected(i));
  printf("\n");

  printf("BBA thr = %d, BBC thr = %d\n",GetBBAThreshold(),GetBBCThreshold());
  printf("BGA thr = %d, BGC thr = %d\n",GetBGAThreshold(),GetBGCThreshold());
  printf("BBAforBG thr = %d, BBCforBG thr = %d\n",GetBBAForBGThreshold(),GetBBCForBGThreshold());

  
  printf("\n"); 
  printf("======================================================\n");
  printf("=======================CIUs config====================\n");
  printf("======================================================\n");
  printf("\n");
 
  for(Int_t CIUNumber = 0; CIUNumber < 2; ++CIUNumber) {
    printf("CIU= %d, Time Res= %f, Width Res= %f, RollOver= %d, TrigCountOffset= %d, SearchWin= %d, MatchWin= %d\n",
	   CIUNumber,
	   GetTimeResolution(CIUNumber),
	   GetWidthResolution(CIUNumber),
	   GetRollOver(CIUNumber),
	   GetTriggerCountOffset(CIUNumber),
	   GetSearchWindow(CIUNumber),
	   GetMatchWindow(CIUNumber));
	   
	  
  }

  printf("\n");
  printf("=======================Trigger windows================\n"); 
  printf("\n");
  Bool_t gatesOpen = kTRUE;
  for (Int_t CIUNumber = 0; CIUNumber < 2; ++CIUNumber) {
    if (GetDelayClk1Win1(CIUNumber)!=0 || GetDelayClk2Win1(CIUNumber)!=0 || GetDelayClk1Win2(CIUNumber)!=0 || GetDelayClk2Win2(CIUNumber)!=0) gatesOpen = kFALSE;
  }
  if(!gatesOpen){
    for(Int_t CIUNumber = 0; CIUNumber < 2; ++CIUNumber) {
      std::cout <<"CIU = "<<CIUNumber<<std::endl;
      std::cout << "P1 BB " << std::bitset<5>(GetClk1Win1(CIUNumber))<<  "   P1 BG " << std::bitset<5>(GetClk1Win2(CIUNumber))<< std::endl;
      std::cout << "P2 BB " << std::bitset<5>(GetClk2Win1(CIUNumber))<<  "   P2 BG " << std::bitset<5>(GetClk2Win2(CIUNumber))<< std::endl;
      std::cout << "  L   " << std::bitset<5>(GetLatchWin1(CIUNumber))<< "     L   " << std::bitset<5>(GetLatchWin2(CIUNumber))<< std::endl;
      std::cout << "  R   " << std::bitset<5>(GetResetWin1(CIUNumber))<< "     R   " << std::bitset<5>(GetResetWin2(CIUNumber))<< std::endl;
      printf("\n");
      std::cout << "BB Delay1 = " <<GetDelayClk1Win1(CIUNumber)<< "  BG Delay1 = " <<GetDelayClk1Win2(CIUNumber)<<std::endl;
      std::cout << "BB Delay2 = " <<GetDelayClk2Win1(CIUNumber)<< "  BG Delay2 = " <<GetDelayClk2Win2(CIUNumber)<<std::endl;  
    }
  }
  else printf("Test window 25ns\n");

  printf("\n");  
  printf("======================================================\n");
  printf("====================Channels config===================\n");
  printf("======================================================\n");
  printf("\n"); 
 
  for(Int_t pmNumber = 0; pmNumber < 16; ++pmNumber) {
    printf("ChOff = %d, ChOn = %d, HV = %.1f, Dead = %s, DelayHit = %.2f, Thr_DCS = %.1f\n",
	   pmNumber,
	   kOnlineChannel[pmNumber],
	   GetMeanHV(pmNumber),
	   IsChannelDead(pmNumber)?"yes":"no",
	   GetTimeOffset(pmNumber),
	   GetDiscriThr(pmNumber)
	   );
  }
  
  printf("\n");
  printf("======================================================\n");
  printf("======================= Pedestal =====================\n");
  printf("======================================================\n");
  printf("\n");

  for(Int_t pmNumber = 0; pmNumber < 16; ++pmNumber) {
    for(Int_t integrator = 0; integrator < 2; ++integrator){
      if(integrator == 0)printf("ChOff = %d, ChOn = %d, Int = %d, Pedestal = %.3f, Width = %3f,", pmNumber, kOnlineChannel[pmNumber],integrator, GetPedestal(pmNumber+16*integrator),GetSigma(pmNumber+16*integrator));
      else printf(" Int = %d, Pedestal = %.3f, Width = %3f\n", integrator, GetPedestal(pmNumber+16*integrator),GetSigma(pmNumber+16*integrator));	
    }
  }

}
