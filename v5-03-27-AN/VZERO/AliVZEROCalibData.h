#ifndef ALIVZEROCALIBDATA_H
#define ALIVZEROCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//                                            // 
//  class for VZERO calibration               //
//                                            //
////////////////////////////////////////////////

#include "TNamed.h"

class AliVZERODataDCS;

class AliVZEROCalibData: public TNamed {

 public:
  enum { kNCIUBoards = 8 };
  
  AliVZEROCalibData();
  AliVZEROCalibData(const char* name);
  AliVZEROCalibData(const AliVZEROCalibData &calibda);
  AliVZEROCalibData& operator= (const AliVZEROCalibData &calibda);
  virtual ~AliVZEROCalibData();
  void Reset();
  void FillDCSData(AliVZERODataDCS * data);

  Float_t  GetPedestal(Int_t channel)   const {return fPedestal[channel];}
  Float_t* GetPedestal()   const {return (float*)fPedestal;}
  Float_t  GetSigma(Int_t channel)   const {return fSigma[channel];}
  Float_t* GetSigma()   const {return (float*)fSigma;}
  Float_t  GetADCmean(Int_t channel)	const {return fADCmean[channel];}
  Float_t* GetADCmean()   const {return (float*)fADCmean;}
  Float_t  GetADCsigma(Int_t channel)	const {return fADCsigma[channel];}
  Float_t* GetADCsigma()   const {return (float*)fADCsigma;}
  Float_t  GetMeanHV(Int_t channel)	const {return fMeanHV[channel];}
  Float_t* GetMeanHV()   const {return (float*)fMeanHV;} 
  Float_t  GetWidthHV(Int_t channel)	const {return fWidthHV[channel];}
  Float_t* GetWidthHV()   const {return (float*)fWidthHV;}
  Bool_t   IsChannelDead(Int_t channel)	const {return fDeadChannel[channel];}
  Bool_t*  GetDeadMap()   const {return (bool*)fDeadChannel;} 
   
  Float_t  GetGain(Int_t channel);
  Float_t  GetTimeOffset(Int_t channel)	const {return fTimeOffset[channel];}
  Float_t* GetTimeOffset()   const {return (float*)fTimeOffset;}
  Float_t  GetTimeGain(Int_t channel)	const {return fTimeGain[channel];}
  Float_t* GetTimeGain()   const {return (float*)fTimeGain;}

  Float_t* GetTimeResolution() const {return (Float_t*) fTimeResolution;};
  Float_t  GetTimeResolution(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fTimeResolution[board]:0);};

  Float_t* GetWidthResolution() const {return (Float_t*) fWidthResolution;};
  Float_t  GetWidthResolution(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fWidthResolution[board]:0);};

  const UInt_t*  GetMatchWindow() const { return fMatchWindow; }
  UInt_t   GetMatchWindow(Int_t board) const { return ((board>=0 && board<kNCIUBoards)?fMatchWindow[board]:0); }
  const UInt_t*  GetSearchWindow() const { return fSearchWindow; }
  UInt_t   GetSearchWindow(Int_t board) const { return ((board>=0 && board<kNCIUBoards)?fSearchWindow[board]:0); }
  const UInt_t*  GetTriggerCountOffset() const { return fTriggerCountOffset; }
  UInt_t   GetTriggerCountOffset(Int_t board) const { return ((board>=0 && board<kNCIUBoards)?fTriggerCountOffset[board]:0); }
  const UInt_t*  GetRollOver() const { return fRollOver; }
  UInt_t   GetRollOver(Int_t board) const { return ((board>=0 && board<kNCIUBoards)?fRollOver[board]:0); }

  Float_t  GetDiscriThr(Int_t channel)	const {return fDiscriThr[channel];}
  Float_t* GetDiscriThr()   const {return (Float_t*)fDiscriThr;}
    
  void     SetPedestal(Float_t val, Int_t channel) {fPedestal[channel]=val;}
  void     SetPedestal(const Float_t* Pedestal);
  void     SetSigma(Float_t val, Int_t channel) {fSigma[channel]=val;}
  void     SetSigma(const Float_t* Sigma);
  void 	   SetADCmean(Float_t val, Int_t channel) {fADCmean[channel]=val;}
  void 	   SetADCmean(const Float_t* ADCmean);  
  void 	   SetADCsigma(Float_t val, Int_t channel) {fADCsigma[channel]=val;}
  void 	   SetADCsigma(const Float_t* ADCsigma);
  void     SetMeanHV(Float_t val, Int_t channel) {fMeanHV[channel]=val;}
  void     SetMeanHV(const Float_t* MeanHV);  
  void     SetWidthHV(Float_t val, Int_t channel) {fWidthHV[channel]=val;}
  void     SetWidthHV(const Float_t* WidthHV); 
  void     SetDeadChannel(Bool_t val, Int_t channel) {fDeadChannel[channel]=val;}
  void     SetDeadMap(const Bool_t* deadMap);  
   
  void     SetTimeOffset(Float_t val, Int_t board, Int_t channel);
  void     SetTimeOffset(const Float_t* TimeOffset);
  void     SetTimeGain(Float_t val, Int_t channel) {fTimeGain[channel]=val;}
  void     SetTimeGain(const Float_t* TimeGain);
  
  void 	   SetParameter(TString name, Int_t val);
  void     SetTimeResolution(UShort_t *resols);
  void     SetTimeResolution(UShort_t resol, Int_t board);
  void     SetWidthResolution(UShort_t *resols);
  void     SetWidthResolution(UShort_t resol, Int_t board);

  void     SetMatchWindow(UInt_t *windows);
  void     SetMatchWindow(UInt_t window, Int_t board);
  void     SetSearchWindow(UInt_t *windows);
  void     SetSearchWindow(UInt_t window, Int_t board);
  void     SetTriggerCountOffset(UInt_t *offsets);
  void     SetTriggerCountOffset(UInt_t offset, Int_t board);
  void     SetRollOver(UInt_t *offsets);
  void     SetRollOver(UInt_t offset, Int_t board);

  void     SetDiscriThr(Float_t thr, Int_t board, Int_t channel);
  void     SetDiscriThr(const Float_t* thresholds);

  Float_t  GetMIPperADC(Int_t channel);
  Float_t  GetHV(Int_t channel, Float_t adcPerMip);

  static Int_t GetOfflineChannelNumber(Int_t board, Int_t channel);
  static Int_t GetBoardNumber(Int_t channel);
  static Int_t GetFEEChannelNumber(Int_t channel);

  Float_t  GetLightYields(Int_t channel);

  Float_t *GetPMGainsA() const { return fPMGainsA; }
  Float_t *GetPMGainsB() const { return fPMGainsB; }

  Float_t  GetCalibDiscriThr(Int_t channel, Bool_t scaled);

 protected:
  void     InitLightYields();
  void     InitPMGains();

  Float_t  fPedestal[128];     // Mean pedestal values
  Float_t  fSigma[128];        // Sigmas of pedestal peaks
  Float_t  fADCmean[128];      // ADC mean values
  Float_t  fADCsigma[128];     // ADC sigma values
  Float_t  fMeanHV[64];        // Mean PMT HV needed to compute MIP value
  Float_t  fWidthHV[64];       // Width of the PMT HV
  
  Float_t  fTimeOffset[64];    // Time offsets of the TDC
  Float_t  fTimeGain[64];      // Gain factors of the TDC
  Bool_t   fDeadChannel[64];   // List of dead channels
  Float_t  fTimeResolution[kNCIUBoards]; // Time Resolution of the TDC (ns / channel)
  Float_t  fWidthResolution[kNCIUBoards]; // Time Width Resolution of the TDC (ns / channel)

  UInt_t   fMatchWindow[kNCIUBoards]; // HPTDC matching window (25ns units)
  UInt_t   fSearchWindow[kNCIUBoards];// HPTDC search window (25ns units)
  UInt_t   fTriggerCountOffset[kNCIUBoards]; // HPTDC trigger count offset (25ns units)
  UInt_t   fRollOver[kNCIUBoards]; // HPTDC roll-over (25ns units)

  Float_t  fDiscriThr[64];     // Discriminator thresholds

  Float_t *fLightYields;       //! Light Yields channel by channel (read from separate OCDB entry)
  Float_t *fPMGainsA;          //! PM gain factors channel by channel (read from separate OCDB entry)
  Float_t *fPMGainsB;          //! PM gain factors channel by channel (read from separate OCDB entry)

  ClassDef(AliVZEROCalibData,8)    // VZERO Calibration data
};

#endif
