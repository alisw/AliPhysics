#ifndef ALIADCALIBDATA_H
#define ALIADCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
#include "AliADConst.h"
class AliADCalibData: public TNamed {

 public:
  AliADCalibData();
  AliADCalibData(const char* name);
  
  AliADCalibData(const AliADCalibData &calibda);
  AliADCalibData& operator= (const AliADCalibData &calibda);
  virtual ~AliADCalibData();
  void Reset();

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

  static Int_t GetBoardNumber(Int_t channel);
  
  Float_t  GetLightYields(Int_t channel);

  Float_t *GetPMGainsA() const { return fPMGainsA; }
  Float_t *GetPMGainsB() const { return fPMGainsB; }

  Float_t  GetCalibDiscriThr(Int_t channel, Bool_t scaled);


 protected:
  void     InitLightYields();
  void     InitPMGains();

  Float_t  fPedestal[32];     // Mean pedestal values
  Float_t  fSigma[32];        // Sigmas of pedestal peaks
  Float_t  fADCmean[32];      // ADC mean values
  Float_t  fADCsigma[32];     // ADC sigma values
  Float_t  fMeanHV[16];        // Mean PMT HV needed to compute MIP value
  Float_t  fWidthHV[16];       // Width of the PMT HV
  
  Float_t  fTimeOffset[16];    // Time offsets of the TDC
  Float_t  fTimeGain[16];      // Gain factors of the TDC
  Bool_t   fDeadChannel[16];   // List of dead channels
  Float_t  fTimeResolution[kNCIUBoards]; // Time Resolution of the TDC (ns / channel)
  Float_t  fWidthResolution[kNCIUBoards]; // Time Width Resolution of the TDC (ns / channel)

  UInt_t   fMatchWindow[kNCIUBoards]; // HPTDC matching window (25ns units)
  UInt_t   fSearchWindow[kNCIUBoards];// HPTDC search window (25ns units)
  UInt_t   fTriggerCountOffset[kNCIUBoards]; // HPTDC trigger count offset (25ns units)
  UInt_t   fRollOver[kNCIUBoards]; // HPTDC roll-over (25ns units)

  Float_t  fDiscriThr[16];     // Discriminator thresholds

  Float_t *fLightYields;       //! Light Yields channel by channel (read from separate OCDB entry)
  Float_t *fPMGainsA;          //! PM gain factors channel by channel (read from separate OCDB entry)
  Float_t *fPMGainsB;          //! PM gain factors channel by channel (read from separate OCDB entry)

  ClassDef(AliADCalibData,1)    // AD Calibration data
};

#endif

