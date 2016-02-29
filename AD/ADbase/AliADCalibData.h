#ifndef ALIADCALIBDATA_H
#define ALIADCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
#include "AliADConst.h"

class AliADDataDCS;


class AliADCalibData: public TNamed {

 public:
  AliADCalibData();
  AliADCalibData(const char* name);
  
  AliADCalibData(const AliADCalibData &calibda);
  AliADCalibData& operator= (const AliADCalibData &calibda);
  virtual ~AliADCalibData();
  void Reset();
  void PrintConfig();
  void PrintConfigShuttle();
  void FillDCSData(AliADDataDCS * data);

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
  Float_t  GetADCperMIP(Int_t channel);
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
  Float_t  GetCalibDiscriThr(Int_t channel);
  
  UShort_t * GetClk1Win1() const {return (UShort_t*)fClk1Win1;};
  UShort_t GetClk1Win1(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fClk1Win1[board]:0);};
  UShort_t * GetClk2Win1() const {return (UShort_t*)fClk2Win1;};
  UShort_t GetClk2Win1(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fClk2Win1[board]:0);};

  UShort_t * GetClk1Win2() const {return (UShort_t*)fClk1Win2;};
  UShort_t GetClk1Win2(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fClk1Win2[board]:0);};
  UShort_t * GetClk2Win2() const {return (UShort_t*)fClk2Win2;};
  UShort_t GetClk2Win2(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fClk2Win2[board]:0);};

  UShort_t * GetDelayClk1Win1() const {return (UShort_t*)fDelayClk1Win1;};
  UShort_t GetDelayClk1Win1(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fDelayClk1Win1[board]:0);};
  UShort_t * GetDelayClk2Win1() const {return (UShort_t*)fDelayClk2Win1;};
  UShort_t GetDelayClk2Win1(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fDelayClk2Win1[board]:0);};
  
  UShort_t * GetDelayClk1Win2() const {return (UShort_t*)fDelayClk1Win2;};
  UShort_t GetDelayClk1Win2(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fDelayClk1Win2[board]:0);};
  UShort_t * GetDelayClk2Win2() const {return (UShort_t*)fDelayClk2Win2;};
  UShort_t GetDelayClk2Win2(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fDelayClk2Win2[board]:0);};
  
  UShort_t * GetLatchWin1() const {return (UShort_t*)fLatchWin1;};
  UShort_t GetLatchWin1(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fLatchWin1[board]:0);};
  UShort_t * GetLatchWin2() const {return (UShort_t*)fLatchWin2;};
  UShort_t GetLatchWin2(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fLatchWin2[board]:0);};
  
  UShort_t * GetResetWin1() const {return (UShort_t*)fResetWin1;};
  UShort_t GetResetWin1(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fResetWin1[board]:0);};
  UShort_t * GetResetWin2() const {return (UShort_t*)fResetWin2;};
  UShort_t GetResetWin2(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fResetWin2[board]:0);};
  
  Bool_t * GetPedestalSubtraction() const {return (Bool_t*) fPedestalSubtraction;};
  Bool_t GetPedestalSubtraction(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fPedestalSubtraction[board]:0);};

  UShort_t GetBBAThreshold() const {return fBBAThreshold;};
  UShort_t GetBBCThreshold() const {return fBBCThreshold;};

  UShort_t GetBGAThreshold() const {return fBGAThreshold;};
  UShort_t GetBGCThreshold() const {return fBGCThreshold;};

  UShort_t GetBBAForBGThreshold() const {return fBBAForBGThreshold;};
  UShort_t GetBBCForBGThreshold() const {return fBBCForBGThreshold;};
  
  UShort_t GetMultADAThrLow() const {return fMultADAThrLow;};
  UShort_t GetMultADAThrHigh() const {return fMultADAThrHigh;};

  UShort_t GetMultADCThrLow() const {return fMultADCThrLow;};
  UShort_t GetMultADCThrHigh() const {return fMultADCThrHigh;};

  UShort_t GetTriggerSelected(Int_t output) const {return ((output>=0 && output<5)?fTriggerSelected[output]:0);};
  
  Bool_t GetEnableCharge(Int_t channel) const;
  Bool_t GetEnableTiming(Int_t channel) const;
  UShort_t GetOnlinePedestal(Int_t integrator, Int_t channel) const;
  UShort_t GetOnlinePedestalCut(Int_t integrator, Int_t channel) const;


  static Int_t GetBoardNumber(Int_t channel);
  static Int_t GetFEEChannelNumber(Int_t channel);
  static Int_t GetOfflineChannelNumber(Int_t board, Int_t channel);
  
  Float_t  GetLightYields(Int_t channel);

  Float_t *GetPMGainsA() const { return fPMGainsA; }
  Float_t *GetPMGainsB() const { return fPMGainsB; }

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
  
  void	   SetADCperMIP(Int_t nADCperMIP);

  void SetClk1Win1(UShort_t* clks);
  void SetClk1Win1(UShort_t clk, Int_t board);
  void SetClk2Win1(UShort_t* clks);
  void SetClk2Win1(UShort_t clk, Int_t board);
  
  void SetClk1Win2(UShort_t* clks);
  void SetClk1Win2(UShort_t clk, Int_t board);
  void SetClk2Win2(UShort_t* clks);
  void SetClk2Win2(UShort_t clk, Int_t board);
  
  void SetDelayClk1Win1(UShort_t* delays);
  void SetDelayClk1Win1(UShort_t delay, Int_t board);
  void SetDelayClk2Win1(UShort_t* delays);
  void SetDelayClk2Win1(UShort_t delay, Int_t board);
  
  void SetDelayClk1Win2(UShort_t* delays);
  void SetDelayClk1Win2(UShort_t delay, Int_t board);
  void SetDelayClk2Win2(UShort_t* delays);
  void SetDelayClk2Win2(UShort_t delay, Int_t board);
  
  void SetLatchWin1(UShort_t *latchs);
  void SetLatchWin1(UShort_t latch, Int_t board);
  void SetLatchWin2(UShort_t *latchs);
  void SetLatchWin2(UShort_t latch, Int_t board);
  
  void SetResetWin1(UShort_t *resets);
  void SetResetWin1(UShort_t reset, Int_t board);
  void SetResetWin2(UShort_t *resets);
  void SetResetWin2(UShort_t reset, Int_t board);
  
  void SetPedestalSubtraction(Bool_t *peds);
  void SetPedestalSubtraction(Bool_t ped, Int_t board);
  
  void SetBBAThreshold(UShort_t th) {fBBAThreshold = th;};
  void SetBBCThreshold(UShort_t th) {fBBCThreshold = th;};

  void SetBGAThreshold(UShort_t th) {fBGAThreshold = th;};
  void SetBGCThreshold(UShort_t th) {fBGCThreshold = th;};

  void SetBBAForBGThreshold(UShort_t th) {fBBAForBGThreshold = th;};
  void SetBBCForBGThreshold(UShort_t th) {fBBCForBGThreshold = th;};
  
  
  void SetMultADAThrLow(UShort_t th) {fMultADAThrLow = th;};
  void SetMultADAThrHigh(UShort_t th) {fMultADAThrHigh = th;};
  
  void SetMultADCThrLow(UShort_t th) {fMultADCThrLow = th;};
  void SetMultADCThrHigh(UShort_t th) {fMultADCThrHigh = th;};
  
  void SetTriggerSelected(UShort_t trigger, Int_t output);
  
  void SetEnableCharge(Bool_t val, Int_t board, Int_t channel);
  void SetEnableTiming(Bool_t val, Int_t board, Int_t channel);
  void SetEnableCharge(Bool_t val, Int_t channel);
  void SetEnableTiming(Bool_t val, Int_t channel);
  void SetOnlinePedestal(UShort_t val, Int_t integrator, Int_t board, Int_t channel);
  void SetOnlinePedestalCut(UShort_t val, Int_t integrator, Int_t board, Int_t channel);
  void SetOnlinePedestal(UShort_t val, Int_t integrator, Int_t channel);
  void SetOnlinePedestalCut(UShort_t val, Int_t integrator, Int_t channel);

 protected:
  void     InitLightYields();
  void     InitPMGains();
  void     InitCalibThresholds();
  Bool_t   IsClkValid(UShort_t clock) const;

  Float_t  fPedestal[32];     // Mean pedestal values - used offline
  Float_t  fSigma[32];        // Sigmas of pedestal peaks - used offline
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
  
  UShort_t fClk1Win1[kNCIUBoards]; //Profil of the Clock 1  of the Window 1 (BB window)
  UShort_t fClk2Win1[kNCIUBoards]; //Profil of the Clock 2  of the Window 1 (BB window)
  UShort_t fClk1Win2[kNCIUBoards]; //Profil of the Clock 1  of the Window 2 (BG window)
  UShort_t fClk2Win2[kNCIUBoards]; //Profil of the Clock 2  of the Window 2 (BG window)
  UShort_t fDelayClk1Win1[kNCIUBoards]; // Delays of the Clock 1  of the Window 1 (BB window)
  UShort_t fDelayClk2Win1[kNCIUBoards]; // Delays of the Clock 2 of the Window 1 (BB window)
  UShort_t fDelayClk1Win2[kNCIUBoards]; // Delays of the Clock 1  of the Window 2 (BG window)
  UShort_t fDelayClk2Win2[kNCIUBoards]; // Delays of the Clock 2 of the Window 2 (BG window)
  UShort_t fLatchWin1[kNCIUBoards]; //Profil of the Clock of the Latch signal of Window 1 (BB window)
  UShort_t fLatchWin2[kNCIUBoards]; //Profil of the Clock of the Latch signal of Window 2 (BG window)
  UShort_t fResetWin1[kNCIUBoards]; //Profil of the Clock of the Reset signal of Window 1 (BB window)
  UShort_t fResetWin2[kNCIUBoards]; //Profil of the Clock of the Reset signal of Window 2 (BG window)
  Bool_t   fPedestalSubtraction[kNCIUBoards]; // Flag to en(dis)able pedestal subtraction before centrality trigger calculation
  UShort_t fBBAThreshold;  // Minimum bias Threshold in number of channel hit for ADA
  UShort_t fBBCThreshold;  // Minimum bias Threshold in number of channel hit for ADC
  UShort_t fBGAThreshold;  // Beam Gas Threshold in number of channel hit for ADA
  UShort_t fBGCThreshold;  // Beam Gas Threshold in number of channel hit for ADC
  UShort_t fBBAForBGThreshold;  // BBA threshold for Beam Gas triggers (i.e. BBA and BGC)
  UShort_t fBBCForBGThreshold;  // BBC threshold for Beam Gas triggers (i.e. BBC and BGA)
  UShort_t fMultADAThrLow;  // Threshold used for multiplicity triggers (i.e. MTA and MTC)
  UShort_t fMultADAThrHigh; // Threshold used for multiplicity triggers (i.e. MTA and MTC)
  UShort_t fMultADCThrLow;  // Threshold used for multiplicity triggers (i.e. MTA and MTC)
  UShort_t fMultADCThrHigh; // Threshold used for multiplicity triggers (i.e. MTA and MTC)
  UShort_t fTriggerSelected[5]; // Triggers selected on the 5 outputs to CTP
  Bool_t   fEnableCharge[16]; // Flag to know is a channel is participating to the Charge triggers
  Bool_t   fEnableTiming[16]; // Flag to know is a channel is participating to the Timing triggers
  UShort_t fPedestalOdd[16]; // Pedestals for the Odd integrators
  UShort_t fPedestalEven[16]; // Pedestals for the Even integrators
  UShort_t fPedestalCutOdd[16]; // Pedestals Cut for the Odd integrators
  UShort_t fPedestalCutEven[16]; // Pedestals Cut for the Even integrators

  Float_t  fDiscriThr[16];     // Discriminator thresholds

  Float_t *fLightYields;       //! Light Yields channel by channel (read from separate OCDB entry)
  Float_t *fPMGainsA;          //! PM gain factors channel by channel (read from separate OCDB entry)
  Float_t *fPMGainsB;          //! PM gain factors channel by channel (read from separate OCDB entry)
  Float_t *fThrCalibA;	       //! Thereshold calibration parameters channel by channel (read from separate OCDB entry)
  Float_t *fThrCalibB;	       //! Thereshold calibration parameters channel by channel (read from separate OCDB entry)

  ClassDef(AliADCalibData,2)    // AD Calibration data
};

#endif

