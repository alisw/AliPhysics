#ifndef ALIITSCALIBRATIONSSD_H
#define ALIITSCALIBRATIONSSD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 
#include "AliITSCalibration.h"
#include "AliITSNoiseSSDv2.h"
#include "AliITSPedestalSSDv2.h"
#include "AliITSGainSSDv2.h"
#include "AliITSBadChannelsSSDv2.h"
#include "TArrayF.h"
#include "TArrayI.h"

/* $Id$ */
//////////////////////////////////////////////
// Calibration class for SSD                   //
//                                          //
//////////////////////////////////////////////
class AliITSCalibrationSSD : public AliITSCalibration {

 public:
    AliITSCalibrationSSD();
    AliITSCalibrationSSD(const char *dataType);
    virtual ~AliITSCalibrationSSD();

    Float_t GetNoiseP(Int_t n) {return fNoise->GetNoiseP(fModule,n); }
    Float_t GetNoiseN(Int_t n) {return fNoise->GetNoiseN(fModule,n); }
    void SetNoise( AliITSNoiseSSDv2* noise) {fNoise=noise;}

    Float_t GetPedestalP(Int_t n) {return fPedestal->GetPedestalP(fModule,n); }
    Float_t GetPedestalN(Int_t n) {return fPedestal->GetPedestalN(fModule,n); }
    void SetPedestal( AliITSPedestalSSDv2* pedestal) {fPedestal=pedestal;}

    Float_t GetGainP(Int_t n) {return fGain->GetGainP(fModule,n); }
    Float_t GetGainN(Int_t n) {return fGain->GetGainN(fModule,n); }
    void SetGainP(Int_t n, Float_t value) {fGain->AddGainP(fModule,n,value);}
    void SetGainN(Int_t n, Float_t value) {fGain->AddGainN(fModule,n,value);}
    void SetGain( AliITSGainSSDv2* gain) {fGain=gain;}

    void   SetBad() {
      fIsBad = kTRUE;
      for(Int_t i=0;i<fgkChipsPerModule;i++) fIsChipBad[i]=kTRUE;
    }
    virtual Bool_t IsBad() const { return fIsBad; }
    void   SetChipBad(Int_t nChip) {
      fIsChipBad[nChip] = kTRUE;
    }
    void FillBadChipMap();
    virtual Bool_t IsChipBad(Int_t nChip) const {
      return fIsChipBad[nChip];
    }
    Int_t ChipsPerModule() const{return fgkChipsPerModule;} // # chips/module
    Int_t ChannelsPerChip() const{ return fgkChannelsPerChip;}// #channels/chip

    void SetBadChannels( AliITSBadChannelsSSDv2* badchannels) {
      fBadChannels=badchannels;}
    Char_t GetBadPChannel(Int_t n) {
      return fBadChannels->GetBadChannelP(fModule,n); }
    Char_t GetBadNChannel(Int_t n) {
      return fBadChannels->GetBadChannelN(fModule,n); }
    Bool_t IsPChannelBad(Int_t n) {
      return fBadChannels->GetBadChannelP(fModule,n)&1; }
    Bool_t IsNChannelBad(Int_t n) {
      return fBadChannels->GetBadChannelN(fModule,n)&1; }

    //
    virtual  void   SetNDetParam(Int_t npar) {
	// set number of param
	fNPar=npar;
    }
    virtual  void   SetDetParam(Double_t *par);

    // Parameters options
    virtual Int_t   NDetParam() const {
	// number of param
	return fNPar;
    }
    virtual void    GetDetParam(Double_t *dpar) const; 

    virtual void    SetSigmaSpread(Double_t, Double_t) {
      NotImplemented("SetSigmaSpread");}
    
    virtual void    SigmaSpread(Double_t &, Double_t &) const {
      NotImplemented("SetSigmaSpread");}
    
    void SetModule(Int_t mod){fModule = mod;}
    void SetModuleIndex(Int_t modId){
      fModule=modId-500; // temporary patch, 500 is n. of SPD+SDD modules
      FillBadChipMap();
    }
    void SetKeVperADC(Double_t a=86.4/120.){fKeVperADC = a;}
    Double_t ADCToKeV(Double_t adc) const {return adc*fKeVperADC;}

    void SetSSDADCpereV(Double_t a=(120./24888.9)*(85.0/77.0)){fSSDADCpereV = a;}
    Double_t GetSSDDEvToADC(Double_t eV) const {return eV*fSSDADCpereV;}
    Int_t GetSSDIEvToADC(Double_t eV) const { 
                                  return ((Int_t) GetSSDDEvToADC(eV)); }


 protected:

    static const Int_t fgkChipsPerModule  = 12;    // Number of chips/module
    static const Int_t fgkChannelsPerChip = 128;   // Number of channels/chip

    static const Int_t fgkNParDefault; // default for fNPar

    Int_t fModule;          //! module number (range 0 -> 1697)
    
    Int_t   fNPar;            // Number of detector param 
    Double_t *fDetPar;         //[fNPar] Array of parameters

    AliITSNoiseSSDv2 *fNoise;
    AliITSPedestalSSDv2 *fPedestal;
    AliITSGainSSDv2 *fGain;
    AliITSBadChannelsSSDv2 *fBadChannels;

    Bool_t   fIsBad;                         // module is dead or alive ?
    Bool_t   fIsChipBad[fgkChipsPerModule];  // chip is dead or alive ?

    Double_t fSSDADCpereV;    // Constant to convert eV to ADC for SSD.
    Double_t fKeVperADC;       // Constant to convert ADC to keV

 private:
    AliITSCalibrationSSD(const AliITSCalibrationSSD &source); // copy constructor
    AliITSCalibrationSSD& operator=(const AliITSCalibrationSSD &source); // ass. op.

    ClassDef(AliITSCalibrationSSD,6) //Response class for SSD
      };
#endif
