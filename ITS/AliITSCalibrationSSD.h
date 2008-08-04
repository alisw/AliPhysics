#ifndef ALIITSCALIBRATIONSSD_H
#define ALIITSCALIBRATIONSSD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 
#include "AliITSCalibration.h"
#include "AliITSNoiseSSDv2.h"
#include "AliITSPedestalSSDv2.h"
#include "AliITSGainSSDv2.h"
#include "AliITSBadChannelsSSDv2.h"
#include "AliITSresponseSSD.h"
#include "TArrayF.h"
#include "TArrayI.h"

/* $Id$ */
//////////////////////////////////////////////
// Response class for SSD                   //
//                                          //
//////////////////////////////////////////////
class AliITSCalibrationSSD : public AliITSCalibration {

 public:
    AliITSCalibrationSSD();
    AliITSCalibrationSSD(const char *dataType);
    virtual ~AliITSCalibrationSSD();

    virtual  void   SetNoiseParam(Double_t , Double_t ) {
      NotImplemented("SetNoiseParam");}
    
    virtual void    GetNoiseParam(Double_t &, Double_t &) const {
      NotImplemented("SetNoiseParam");}

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
    
    
    virtual void   SetThresholds(Double_t /* a */, Double_t /* b */)
      {NotImplemented("SetThresholds");}
    virtual void   Thresholds(Double_t & /* a */, Double_t & /* b */) const 
      {NotImplemented("Thresholds");}
  
    virtual void SetADCpereV(Double_t a=120./24888.9) {((AliITSresponseSSD*)fResponse)->SetADCpereV(a);}
    virtual Double_t GetDEvToADC(Double_t eV) const {return ((AliITSresponseSSD*)fResponse)->DEvToADC(eV);}
    virtual Int_t IEvToADC(Double_t eV) const {return ((AliITSresponseSSD*)fResponse)->IEvToADC(eV);} 

    virtual void SetKeVperADC(Double_t a=86.4/120.0) {((AliITSresponseSSD*)fResponse)->SetKeVperADC(a);}
    virtual Double_t ADCToKeV(Double_t adc) const {return ((AliITSresponseSSD*)fResponse)->ADCToKeV(adc);}

  virtual Double_t GetCouplingPR() const {return ((AliITSresponseSSD*)fResponse)->GetCouplingPR();}
  virtual Double_t GetCouplingPL() const {return ((AliITSresponseSSD*)fResponse)->GetCouplingPL();}
  virtual Double_t GetCouplingNR() const {return ((AliITSresponseSSD*)fResponse)->GetCouplingNR();}
  virtual Double_t GetCouplingNL() const {return ((AliITSresponseSSD*)fResponse)->GetCouplingNL();}
  virtual void SetCouplings(Double_t pr, Double_t pl, Double_t nr, Double_t nl) 
    { ((AliITSresponseSSD*)fResponse)->SetCouplings(pr,pl,nr,nl);}

  virtual Int_t GetZSThreshold() const {return ((AliITSresponseSSD*)fResponse)->GetZSThreshold();}
  virtual void SetZSThreshold(Int_t zsth) 
    { ((AliITSresponseSSD*)fResponse)->SetZSThreshold(zsth);}

 void SetModule(Int_t mod){fModule = mod;}

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

 private:
    AliITSCalibrationSSD(const AliITSCalibrationSSD &source); // copy constructor
    AliITSCalibrationSSD& operator=(const AliITSCalibrationSSD &source); // ass. op.

    ClassDef(AliITSCalibrationSSD,4) //Response class for SSD
      };
#endif
