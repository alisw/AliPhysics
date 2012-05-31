/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALITPCEVENTINFO_H
#define ALITPCEVENTINFO_H

#include <TObject.h>

class AliTPCeventInfo:public TObject {
public:
  AliTPCeventInfo();
  virtual ~AliTPCeventInfo();

  Int_t GetNSamples() const {return fNSamples;}
  void SetNSamples(Int_t NSamples) {fNSamples=NSamples;};
  Int_t GetSamplingRate() const {return fSamplingRate;}
  void SetSamplingRate(Int_t SamplingRate) {fSamplingRate=SamplingRate;}
  Int_t GetSamplingPhase() const {return fSamplingPhase;}
  void SetSamplingPhase(Int_t SamplingPhase) {fSamplingPhase=SamplingPhase;}

  Bool_t GetBC() const {return fBC;}
  void SetBC(Bool_t BC) {fBC=BC;}
  Int_t GetBCMode1() const {return fBCMode1;}
  void SetBCMode1(Int_t BCMode1) {fBCMode1=BCMode1;}
  Int_t GetBCMode2() const {return fBCMode2;}
  void SetBCMode2(Int_t BCMode2) {fBCMode2=BCMode2;}
  Int_t GetBCPre() const {return fBCPre;}
  void SetBCPre(Int_t BCPre) {fBCPre=BCPre;}
  Int_t GetBCPost() const {return fBCPost;}
  void SetBCPost(Int_t BCPost) {fBCPost=BCPost;}

  Int_t GetNPreTrigger() const {return fNPreTrigger;}
  void SetNPreTrigger(Int_t NPreTrigger) {fNPreTrigger=NPreTrigger;}
  Bool_t GetTC() const {return fTC;}
  void SetTC(Bool_t TC) {fTC=TC;}

  Bool_t GetZS() const {return fZS;}
  void SetZS(Bool_t ZS) {fZS=ZS;}
  Int_t GetZSGlitchFilter() const {return fZSGlitchFilter;}
  void SetZSGlitchFilter(Int_t ZSGlitchFilter) {fZSGlitchFilter=ZSGlitchFilter;}
  Int_t GetZSPre() const {return fZSPre;}
  void SetZSPre(Int_t ZSPre) {fZSPre=ZSPre;}
  Int_t GetZSPost() const {return fZSPost;}
  void SetZSPost(Int_t ZSPost) {fZSPost=ZSPost;}
  Int_t GetZSOffset() const {return fZSOffset;}
  void SetZSOffset(Int_t ZSOffset) {fZSOffset=ZSOffset;}

private:
  Int_t fNSamples; // Number of samples per channel
  Int_t fSamplingRate; // ?
  Int_t fSamplingPhase; // ?

  Bool_t fBC; // Baseline correction enabled?
  Int_t fBCMode1; // 1st baseline correction mode
  Int_t fBCMode2; // 2nd baseline correction mode
  Int_t fBCPre; // Number of presamples excluded from 2nd baseline correction
  Int_t fBCPost; // Number of postsamples excluded from 2nd baseline correction
  Int_t fNPreTrigger; // Number of pre-trigger samples
  Bool_t fTC; // Tail cancelation filter?

  Bool_t fZS; // Zero suppresion enabled?
  Int_t fZSGlitchFilter; // Glitch filter configuration for zero suppression
  Int_t fZSPre; // Number of presamples excluded from zero suppression
  Int_t fZSPost; // Number of postsamples excluded from zero suppression
  Int_t fZSOffset; // Zero suppression offset

  ClassDef(AliTPCeventInfo,1)
};

#endif
