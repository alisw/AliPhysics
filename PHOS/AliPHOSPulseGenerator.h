#ifndef ALIPHOSPULSEGENERATOR_H
#define ALIPHOSPULSEGENERATOR_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */

// The class which simulates the pulse shape from the PHOS FEE shaper,
// make sampled amplitudes, digitize them.
// The shape is described by the function RawResponseFunction
// The input parameters for the shape function (time and aplitude) are passed
// to the class via constructor.
// Other parameters related to the shaper are hard-coded in this class

#include <Rtypes.h>

class AliPHOSPulseGenerator : public TObject
{
public:
  AliPHOSPulseGenerator(Double_t a=0, Double_t t0=0);
  AliPHOSPulseGenerator(const AliPHOSPulseGenerator & pulse);
  virtual  ~AliPHOSPulseGenerator();

  void      AddBaseline(Double_t baselineLevel);
  void      AddNoise   (Double_t sigma);
  void      AddNoise   (Double_t *sigma, Double_t cutoff);
  void      AddPretriggerSamples(Int_t nPresamples);
  void      GetSamples(Int_t *adcHG, Int_t *adcLG) const;
  Bool_t    MakeSamples();
  void      Digitize();
  Bool_t    GetDigitize() {return fDigitize;}
  void      SetDigitise (Bool_t flag) {fDigitize  = flag;}
  void      SetAmplitude(Double_t  a) {fAmplitude = a   ; Reset();}
  void      SetTZero    (Double_t t0) {fTZero     = t0  ; Reset();}
  void      Reset();

  // Raw Read Out
  Double_t        GetRawFormatCapa()        const { return fgCapa ; }
  static Double_t GetRawFormatHighCharge()        { return fgHighCharge ; }
  static Double_t GetRawFormatHighGain()          { return fgHighGain ; }
  static Double_t GetRawFormatHighLowGainFactor() { return fgHighLowGainFactor ; }
  static Double_t GetRawFormatLowCharge()         { return ( fgHighCharge *  fgHighLowGainFactor ) ; }
  static Double_t GetRawFormatLowGain()           { return ( fgHighGain / fgHighLowGainFactor ) ; }
  Int_t           GetRawFormatOrder()       const { return fgOrder ; }
  static Int_t    GetRawFormatTimeBins()          { return fkTimeBins ; }
  Double_t        GetRawFormatTimeMax()     const { return fgTimeMax ; }
  Double_t        GetRawFormatTimePeak()    const { return fgTimePeak ; }
  Double_t        GetRawFormatTimeTrigger() const { return fgTimeTrigger ; }
  static Double_t RawResponseFunctionMax(Double_t charge, Double_t gain) ;
  static Double_t RawResponseFunction   (Double_t *x, Double_t *par) ;

  virtual void Print(Option_t*) const;
  virtual void Draw (Option_t*);

  AliPHOSPulseGenerator& operator = (const AliPHOSPulseGenerator &) {
    Fatal("operator =", "not implemented") ;
    return *this;
  }

private:
  static Double_t fgCapa ;              // capacitor of the preamplifier
  static Double_t fgHighCharge ;        // high charge (to convert energy to charge)
  static Double_t fgHighGain ;          // high gain
  static Double_t fgHighLowGainFactor ; // high to low gain factor
  static Int_t    fgOrder ;             // order of the gamma function
  static const Int_t fkTimeBins = 64 ;  // number of sampling bins
  static Double_t fgTimeMax ;           // maximum sampled time
  static Double_t fgTimePeak ;          // peaking time
  static Double_t fgTimeTrigger ;       // time of the trigger for the RO signal 
  
private:
  Double_t  fAmplitude;    // signal amplitude in GeV
  Double_t  fTZero;        // signal start time in ns
  Double_t *fDataHG;       // samples array for high gain
  Double_t *fDataLG;       // samples array for low  gain
  Bool_t    fDigitize;     // true is samples should be rounded to integers
  
  ClassDef(AliPHOSPulseGenerator,1)

};

#endif

