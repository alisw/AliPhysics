#ifndef ALIT0RECOPARAM_H
#define ALIT0RECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with T0 reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliDetectorRecoParam.h"

class AliT0RecoParam : public AliDetectorRecoParam
{
 public: 
  AliT0RecoParam();
  AliT0RecoParam(const AliT0RecoParam &p); //copy constructor
  AliT0RecoParam& operator=(const AliT0RecoParam &p);
  virtual ~AliT0RecoParam();
  
  static   AliT0RecoParam *GetLowFluxParam();        // make reco parameters for low  flux env
  static   AliT0RecoParam *GetHighFluxParam();       // make reco parameters for high flux env 
  static   AliT0RecoParam *GetLaserTestParam();  // special setting for laser SetLaserTestParam 
  //for monitoring
  //  static   AliT0RecoParam *GetHistRange();  //  limit of monitoring histograms
  
  //old staff  
  Float_t GetRefAmp()  const  {return fRefAmp;}
  void    SetRefAmp(Float_t amp)   { fRefAmp = amp;}
  //now number of bad channel
  Int_t   GetRefPoint() const {return fRefPoint;}
  void    SetRefPoint(Int_t ref) {fRefPoint = ref;}
  
  //now low and high limit for multi-bunch recontruction
  Float_t   GetLow(Int_t numhist) const {return fLow[numhist];}
  //  Float_t   GetLow() {return *fLow;}
  void      SetLow(Int_t numhist, Float_t low) {fLow[numhist] = low;}
  
  Float_t   GetHigh(Int_t numhist) const  {return fHigh[numhist];}
  //  Float_t   GetHigh()  {return *fHigh;}
  void      SetHigh(Int_t numhist, Float_t high) {fHigh[numhist] = high;}
  Float_t   GetLatencyL1() const {return fLatencyL1;}
  void      SetLatencyL1(Float_t lat) {fLatencyL1 = lat;}
  Float_t   GetLatencyL1A() const {return fLatencyL1A;}
  void      SetLatencyL1A(Float_t lat) {fLatencyL1A = lat;}
  Float_t   GetLatencyL1C() const {return fLatencyL1C;}
  void      SetLatencyL1C(Float_t lat) {fLatencyL1C = lat;}
  Float_t   GetLatencyHPTDC() const {return fLatencyHPTDC;}
  void      SetLatencyHPTDC(Float_t lat) {fLatencyHPTDC = lat;}
  Float_t   GetVertexShift() const {return fVertexShift;}
  void      SetVertexShift(Float_t sh) {fVertexShift = sh;}
  
  //new staff
  Int_t  GetBadChannels(Int_t i) const {return fBadChannels[i];}
  void SetBadChannels(Int_t i, Int_t value) {fBadChannels[i] = value;}
  Float_t GetAmpLowThreshold() const {return fLow[200];}
  Float_t GetAmpHighThreshold() const {return fHigh[200];}
  
  void SetSatelliteThresholds(Float_t low, Float_t high) 
  {fSatelliteThresholds[0]=low;  fSatelliteThresholds[1]=high;}
  Float_t GetLowSatelliteThreshold() const {return fSatelliteThresholds[0];}
  Float_t GetHighSatelliteThreshold() const {return fSatelliteThresholds[1];}
  void SetEq (Int_t eq) {fEqualised = eq; };
  Int_t GetEq () const { return fEqualised;}

  void PrintParameters() const;
  
 protected:
  Float_t   fRefAmp;            // for slewing correcton
  Int_t     fRefPoint;          // #channel for RefPoint
  Float_t   fLow[500];          //low limit of monitoring histograms
  Float_t   fHigh[500];         //high limit of monitoring histograms
  Float_t   fLatencyL1;         //Latency L1
  Float_t   fLatencyL1A;         //Latency L1 for OrA
  Float_t   fLatencyL1C;         //Latency L1 for orC
  Float_t   fLatencyHPTDC;      //Latency HPTDC
  Float_t   fVertexShift;       // for slewing correcton
  Int_t  fBadChannels[24];   // bad channels map
  Float_t  fSatelliteThresholds[2];   // what we call satellite
  Int_t fEqualised;                   // do we write pure CFD or equalized 
  
  ClassDef(AliT0RecoParam, 6);
 
};
#endif
