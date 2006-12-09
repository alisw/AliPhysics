#ifndef ALITPCRECOPARAM_H
#define ALITPCRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with TPC reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TObject.h"

class AliTPCRecoParam : public TObject
{
 public: 
  AliTPCRecoParam();
  virtual ~AliTPCRecoParam();
  Double_t GetCtgRange() const     { return fCtgRange;}
  Double_t GetMaxSnpTracker() const{ return fMaxSnpTracker;}
  Double_t GetMaxSnpTrack() const  { return fMaxSnpTrack;}
  //
  Int_t    GetFirstBin() const     { return fFirstBin;}
  Int_t    GetLastBin() const      { return fLastBin;}
  void     SetTimeBinRange(Int_t first, Int_t last){ fFirstBin = first; fLastBin = last;}
  Bool_t   GetCalcPedestal()       const  { return fBCalcPedestal;}
  Bool_t   GetDoUnfold()           const  { return fBDoUnfold;}
  Float_t  GetDumpAmplitudeMin()   const  { return fDumpAmplitudeMin;}
  Float_t  GetMaxNoise()           const  { return fMaxNoise;}  
  Float_t  GetMinMaxCutAbs()       const  { return fMinMaxCutAbs; }
  Float_t  GetMinLeftRightCutAbs() const  { return fMinLeftRightCutAbs;}  // minimal amplitude left right - PRF
  Float_t  GetMinUpDownCutAbs()    const  { return fMinUpDownCutAbs;}  // minimal amplitude up-down - TRF 
  Float_t  GetMinMaxCutSigma()       const  { return fMinMaxCutSigma; }
  Float_t  GetMinLeftRightCutSigma() const  { return fMinLeftRightCutSigma;}  // minimal amplitude left right - PRF
  Float_t  GetMinUpDownCutSigma()    const  { return fMinUpDownCutSigma;}  // minimal amplitude up-down - TRF 

  //
  Bool_t   GetDoKinks() const      { return fBKinkFinder;}
  Float_t  GetMaxC()    const      { return fMaxC;}
  Bool_t   GetSpecialSeeding() const { return fBSpecialSeeding;}
  Bool_t   GetBYMirror() const { return fBYMirror;}
  static   AliTPCRecoParam *GetLowFluxParam();        // make reco parameters for low  flux env.
  static   AliTPCRecoParam *GetHighFluxParam();       // make reco parameters for high flux env. 
  static   AliTPCRecoParam *GetLaserTestParam(Bool_t bPedestal);  // special setting for laser 
  static   AliTPCRecoParam *GetCosmicTestParam(Bool_t bPedestal); // special setting for cosmic  
  //
 protected:
  Double_t fCtgRange;        // +-fCtgRange is the ctg(Theta) window used for clusterization and tracking (MI) 
  Double_t fMaxSnpTracker;   // max sin of local angle  - for TPC tracker
  Double_t fMaxSnpTrack;     // max sin of local angle  - for track 
  Bool_t   fBYMirror;        // mirror of the y - pad coordinate 
  //
  //   clusterer parameters
  //
  Int_t    fFirstBin;        // first time bin used by cluster finder
  Int_t    fLastBin;         // last time bin  used by cluster finder 
  Bool_t   fBCalcPedestal;   // calculate Pedestal
  Bool_t   fBDoUnfold;       // do unfolding of clusters
  Float_t  fDumpAmplitudeMin; // minimal amplitude of signal to be dumped 
  Float_t  fMaxNoise;        // maximal noise sigma on pad to be used in cluster finder
  Float_t  fMinMaxCutAbs;    // minimal amplitude at cluster maxima
  Float_t  fMinLeftRightCutAbs;  // minimal amplitude left right - PRF
  Float_t  fMinUpDownCutAbs;  // minimal amplitude up-down - TRF 
  Float_t  fMinMaxCutSigma;    // minimal amplitude at cluster maxima
  Float_t  fMinLeftRightCutSigma;  // minimal amplitude left right - PRF
  Float_t  fMinUpDownCutSigma;  // minimal amplitude up-down - TRF 
  //
  //
  Float_t  fMaxC;            // maximal curvature for tracking
  Bool_t   fBSpecialSeeding; // special seeding with big inclination angles allowed (for Cosmic and laser)
  Bool_t   fBKinkFinder;     // do kink finder reconstruction
  ClassDef(AliTPCRecoParam, 1)
};


#endif
