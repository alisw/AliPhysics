#ifndef ALITRACKFITTERKALMAN_H
#define ALITRACKFITTERKALMAN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//
//                        Kalman-Filter-like fit 
//   to a straight-line crossing a set of arbitrarily oriented planes.
//           (See AliTrackFitterKalman.cxx for the details)
//
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include "AliTrackFitter.h"

class AliTrackFitterKalman : public AliTrackFitter {
public:
  AliTrackFitterKalman() : 
     AliTrackFitter(), 
     fMaxChi2(fgkMaxChi2) {}

  AliTrackFitterKalman(AliTrackPointArray *array, Bool_t owner = kTRUE);
  virtual ~AliTrackFitterKalman() {}

  void SetMaxChi2(Double_t chi2) {fMaxChi2=chi2;}

  void   SetSeed(const Double_t par[6], const Double_t cov[15]);
  Bool_t MakeSeed(const AliTrackPoint *p, const AliTrackPoint *p2);

  Bool_t GetPCA(const AliTrackPoint &p, AliTrackPoint &p2) const;

  Bool_t Begin(Int_t first, Int_t last);
  Bool_t AddPoint(const AliTrackPoint *p);
  Bool_t Update() {return kTRUE;}

private:
  AliTrackFitterKalman(const AliTrackFitterKalman &kalman);
  AliTrackFitterKalman &operator=(const AliTrackFitterKalman& kalman);

  Bool_t Propagate(const AliTrackPoint *p);
  Double_t GetPredictedChi2(const AliTrackPoint *p) const;
  Bool_t Update(const AliTrackPoint *p,Double_t chi2);


  static const Double_t fgkMaxChi2;  // Default maximal allowed chi2 

  Double_t fMaxChi2;                 // A point is added if chi2 < fMaxChi2 

  ClassDef(AliTrackFitterKalman,3)   // Kalman-Filter fit to a straight line

};

#endif
