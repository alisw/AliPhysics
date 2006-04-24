#ifndef ALITRDSEED_H
#define ALITRDSEED_H   

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

#include "TObject.h" 

class AliTRDcluster;

class AliTRDseed : public TObject {

  friend class AliTRDtracker;

 public:

  AliTRDseed(); 
  ~AliTRDseed() {};                 

  static void    EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh);
  static Float_t FitRiemanTilt(AliTRDseed *seed, Bool_t error);
  void           UseClusters(); // use clusters
  void           Update();      // update information - without tilt correction
  void           CookLabels();  // cook label
  void           UpdateUsed();
  void           Reset();       // reset seed
  Bool_t         IsOK() const { return fN2 > 8;}

 private:

  Float_t        fTilt;         // tilting angle
  Float_t        fPadLength;    // pad length
  Float_t        fX0;           // x0 position
  Float_t        fX[25];        // !x position
  Float_t        fY[25];        // !y position
  Float_t        fZ[25];        // !z position
  Int_t          fIndexes[25];  // !indexes
  AliTRDcluster *fClusters[25]; // !clusters
  Bool_t         fUsable[25];   // !indication  - usable cluster
  Float_t        fYref[2];      // reference y
  Float_t        fZref[2];      // reference z
  Float_t        fYfit[2];      // y fit position +derivation
  Float_t        fYfitR[2];     // y fit position +derivation
  Float_t        fZfit[2];      // z fit position
  Float_t        fZfitR[2];     // z fit position
  Float_t        fSigmaY;       // "robust" sigma in Y - constant fit
  Float_t        fSigmaY2;      // "robust" sigma in Y - line fit
  Float_t        fMeanz;        // mean vaue of z
  Float_t        fZProb;        // max probbable z
  Int_t          fLabels[2];    // labels
  Int_t          fN;            // number of associated clusters
  Int_t          fN2;           // number of not crossed
  Int_t          fNUsed;        // number of used clusters
  Int_t          fFreq;         // freq
  Int_t          fNChange;      // change z counter
  Float_t        fMPads;        // mean number of pads per cluster
  // global
  //
  Float_t        fC;            // curvature
  Float_t        fCC;           // curvature with constrain
  Float_t        fChi2;         // global chi2
  Float_t        fChi2Z;        // global chi2

  ClassDef(AliTRDseed,1)        // Seed for a local TRD track

};
#endif 
