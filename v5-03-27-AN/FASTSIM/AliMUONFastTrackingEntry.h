#ifndef ALIMUONFASTTRACKINGENTRY_H
#define ALIMUONFASTTRACKINGENTRY_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include <TClassTable.h>

static const Int_t kSplitP = 5; 
static const Int_t kSplitTheta = 3; 

class AliMUONFastTrackingEntry {
 public:
  AliMUONFastTrackingEntry();
  virtual ~AliMUONFastTrackingEntry(){;}
  Float_t GetP()const {return fP;}
  Float_t GetTheta()const {return fTheta;}
  Float_t GetPhi()const {return fPhi;}
  Float_t GetMeanp()const {return fMeanp;}
  Float_t GetMeantheta()const {return fMeantheta;}
  Float_t GetMeanphi()const {return fMeanphi;}
  Float_t GetSigmap()const {return fSigmap;}
  Float_t GetSigmatheta()const {return fSigmatheta;}
  Float_t GetSigmaphi()const {return fSigmaphi;}
  Float_t GetSigma1p()const {return fSigma1p;}
  Float_t GetChi2p()const {return fChi2p;}
  Float_t GetChi2theta()const {return fChi2theta;}
  Float_t GetChi2phi()const {return fChi2phi;}
  Float_t GetAcc(Int_t i, Int_t j)const {return fAcc[i][j];}
  Float_t GetEff(Int_t i, Int_t j) const {return fEff[i][j];}
  Float_t GetNormG2()const {return fNormG2;}
  Float_t GetMeanG2()const {return fMeanG2;}
  Float_t GetSigmaG2()const {return fSigmaG2;}

  void SetP(Float_t p){fP = p;}
  void SetTheta(Float_t theta){fTheta = theta;}
  void SetPhi(Float_t phi){fPhi = phi;}
  void SetMeanp(Float_t meanp){fMeanp = meanp;}
  void SetMeantheta(Float_t meantheta){fMeantheta = meantheta;}
  void SetMeanphi(Float_t meanphi){fMeanphi = meanphi;}
  void SetSigmap(Float_t sigmap){fSigmap = sigmap;}
  void SetSigmatheta(Float_t sigmatheta){fSigmatheta = sigmatheta;}
  void SetSigmaphi(Float_t sigmaphi){fSigmaphi = sigmaphi;}
  void SetSigma1p(Float_t sigma1p){fSigma1p = sigma1p;}
  void SetChi2p(Float_t chi2p){fChi2p = chi2p;}
  void SetChi2theta(Float_t chi2theta){fChi2theta = chi2theta;}
  void SetChi2phi(Float_t chi2phi){fChi2phi = chi2phi;}
  void SetAcc(Int_t i, Int_t j, Float_t acc) {fAcc[i][j] = acc;}
  void SetEff(Int_t i, Int_t j, Float_t eff) {fEff[i][j] = eff;}
  void SetNormG2(Float_t normG2){fNormG2 = normG2;}
  void SetMeanG2(Float_t meanG2){fMeanG2 = meanG2;}
  void SetSigmaG2(Float_t sigmaG2){fSigmaG2 = sigmaG2;}

 protected:
  Float_t fP;              // momentum 
  Float_t fTheta;          // polar angle 
  Float_t fPhi;            // azimuth 
  Float_t fMeanp;          // mean value of p distribution in current LUT cell
  Float_t fMeantheta;      // mean value of theta distr. in current LUT cell
  Float_t fMeanphi;        // mean value of phi distr. in current LUT cell
  Float_t fSigmap;         // sigma of p distr. in current LUT cell
  Float_t fSigmatheta;     // sigma of theta distr. in current LUT cell
  Float_t fSigmaphi;       // sigma of phi distr. in current LUT cell
  Float_t fSigma1p;        // param. for asymmetry in p distribution
  Float_t fChi2p;          // chi2 for p 
  Float_t fChi2theta;      // chi2 for theta
  Float_t fChi2phi;        // chi2 for phi
  Float_t fAcc[5][3];      // acceptance (subdivided in narrower cells in p and theta for low momenta) 
  Float_t fEff[5][3];      // efficiency (subdivided in narrower cells in p and theta for low momenta) 
  Float_t fNormG2;         // params for momentum gaussian smearing due to BKG
  Float_t fMeanG2;         // params for momentum gaussian smearing due to BKG
  Float_t fSigmaG2;        // params for momentum gaussian smearing due to BKG
  ClassDef(AliMUONFastTrackingEntry,1)       
};


#endif
