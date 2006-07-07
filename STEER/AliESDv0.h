#ifndef ALIESDV0_H
#define ALIESDV0_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          ESD V0 Vertex Class
//          This class is part of the Event Summary Data set of classes
//    Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TPDGCode.h>

class AliExternalTrackParam;

class AliESDv0 : public TObject {
public:
  AliESDv0();
  AliESDv0(const AliExternalTrackParam &t1, Int_t i1,
           const AliExternalTrackParam &t2, Int_t i2);

  AliESDv0(const AliESDv0&);
  virtual ~AliESDv0();
  AliESDv0& operator=(const AliESDv0&);

  Double_t ChangeMassHypothesis(Int_t code=kK0Short); 

  Int_t    GetPdgCode() const {return fPdgCode;}
  Double_t GetEffMass() const {return fEffMass;}
  Double_t GetChi2V0() const {return fChi2V0;}
  void     GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetNPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetXYZ(Double_t &x, Double_t &y, Double_t &z) const;
  Double_t GetD(Double_t x0=0.,Double_t y0=0.,Double_t z0=0.) const;
  Int_t    GetNindex() const {return fNidx;}
  Int_t    GetPindex() const {return fPidx;}
  void     SetESDindexes(Int_t ip, Int_t im){fNidx=ip;fPidx=im;}
  void     SetDcaV0Daughters(Double_t rDcaV0Daughters=0.);
  Double_t GetDcaV0Daughters() {return fDcaV0Daughters;}
  Double_t GetV0CosineOfPointingAngle(Double_t&, Double_t&, Double_t&) const;

protected: 
  Int_t    fPdgCode;        // reconstructed V0's type (PDG code)
  Double_t fEffMass;        // reconstructed V0's effective mass
  Double_t fDcaV0Daughters; // dca between V0's daughters
  Double_t fChi2V0;         // V0's chi2 value
  Double_t fPos[3];         // V0's position (global)
  Double_t fPosCov[6];      // covariance matrix of the vertex position

  Int_t    fNidx;           // index of the negative daughter
  Double_t fNmom[3];        // momentum of the negative daughter (global)
  Double_t fNmomCov[6];     // covariance matrix of the negative daughter mom.

  Int_t    fPidx;           // index of the positive daughter
  Double_t fPmom[3];        // momentum of the positive daughter (global)
  Double_t fPmomCov[6];     // covariance matrix of the positive daughter mom.

  ClassDef(AliESDv0,1)      // ESD V0 vertex
};

inline 
void AliESDv0::GetNPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
px=fNmom[0]; py=fNmom[1]; pz=fNmom[2];
}

inline 
void AliESDv0::GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
px=fPmom[0]; py=fPmom[1]; pz=fPmom[2];
}

inline
void AliESDv0::SetDcaV0Daughters(Double_t rDcaV0Daughters){
  fDcaV0Daughters=rDcaV0Daughters;
}

#endif
