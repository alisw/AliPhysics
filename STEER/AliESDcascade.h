#ifndef ALIESDCASCADE_H
#define ALIESDCASCADE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                        ESD Cascade Vertex Class
//               Implementation of the cascade vertex class
//    Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TPDGCode.h>
#include "AliESDv0.h"

class AliExternalTrackParam;

#define kXiMinus       3312
#define kXiPlusBar    -3312
#define kOmegaMinus    3334
#define kOmegaPlusBar -3334

class AliESDcascade : public AliESDv0 {

public:
  AliESDcascade();
  AliESDcascade(const AliESDcascade&);
  AliESDcascade(const AliESDv0 &v0,
                const AliExternalTrackParam &t, Int_t i);
  ~AliESDcascade();

  Double_t ChangeMassHypothesis(Double_t &v0q, Int_t code=kXiMinus); 

  Int_t    GetPdgCodeXi() const {return fPdgCodeXi;}
  Double_t GetEffMassXi() const {return fEffMassXi;}
  Double_t GetChi2Xi()  const {return fChi2Xi;}
  void     GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetXYZcascade(Double_t &x, Double_t &y, Double_t &z) const;
  Double_t GetDcascade(Double_t x0=0.,Double_t y0=0.,Double_t z0=0.) const;

  void     GetBPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
     px=fBachMom[0]; py=fBachMom[1]; pz=fBachMom[2];
  }

  Int_t    GetBindex() const {return fBachIdx;}
  void     SetDcaXiDaughters(Double_t rDcaXiDaughters=0.);
  Double_t GetDcaXiDaughters() const {return fDcaXiDaughters;}
  Double_t GetCascadeCosineOfPointingAngle(Double_t&, Double_t&, Double_t&) const;

  void GetPosCovXi(Double_t cov[6]) const;

protected: 
  Int_t      fPdgCodeXi;      // reconstructed cascade type (PDG code)
  Double32_t fEffMassXi;      // reconstructed cascade effective mass
  Double32_t fChi2Xi;         // chi2 value
  Double32_t fDcaXiDaughters; // dca between Xi's daughters
  Double32_t fPosXi[3];       // cascade vertex position (global)
  Double32_t fPosCovXi[6];    // covariance matrix of the vertex position

  Int_t      fBachIdx;        // label of the bachelor track
  Double32_t fBachMom[3];     // bachelor momentum (global)
  Double32_t fBachMomCov[6];  // covariance matrix of the bachelor momentum.

private:
  AliESDcascade& operator=(const AliESDcascade&);

  ClassDef(AliESDcascade,4) // reconstructed cascade vertex
};

inline
void AliESDcascade::SetDcaXiDaughters(Double_t rDcaXiDaughters){
  fDcaXiDaughters=rDcaXiDaughters;
}

#endif
