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

class AliESDtrack;
class AliESDv0;

#define kXiMinus       3312
#define kXiPlusBar    -3312
#define kOmegaMinus    3334
#define kOmegaPlusBar -3334

class AliESDcascade : public TObject {
public:
  AliESDcascade();

  Double_t ChangeMassHypothesis(Double_t &v0q, Int_t code=kXiMinus); 

  Int_t GetPdgCode() const {return fPdgCode;}
  Double_t GetEffMass() const {return fEffMass;}
  Double_t GetChi2() const {return fChi2;}
  void GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void GetXYZ(Double_t &x, Double_t &y, Double_t &z) const;
  Double_t GetD(Double_t x0=0.,Double_t y0=0.,Double_t z0=0.) const;

  void GetNPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
     px=fV0mom[0][0]; py=fV0mom[0][1]; pz=fV0mom[0][2];
  }
  Int_t GetNindex() const {return fV0idx[0];}
  void GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
     px=fV0mom[1][0]; py=fV0mom[1][1]; pz=fV0mom[1][2];
  }
  Int_t GetPindex() const {return fV0idx[1];}
  void GetBPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
     px=fBachMom[0]; py=fBachMom[1]; pz=fBachMom[2];
  }
  Int_t GetBindex() const {return fBachIdx;}

protected: 
  Int_t fPdgCode;           // reconstructed cascade type (PDG code)
  Double_t fEffMass;        // reconstructed cascade effective mass
  Double_t fChi2;           // chi2 value
  Double_t fPos[3];         // cascade vertex position (global)
  Double_t fPosCov[6];      // covariance matrix of the vertex position

  Int_t fV0idx[2];          // indeices of the V0 daughter tracks
  Double_t fV0mom[2][3];    // V0 daughters' momenta (global)
  Double_t fV0momCov[6];    // covariance matrix of the V0 momentum.

  Int_t fBachIdx;           // label of the bachelor track
  Double_t fBachMom[3];     // bachelor momentum (global)
  Double_t fBachMomCov[6];  // covariance matrix of the bachelor momentum.

  ClassDef(AliESDcascade,1) // reconstructed cascade vertex
};

#endif


