#ifndef ALICASCADEVERTEX_H
#define ALICASCADEVERTEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          Cascade Vertex Class
//
//    Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr
//-------------------------------------------------------------------------

#include <TObject.h>
#include "AliPDG.h"

class AliITStrackV2;
class AliV0vertex;

#define kXiMinus       3312
#define kXiPlusBar    -3312
#define kOmegaMinus    3334
#define kOmegaPlusBar -3334

class AliCascadeVertex : public TObject {
public:
  AliCascadeVertex();
  AliCascadeVertex(const AliV0vertex &vtx, const AliITStrackV2 &trk);

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
  Int_t GetNlabel() const {return fV0lab[0];}
  void GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
     px=fV0mom[1][0]; py=fV0mom[1][1]; pz=fV0mom[1][2];
  }
  Int_t GetPlabel() const {return fV0lab[1];}
  void GetBPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
     px=fBachMom[0]; py=fBachMom[1]; pz=fBachMom[2];
  }
  Int_t GetBlabel() const {return fBachLab;}

private: 
  Int_t fPdgCode;           // reconstructed cascade type (PDG code)
  Double_t fEffMass;        // reconstructed cascade effective mass
  Double_t fChi2;           // chi2 value
  Double_t fPos[3];         // cascade vertex position (global)
  Double_t fPosCov[6];      // covariance matrix of the vertex position

  Int_t fV0lab[2];          // labels of the V0 daughter tracks
  Double_t fV0mom[2][3];    // V0 daughters' momenta (global)
  Double_t fV0momCov[6];    // covariance matrix of the V0 momentum.

  Int_t fBachLab;           // label of the bachelor track
  Double_t fBachMom[3];      // bachelor momentum (global)
  Double_t fBachMomCov[6];   // covariance matrix of the bachelor momentum.

  ClassDef(AliCascadeVertex,1)   // reconstructed cascade vertex
};

#endif


