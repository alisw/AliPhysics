#ifndef ALIV0VERTEX_H
#define ALIV0VERTEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          V0 Vertex Class
//
//    Origin: Iouri Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TPDGCode.h>

class AliITStrackV2;

class AliV0vertex : public TObject {
public:
  AliV0vertex();
  AliV0vertex(const AliITStrackV2 &neg, const AliITStrackV2 &pos);

  Double_t ChangeMassHypothesis(Int_t code=kK0Short); 

  Int_t GetPdgCode() const {return fPdgCode;}
  Double_t GetEffMass() const {return fEffMass;}
  Double_t GetChi2() const {return fChi2;}
  void GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void GetNPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void GetXYZ(Double_t &x, Double_t &y, Double_t &z) const;
  Double_t GetD(Double_t x0=0.,Double_t y0=0.,Double_t z0=0.) const;
  Int_t GetNlabel() const {return fNlab;}
  Int_t GetPlabel() const {return fPlab;}

private: 
  Int_t fPdgCode;           // reconstructed V0's type (PDG code)
  Double_t fEffMass;        // reconstructed V0's effective mass
  Double_t fChi2;           // V0's chi2 value
  Double_t fPos[3];         // V0's position (global)
  Double_t fPosCov[6];      // covariance matrix of the vertex position

  Int_t fNlab;              // label of the negative daughter
  Double_t fNmom[3];        // momentum of the negative daughter (global)
  Double_t fNmomCov[6];     // covariance matrix of the negative daughter mom.

  Int_t fPlab;              // label of the positive daughter
  Double_t fPmom[3];        // momentum of the positive daughter (global)
  Double_t fPmomCov[6];     // covariance matrix of the positive daughter mom.

  ClassDef(AliV0vertex,1)   // reconstructed V0 vertex
};

inline 
void AliV0vertex::GetNPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
px=fNmom[0]; py=fNmom[1]; pz=fNmom[2];
}

inline 
void AliV0vertex::GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
px=fPmom[0]; py=fPmom[1]; pz=fPmom[2];
}

#endif


