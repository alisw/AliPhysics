#ifndef ALIGEVSIMPARTICLE_H
#define ALIGEVSIMPARTICLE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////////
//
// AliGeVSimParticle is a helper class for GeVSim (AliGenGeVSim) event generator.
// An object of this class represents one particle type and contain 
// information about particle type thermal parameters.
//
//////////////////////////////////////////////////////////////////////////////
//
// For examples, parameters and testing macros refer to:
// http:/home.cern.ch/radomski
// 
// for more detailed description refer to ALICE NOTE
// "GeVSim Monte-Carlo Event Generator"
// S.Radosmki, P. Foka.
//  
// Author:
// Sylwester Radomski,
// GSI, March 2002
//  
// S.Radomski@gsi.de
//
//////////////////////////////////////////////////////////////////////////////
//
// Updated and revised: September 2002, S. Radomski, GSI
//
////////////////////////////////////////////////////////////////////////////////


#include "TObject.h"

class AliGeVSimParticle : public TObject {

 public:
  
  ////////////////////////////////////////////////////////////////////////////
  
  AliGeVSimParticle() {}
  AliGeVSimParticle(Int_t pdg, Int_t model, Float_t multiplicity); 
  AliGeVSimParticle(Int_t pdg, Int_t model, Float_t multiplicity, 
		    Float_t T, Float_t dY = 1., Float_t exp=0.);
  
  ~AliGeVSimParticle() {}
  
  ////////////////////////////////////////////////////////////////////////////
  
  Int_t GetPdgCode() const {return fPDG;}
  Int_t GetModel() const {return fModel;}

  Float_t GetTemperature() const {return fT;}
  Float_t GetSigmaY() const {return fSigmaY;}
  Float_t GetExpansionVelocity() const {return fExpansion;}

  void SetModel(Int_t model);
  void SetTemperature(Float_t T) {fT = T;}
  void SetSigmaY(Float_t sigma) {fSigmaY = sigma;}
  void SetExpansionVelocity(Float_t vel) {fExpansion = vel;}
  

  // Multiplicity

  void    SetMultiplicity(Float_t mult);
  Float_t GetMultiplicity() const {return fN;}

  void   SetMultTotal(Bool_t isTotal = kTRUE);

  Bool_t IsMultTotal() const {return fMultTotal;}
  Bool_t IsMultForced() const {return fIsSetMult;}
  
  // Flow
  
  void SetDirectedSimple(Float_t v1);
  void SetEllipticSimple(Float_t v2);

  void SetDirectedParam(Float_t v11, Float_t v12=0, Float_t v13=1, Float_t v14=0);
  void SetEllipticParam(Float_t v21, Float_t pTmax, Float_t v22=0.);
  void SetEllipticOld(Float_t v21, Float_t v22, Float_t v23);

  Bool_t IsFlowSimple() const;

  Float_t GetDirectedFlow(Float_t pt, Float_t y);
  Float_t GetEllipticFlow(Float_t pt, Float_t y);
  
  ////////////////////////////////////////////////////////////////////////////
  
 private:
  
  Int_t fPDG;            // Particle type code
  Int_t fModel;          // Transverse momentum model

  Float_t fN;            // Multiplicity (subject to scalling)
  Bool_t  fMultTotal;    // multiplicity mode: Total or dN/dY
  Bool_t  fIsSetMult;    // force multiplicity mode or use from AliGenGeVSim

  Float_t fT;            // Slope Parameter (subject to scalling)
  Float_t fSigmaY;       // Rapidity Width
  Float_t fExpansion;    // Expansion Velocity in c units (subject to scalling)
  
  Float_t fV1[4];        // Directed Flow coefficient parameters
  Float_t fV2[3];        // Elliptic Flow coefficient parameters
  
  Bool_t fIsDirectedSimple;  // indicate use constant value for directed (v1) 
  Bool_t fIsEllipticSimple;  // indicate use constant value for elliptic (v2)
  Bool_t fIsEllipticOld;     // linear / quadratical pT parametrisation

  ClassDef(AliGeVSimParticle, 3)
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
