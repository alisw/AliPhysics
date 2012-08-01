#ifndef ALIITSRESPONSESPD_H
#define ALIITSRESPONSESPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliITSresponse.h"
///////////////////////////////////////////
//                                       //
// ITS response class for Pixels         //
///////////////////////////////////////////
  
class AliITSresponsePixUpg :  public AliITSresponse {
 public:
  AliITSresponsePixUpg(); // default constructor
  virtual ~AliITSresponsePixUpg() {;} // destructror
  //
  virtual  void   SetSigmaDiffusionAsymmetry(Double_t ecc)        {fEccDiff=ecc;}   
  virtual  void   GetSigmaDiffusionAsymmetry(Double_t &ecc) const {ecc=fEccDiff;}
  //
 protected:
  //
  static const Float_t fgkDiffCoeffDefault;  //default for fDiffCoeff
  static const TString fgkCouplingOptDefault;  // type of pixel Coupling (old or new)
  static const Float_t fgkEccentricityDiffDefault;//default for fCouplRow 
  
    TString fCouplOpt;        // Coupling Option
    Float_t fEccDiff;         // Eccentricity (i.e. asymmetry parameter) in the Gaussian Diffusion


    ClassDef(AliITSresponsePixUpg,1) // pixel upgrade base response class
};

#endif
