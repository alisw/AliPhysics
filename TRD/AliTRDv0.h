#ifndef TRDv0_H
#define TRDv0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 0    //
////////////////////////////////////////////////////////
 
#include "AliTRD.h"
 
class AliTRDv0 : public AliTRD {

public:
  AliTRDv0() {}
  AliTRDv0(const char *name, const char *title);
  virtual        ~AliTRDv0() {}
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual Int_t   IsVersion() const           { return 0;           };
  virtual void    Hits2Clusters();
  virtual void    SetHits(Int_t ihit = 1)     { fHitsOn = ihit;     };
  virtual void    StepManager();
  virtual void    Init();
  
  virtual void    SetRphiSigma(Float_t sigma) { fRphiSigma = sigma; };
  virtual void    SetRphiDist(Float_t dist)   { fRphiDist  = dist;  };

  virtual Float_t GetRphiSigma()              { return fRphiSigma;  };
  virtual Float_t GetRphiDist()               { return fRphiDist;   };

protected:
  Int_t        fIdSens;     // Sensitive volume identifier

  Int_t        fIdChamber1; // Driftchamber volume identifier
  Int_t        fIdChamber2; // 
  Int_t        fIdChamber3; // 

  Int_t        fHitsOn;     // Used to switch hits on

  Float_t      fRphiSigma;  // Gaussian position smearing in rphi-direction
  Float_t      fRphiDist;   // Maximum distnace for non-overlapping cluster

  ClassDef(AliTRDv0,1)      // Transition Radiation Detector version 0 (fast simulator)

};

#endif
