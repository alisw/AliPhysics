#ifndef ALITRDCLUSTERIZERV0_H
#define ALITRDCLUSTERIZERV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRD.h"
#include "AliTRDclusterizer.h"

///////////////////////////////////////////////////////
//  Finds and handles cluster (fast simulation)      //
///////////////////////////////////////////////////////

class AliTRDclusterizerV0 : public AliTRDclusterizer {

 public:

  AliTRDclusterizerV0();
  AliTRDclusterizerV0(const Text_t* name, const Text_t* title);
  virtual ~AliTRDclusterizerV0();

  virtual void    Init();
  virtual Bool_t  MakeCluster();
  
  virtual void    SetRphiSigma(Float_t sigma) { fRphiSigma = sigma; };
  virtual void    SetRphiDist(Float_t dist)   { fRphiDist  = dist;  };

  virtual Float_t GetRphiSigma()              { return fRphiSigma;  };
  virtual Float_t GetRphiDist()               { return fRphiDist;   };

 protected:

  Float_t      fRphiSigma;           // Gaussian position smearing in rphi-direction
  Float_t      fRphiDist;            // Maximum distance for non-overlapping cluster

  ClassDef(AliTRDclusterizerV0,1)    // TRD-Cluster manager, fast simulator

};

#endif
