#ifndef ALIPMDRECONSTRUCTOR_H
#define ALIPMDRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliReconstructor.h"

class AliPMDReconstructor: public AliReconstructor {
public:
  virtual void   Reconstruct(AliRawReader* rawReader,
			     TTree* clustersTree) const;
  virtual void   Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  virtual void   FillESD(AliRawReader* /*rawReader*/, TTree* clustersTree, 
			 AliESDEvent* esd) const;

  virtual void   FillESD(TTree* /*digitsTree*/, TTree* clustersTree, 
			 AliESDEvent* esd) const;

private:

  ClassDef(AliPMDReconstructor, 5)   // class for the PMD reconstruction
};

#endif
