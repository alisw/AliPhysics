#ifndef ALIPMDRECONSTRUCTOR_H
#define ALIPMDRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliReconstructor.h"
#include "AliPMDRecoParam.h"

class AliPMDReconstructor: public AliReconstructor {
public:
 AliPMDReconstructor() : AliReconstructor() {}

  virtual void   Reconstruct(AliRawReader* rawReader,
			     TTree* clustersTree) const;
  virtual void   Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  virtual void   FillESD(AliRawReader* /*rawReader*/, TTree* clustersTree, 
			 AliESDEvent* esd) const;

  virtual void   FillESD(TTree* /*digitsTree*/, TTree* clustersTree, 
			 AliESDEvent* esd) const;

  static const AliPMDRecoParam* GetRecoParam() { return dynamic_cast<const AliPMDRecoParam*>(AliReconstructor::GetRecoParam(10)); } // getting RecoParam obj

private:

  ClassDef(AliPMDReconstructor, 6)   // class for the PMD reconstruction
};

#endif
