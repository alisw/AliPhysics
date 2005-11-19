#ifndef ALIPMDRECONSTRUCTOR_H
#define ALIPMDRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliReconstructor.h"

class AliPMDReconstructor: public AliReconstructor {
public:
  virtual void   Init(AliRunLoader* runLoader);
  virtual void   Reconstruct(AliRunLoader* runLoader) const;
  virtual void   Reconstruct(AliRunLoader* runLoader,
			     AliRawReader *rawReader) const;
  virtual void   Reconstruct(AliRawReader* rawReader,
			     TTree* clustersTree) const;
  virtual void   Reconstruct(TTree* digitsTree, TTree* clustersTree) const {
    AliReconstructor::Reconstruct(digitsTree,clustersTree);
  }
  virtual Bool_t HasLocalReconstruction() const { return kTRUE; }

  //virtual void   FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  virtual void   FillESD(AliRawReader* /*rawReader*/, TTree* clustersTree, 
			 AliESD* esd) const;
  virtual void   FillESD(TTree* digitsTree, TTree* clustersTree, 
			 AliESD* esd) const {
    AliReconstructor::FillESD(digitsTree,clustersTree,esd);
  }
  virtual void   FillESD(AliRunLoader* runLoader, AliESD* esd) const {
    AliReconstructor::FillESD(runLoader,esd);
  }
  virtual void   FillESD(AliRunLoader* runLoader, 
			 AliRawReader* rawReader, AliESD* esd) const {
    AliReconstructor::FillESD(runLoader,rawReader,esd);
  }
 
private:

  ClassDef(AliPMDReconstructor, 3)   // class for the PMD reconstruction
};

#endif
