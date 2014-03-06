#ifndef FITV0_H
#define FITV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  simple geometry  classes for set: FIT  
// Alla.Maevskaya@cern.ch 
////////////////////////////////////////////////
 
#include "AliFIT.h"
#include "TGraph.h" 
class AliFITv0 : public AliFIT {
  
public:

  enum constants {kAir=1, kVac=3,kGlass=6, kSensAir=22};

 
  AliFITv0();
  AliFITv0(const char *name, const char *title);
  AliFITv0(const AliFITv0& o):AliFIT(),
    fIdSens1(0) {((AliFITv0 &) o).Copy(*this);}
  
  AliFITv0& operator=(const AliFITv0&) { return *this; }
  virtual       ~AliFITv0();
  virtual void   CreateGeometry();
  virtual void   AddAlignableVolumes() const;
  virtual void   CreateMaterials() ;
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();

protected:
  Int_t fIdSens1; // Sensetive volume  in T0
 
  ClassDef(AliFITv0,1)  //Class for FIT version 1
};


#endif


