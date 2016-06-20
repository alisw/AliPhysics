#ifndef FITV4_H
#define FITV4_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
// Full geomrtry  hits classes for detector: FIT    //
////////////////////////////////////////////////
 
#include "AliFIT.h"
#include "TGraph.h" 
class AliFITv4 : public AliFIT {
  
public:

  enum constants {kAir=1, kVac=3, kGlass=6, kOpGlass=16, kOpGlassCathode=19,kSensAir=22};

 
  AliFITv4();
  AliFITv4(const char *name, const char *title);
  AliFITv4(const AliFITv4& o):AliFIT(),
    fIdSens1(0),fIdSens2(0),
    fPMTeff(0x0) {((AliFITv4 &) o).Copy(*this);}
  
  AliFITv4& operator=(const AliFITv4&) { return *this; }
  virtual       ~AliFITv4();
  virtual void   CreateGeometry();
  virtual void   DefineOpticalProperties();
  virtual void   AddAlignableVolumes() const;
  virtual void   CreateMaterials() ;
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  Bool_t RegisterPhotoE(Double_t energy);
  virtual void   StepManager();
  void SetPMTeff();

protected:
  Int_t fIdSens1; // Sensetive volume  in T0
  Int_t fIdSens2; // Sensetive volume  in T0
  TGraph *fPMTeff; //pmt registration effeicincy
 
  ClassDef(AliFITv4,1)  //Class for FIT version 2
};


#endif


