#ifndef FITV2_H
#define FITV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
// Full geomrtry  hits classes for detector: FIT    //
////////////////////////////////////////////////
 
#include "AliFIT.h"
#include "TGraph.h" 
class AliFITv2 : public AliFIT {
  
public:

  enum constants {kAir=1, kVac=3, kGlass=6, kOpGlass=16, kOpGlassCathode=19,kSensAir=22};

 
  AliFITv2();
  AliFITv2(const char *name, const char *title);
  AliFITv2(const AliFITv2& o):AliFIT(),
    fIdSens1(0),
    fPMTeff(0x0) {((AliFITv2 &) o).Copy(*this);}
  
  AliFITv2& operator=(const AliFITv2&) { return *this; }
  virtual       ~AliFITv2();
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
 
  ClassDef(AliFITv2,1)  //Class for FIT version 1
};


#endif


