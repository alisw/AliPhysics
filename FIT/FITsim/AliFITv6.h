#ifndef FITV6_H
#define FITV6_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
// Full geomrtry  hits classes for detector: FIT    //
////////////////////////////////////////////////
 
#include "AliFIT.h"
#include "TGraph.h" 
#include <sstream>

class AliFITv6 : public AliFIT {
  
public:

  enum constants {kAir=1, kVac=3, kGlass=6, kOpGlass=16, kOpGlassCathode=19,kSensAir=22};

 
  AliFITv6();
  AliFITv6(const char *name, const char *title);
  AliFITv6(const AliFITv6& o):AliFIT(),
    fIdSens1(0),fIdSens2(0),
    fPMTeff(0x0) {((AliFITv6 &) o).Copy(*this);}
  
  AliFITv6& operator=(const AliFITv6&) { return *this; }
  virtual       ~AliFITv6();
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
  Int_t GetCellId(Int_t *vol);

  //V0+
  Int_t nSectors, nRings;
  Int_t fIdV0Plus[8][5];//Sensitive volumes [nSectors][nRings], if modified then update the construct in .cxx
  Int_t fCellId;//Scintillator cell number
  Int_t fSenseless;//Senseless for T0+

private: 
  // Optical properties reader: e-Energy, abs-AbsorptionLength[cm], n-refractive index
  Int_t ReadOptProperties(const string inputFilePath, Float_t **e,
    Double_t **de, Float_t **abs, Float_t **n, Int_t &kNbins) const;
  void FillOtherOptProperties(Float_t **efficAll, Float_t **rindexAir,
    Float_t **absorAir, Float_t **rindexCathodeNext, Float_t **absorbCathodeNext,
    Double_t **efficMet, Double_t **aReflMet, const Int_t kNbins) const;
  void DeleteOptPropertiesArr(Float_t **e, Double_t **de, Float_t **abs,
    Float_t **n, Float_t **efficAll, Float_t **rindexAir, Float_t **absorAir,
    Float_t **rindexCathodeNext, Float_t **absorbCathodeNext,
    Double_t **efficMet, Double_t **aReflMet) const; // should be called at the very end of the simulation to free the memory

  //V0+ parameters related to geometry
  Double_t fV0PlusR0, fV0PlusR1, fV0PlusR2, fV0PlusR3, fV0PlusR4, fV0PlusR5, fV0PlusR6;
  Double_t fV0PlusSciWd, fV0PlusFraWd, fV0PlusZposition;
  Float_t fV0PlusnMeters; 
  
  //V0+ parameters related to light production:
  Double_t fV0PlusLightYield;       // Lightyield in BC404 (from V0A)
  Double_t fV0PlusLightAttenuation; // LightAttenuation in fibers (from V0A)
  Double_t fV0PlusFibToPhot;        // Loss in Fibers - Photocathode Connection (from V0A)
 
  ClassDef(AliFITv6,1)  //Class for FIT version 6
};


#endif
