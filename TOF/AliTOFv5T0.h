#ifndef ALITOFv5T0_H
#define ALITOFv5T0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________//
//                                                                         //
// Implementation version v5 of TOF Manager class                          //
// FULL COVERAGE VERSION + OPTION FOR PHOS HOLES                           //
//                                                                         //
// -- Authors: G. Cara Romeo, A. De Caro                                   //
//                                                                         //
//_________________________________________________________________________//

#include "AliTOF.h"
 
 
class AliTOFv5T0 : public AliTOF {

public:
  AliTOFv5T0();
  AliTOFv5T0(const char *name, const char *title);
  virtual ~AliTOFv5T0() {};
  virtual void   BuildGeometry();
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 7;}
  virtual void   AddAlignableVolumes() const;
  virtual void   TOFpc(Float_t xtof,  Float_t ytof, Float_t zlenA,
		       Float_t zlenB);
  virtual void   TOFpc(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t) {};
  virtual void   StepManager();
  virtual void   DrawModule() const;
  virtual void   DrawDetectorModules() const;
  virtual void   DrawDetectorStrips() const;
 
 protected:

  void MaterialMixer(Float_t* p,Float_t* a,Float_t* m,Float_t* d,Float_t* s,Int_t n) const;

private:
  Int_t fIdFTOA; // FTOA volume identifier (outer plate A)
  Int_t fIdFTOB; // FTOB volume identifier (outer plate B)
  Int_t fIdFTOC; // FTOC volume identifier (outer plate C)
  Int_t fIdFLTA; // FLTA volume identifier (inner plate A)
  Int_t fIdFLTB; // FLTB volume identifier (inner plate B)
  Int_t fIdFLTC; // FLTC volume identifier (inner plate C)
  Bool_t fTOFHoles; // Selecting Geometry with and w/o holes
 
  ClassDef(AliTOFv5T0,0)  //Time Of Flight version 5
};
 
#endif /* ALITOFv5T0_H */
