#ifndef ZDC_H
#define ZDC_H
////////////////////////////////////////////////
//  Manager and hits classes for set:ZDC      //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"

 
class AliZDC : public AliDetector {

public:
  AliZDC();
  AliZDC(const char *name, const char *title);
  virtual      ~AliZDC() {}
  virtual void  AddHit(Int_t, Int_t*, Float_t*);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials() {}
  Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual Int_t IsVersion() const =0;
  virtual void  SetBeam(Int_t beam, Float_t fx, Float_t fy, Float_t sx, Float_t sy,
			Float_t div, Float_t angle, Int_t cross);
  virtual void  SetHijing(Int_t hij, Int_t hijf, Int_t hijsp, const char *file);
  virtual void  SetVenus(Int_t hiv, Int_t hivf, Int_t hivsp, const char *file);
  virtual void  SetKine(Int_t code, Float_t pmom, Float_t cx, Float_t cy, Float_t cz, Int_t type, Int_t fermi);
  virtual void  StepManager();
 
   ClassDef(AliZDC,1)  // Zero Degree Calorimeter base class
};
 
 
//____________________________________________________________________________ 
class AliZDCv1 : public AliZDC {

public:
  AliZDCv1();
  AliZDCv1(const char *name, const char *title);
  virtual      ~AliZDCv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  DrawModule();
 
   ClassDef(AliZDCv1,1)  // Zero Degree Calorimeter version 1
};
 
 
//_____________________________________________________________________________
class AliZDChit : public AliHit {
public:
  Int_t      fVolume[4];  //array of volumes
  Float_t    fEnergy;     //Total energy deposited in eV
 
public:
  AliZDChit() {}
  AliZDChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliZDChit() {}
 
  ClassDef(AliZDChit,1)  // Hits for the Zero Degree Calorimeter
};
 
#endif
