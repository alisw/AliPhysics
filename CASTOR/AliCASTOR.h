#ifndef CASTOR_H
#define CASTOR_H
////////////////////////////////////////////////
//  Manager and hits classes for set:CASTOR     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"
 
 
class AliCASTOR : public AliDetector {
 
public:
  AliCASTOR();
  AliCASTOR(const char *name, const char *title);
  virtual      ~AliCASTOR() {}
  virtual void  AddHit(Int_t, Int_t*, Float_t*);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials() {}
  virtual Int_t IsVersion() const =0;
  Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual void  DrawModule()=0;
  virtual void  StepManager()=0;
  
  ClassDef(AliCASTOR,1)  // CASTOR detector class
};
 
 
 
class AliCASTORv1 : public AliCASTOR {

private:

  Int_t fOdFiber;           //Tracking medium of fiber
  Int_t fOdCladding;        //Tracking medium of cladding
  Int_t fOdAbsorber;        //Tracking medium absorber
  Int_t fOctants;           //Number of azimutal sectors
  Int_t fLayersEM;          //Number of EM layers
  Int_t fLayersHad;         //Number of hadronic layers
  Float_t fPhiOct;          //Phi opening of azimutal sectors
  Float_t fRadCore;         //Radius of fiber core
  Float_t fRadFactor;       //Normalisation factor for radius
  
public:
  AliCASTORv1();
  AliCASTORv1(const char *name, const char *title);
  virtual      ~AliCASTORv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual void  DrawModule();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  StepManager();
 
   ClassDef(AliCASTORv1,1)  //Class for CASTOR version 1
};
 
 
 
//_____________________________________________________________________________
 
class AliCASTORhit : public AliHit {
public:
  Int_t      fVolume;     //array of volumes
 
public:
  AliCASTORhit() {}
  AliCASTORhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliCASTORhit() {}
 
  ClassDef(AliCASTORhit,1)  //Hits for CASTOR
};
 
#endif
