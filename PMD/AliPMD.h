#ifndef PMD_H
#define PMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:PMD      //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"


class AliPMD : public AliDetector {
  
protected:
  Float_t fPar[4];           // pmdin, pmdout, thgas, thcell
  Float_t fIn[5];            // thmin, thmax, zdist, thlow, thhigh
  Float_t fGeo[3];           // wafer, edge, numqu
  Float_t fPadSize[4];       // size of the pads
  Int_t   fNumPads[4];       // number of the pads
public:
  AliPMD();
  AliPMD(const char *name, const char *title);
  virtual      ~AliPMD() {}
  virtual void  AddHit(Int_t, Int_t*, Float_t*);
   virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials() {}
  Int_t         DistancetoPrimitive(Int_t, Int_t);
  virtual Int_t IsVersion() const =0;
  virtual void  SetPAR(Float_t, Float_t, Float_t, Float_t);
  virtual void  SetIN(Float_t, Float_t, Float_t, Float_t, Float_t);
  virtual void  SetGEO(Float_t, Float_t, Float_t);
  virtual void  SetPadSize(Float_t, Float_t, Float_t, Float_t);
  virtual void  StepManager();
  
  ClassDef(AliPMD,1)  // Base Class for Photon Multiplicity Detector
};

 
 
//___________________________________________
 
class AliPMDhit : public AliHit {
public:
  Int_t      fVolume[5];  //array of volumes
  Float_t    fEnergy;     //Total energy deposited in eV
public:
  AliPMDhit() {}
  AliPMDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  AliPMDhit(AliPMDhit* oldhit) {*this=*oldhit;}
  virtual ~AliPMDhit() {}
  inline virtual Int_t GetVolume(Int_t i) {return fVolume[i];}
  inline virtual Float_t GetEnergy() {return fEnergy;}
  inline int operator == (AliPMDhit &cell) {
    Int_t i;
    if(fTrack!=cell.GetTrack()) return 0;
    for (i=0; i<4; i++) if(fVolume[i]!=cell.GetVolume(i)) return 0;
    return 1;
  }
  inline virtual AliPMDhit& operator + (AliPMDhit &cell) {
    fEnergy+=cell.GetEnergy();
    return *this;
  }
  virtual void Print(Option_t *) {
    printf("PMD Cell %d %d %d %d\n   Primary %d -   Energy %f\n",
	   fVolume[0],fVolume[1],fVolume[2],fVolume[3],fTrack,fEnergy);
  }

 
  ClassDef(AliPMDhit,1)  //Hits object for set:PMD
};
#endif
