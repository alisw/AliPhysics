#ifndef TOF_H
#define TOF_H
////////////////////////////////////////////////
//  Manager and hits classes for set:TOF     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"
 
 
class AliTOF : public AliDetector {

protected:
   Int_t fIdSens;
 
public:
  AliTOF();
  AliTOF(const char *name, const char *title);
  virtual        ~AliTOF() {}
  virtual void    AddHit(Int_t, Int_t*, Float_t*);
  virtual void    BuildGeometry();
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual void    Init();
  virtual Int_t   IsVersion() const =0;
  Int_t           DistancetoPrimitive(Int_t px, Int_t py);
  virtual void    StepManager()=0;
  virtual void    TOFpc(Float_t, Float_t, Float_t, Float_t, Float_t) {}
  virtual void    DrawModule();
  
  ClassDef(AliTOF,1)  // Time Of Flight base class
};
 
//___________________________________________
 
class AliTOFhit : public AliHit {
public:
  Int_t      fVolume[3];  //array of volumes
  Float_t    fPx;         // px in TOF
  Float_t    fPy;         // py in TOF
  Float_t    fPz;         // pz in TOF
  Float_t    fPmom;       // P in TOF
  Float_t    fTof;        // Time of Flight
 
public:
  AliTOFhit() {}
  AliTOFhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliTOFhit() {}
 
  ClassDef(AliTOFhit,1)  // Hits for Time Of Flight
};
 
#endif
