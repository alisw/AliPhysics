#ifndef AliPoints_H
#define AliPoints_H

#include "TPolyMarker3D.h"
#include "AliDetector.h"
#include "TParticle.h"

class AliPoints : public TPolyMarker3D {
protected:
   AliDetector     *fDetector;    //Pointer to AliDetector object
   Int_t            fIndex;       //Particle number in AliRun::fParticles

public:
  AliPoints();
  AliPoints(Int_t nhits);
  virtual ~AliPoints();
  virtual Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual void          ExecuteEvent(Int_t event, Int_t px, Int_t py);
  AliDetector          *GetDetector() {return fDetector;}
  Int_t                 GetIndex() {return fIndex;}
  TParticle            *GetParticle() const;
  virtual const Text_t *GetName() const;
  virtual void          InspectParticle(); // *MENU*
  virtual void          DumpParticle(); // *MENU*
  virtual Text_t       *GetObjectInfo(Int_t px, Int_t py);
  virtual void          Propagate(); // *MENU*
  virtual void          SetDetector(AliDetector *det) {fDetector = det;}
  virtual void          SetParticle(Int_t index) {fIndex = index;}
  
  ClassDef(AliPoints,1) //Class to draw detector hits (is PolyMarker3D)
};
#endif
