#ifndef TRD_H
#define TRD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set: TRD     //
////////////////////////////////////////////////
 
#include "AliRun.h"
#include "AliDetector.h"
#include "AliHit.h" 
#include "AliDigit.h"

#include "AliTRDconst.h"
//#include "AliTRDgeometry.h"

class AliTRDgeometry;

//_____________________________________________________________________________
class AliTRD : public AliDetector {

 public:

  AliTRD();
  AliTRD(const char *name, const char *title);
  virtual           ~AliTRD();
  virtual void       AddHit(Int_t, Int_t, Float_t*);
  virtual void       AddDigit(Int_t*);    
  virtual void       AddRecPoint(Float_t*, Int_t*, Int_t, Float_t);
  virtual void       BuildGeometry();
  virtual void       CreateGeometry();
  virtual void       CreateMaterials();
  virtual void       DrawModule();
  Int_t              DistancetoPrimitive(Int_t px, Int_t py);
  TClonesArray      *RecPoints()           { return fRecPoints;   };
  virtual void       Init();
  virtual Int_t      IsVersion() const = 0;
  virtual void       MakeBranch(Option_t* option);     
  virtual void       ResetRecPoints();
  virtual void       StepManager() = 0; 
  virtual void       SetTreeAddress();

  virtual void       SetGasMix(Int_t imix = 0);
  virtual void       SetHits(Int_t ihit = 1) {};

  AliTRDgeometry    *GetGeometry()         { return fGeometry; };

  virtual Int_t      GetSensChamber() = 0;
  virtual Int_t      GetSensPlane()   = 0;
  virtual Int_t      GetSensSector()  = 0;
 
 protected:

  Int_t              fGasMix;            // Gas mixture. 0: Xe/Isobutane 1: Xe/CO2

  AliTRDgeometry    *fGeometry;          // The TRD geometry

  TClonesArray      *fRecPoints;         // List of reconstructed points
  Int_t              fNRecPoints;        //! Number of reconstructed points

  ClassDef(AliTRD,1)                     // Transition Radiation Detector base class

};

//_____________________________________________________________________________
class AliTRDhit : public AliHit {

public:
  Int_t        fDetector;   // TRD detector number
  Float_t      fQ;          // Charge created by a hit (slow simulator only)
 
public:
  AliTRDhit() {}
  AliTRDhit(Int_t shunt, Int_t track, Int_t det, Float_t *hits);
  virtual ~AliTRDhit() {};
 
  ClassDef(AliTRDhit,2)     // Hits for Transition Radiation Detector

};

#endif
