#ifndef FMD_H
#define FMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:FMD     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"
 
 
class AliFMD : public AliDetector {
 
public:
  AliFMD();
  AliFMD(const char *name, const char *title);
  virtual       ~AliFMD() {}
  virtual void   AddHit(Int_t, Int_t*, Float_t*);
  virtual void   BuildGeometry();
  virtual void   CreateGeometry() {}
  virtual void   CreateMaterials() {}
  Int_t          DistancetoPrimitive(Int_t, Int_t);
  virtual Int_t  IsVersion() const =0;
  virtual void   Init();
  virtual void   DrawModule()=0;
  virtual void   StepManager();
  
  ClassDef(AliFMD,1)  //Class for the FMD detector
};

//_____________________________________________________________________________
 
class AliFMDhit : public AliHit {
public:
  Int_t      fVolume;  //Volume copy identifier
  
public:
  AliFMDhit() {}
  AliFMDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliFMDhit() {}
  
  ClassDef(AliFMDhit,1)  //Hits for detector FMD
};

#endif
