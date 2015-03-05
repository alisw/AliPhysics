#ifndef ALIFRAMEV3_H
#define ALIFRAMEV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and class for detector: FRAME  version 2    //
/////////////////////////////////////////////////////////
 
#include "AliFRAME.h"
#include "TGeoCompositeShape.h"

class AliFRAMEv3 : public AliFRAME {
  
public:
  AliFRAMEv3();
  AliFRAMEv3(const char *name, const char *title);
  virtual       ~AliFRAMEv3() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   AddAlignableVolumes() const;
  virtual void   Init();
  virtual void   StepManager();
  virtual Int_t  IsVersion() const;
  virtual void   SetHoles(Int_t flag=0) {fHoles = flag;}
  virtual Int_t  Holes() const {return fHoles;}
  virtual void   MakeHeatScreen(const char* name, Float_t dyP, Int_t rot1, Int_t rot2);
  virtual void   WebFrame(const char* name, Float_t dHz, Float_t theta0, Float_t phi0);
  virtual TGeoCompositeShape* CreateTOFRail (Float_t y);

 private:
  Int_t  fHoles; // flag fHoles=0 => no holes, with holes otherwise
  
   ClassDef(AliFRAMEv3,2)  //Class for FRAME version 3
};
 
#endif
