#ifndef ITSv1_H
#define ITSv1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 1    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSv1 : public AliITS {

private:
    Int_t fId1N; // The number of layers for geometry version 5
    // The name of the layers as defined in the Geant tree.
    char  **fId1Name;
 
public:
  AliITSv1();
  AliITSv1(const char *name, const char *title);
  virtual       ~AliITSv1() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init(); 
  virtual Int_t  IsVersion() const {return 1;}
  virtual void   DrawModule();
  virtual void   StepManager();
  
   ClassDef(AliITSv1,1)  //Hits manager for set:ITS version 1
};
 
#endif
