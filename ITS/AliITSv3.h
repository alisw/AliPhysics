#ifndef ITSv3_H
#define ITSv3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 3    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSv3 : public AliITS {

protected:
  Int_t fMinorVersionV3;  //Minor version identifier
 
public:
   AliITSv3();
   AliITSv3(const char *name, const char *title);
   virtual       ~AliITSv3() {}
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   Init();   
   virtual Int_t  IsVersion() const {return 3;}
   virtual inline void   SetMinorVersion(Int_t version) {fMinorVersionV3=version;}
   virtual void   StepManager();
   
   ClassDef(AliITSv3,1)  //Hits manager for set:ITS version 3
};
 
#endif
