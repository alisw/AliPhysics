#ifndef ITSv4_H
#define ITSv4_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 4    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSv4 : public AliITS {
 
protected:
  Int_t fMinorVersion;  //Minor version identifier
 
public:
   AliITSv4();
   AliITSv4(const char *name, const char *title);
   virtual       ~AliITSv4() {}
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   Init();   
   virtual Int_t  IsVersion() const {return 4;}
   virtual void   SetMinorVersion(Int_t version) {fMinorVersion=version;}
   virtual void   StepManager();
   
   ClassDef(AliITSv4,1)  //Hits manager for set:ITS version 4
};
 
#endif
