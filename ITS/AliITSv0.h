#ifndef ITSv0_H
#define ITSv0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 0    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"

class AliITSv0 : public AliITS {
 
public:
   AliITSv0();
   AliITSv0(const char *name, const char *title);
   virtual       ~AliITSv0() {}
   virtual void   CreateGeometry();
   virtual Int_t  IsVersion() const {return 0;}
   virtual void   DrawModule();
   virtual void   StepManager();
   
   ClassDef(AliITSv0,1)  //Hits manager for set:ITS version 0
};
 
#endif
