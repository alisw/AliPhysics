#ifndef ALIITSV3_H
#define ALIITSV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 3    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"

class AliITSv3 : public AliITS {

 public:
    AliITSv3();
    AliITSv3(const char *name, const char *title);
    AliITSv3(const AliITSv3 &source); // copy constructor
    AliITSv3& operator=(const AliITSv3 &source); // assignment operator
    virtual       ~AliITSv3() ;
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init();   
    virtual Int_t  IsVersion() const {
      // returns the ITS version number 
      return 3;
    }
    virtual void   SetMinorVersion(Int_t version) {
      // sets the minor version 
      fMinorVersionV3=version;
    }
    virtual void   StepManager();

 protected:
    Int_t fMinorVersionV3;  //Minor version identifier

   
    ClassDef(AliITSv3,1)//Hits manager for set:ITS version 3, TP detailed geometry
};
 
#endif

