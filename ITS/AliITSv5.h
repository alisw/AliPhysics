#ifndef ALIITSV5_H
#define ALIITSV5_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//     Manager and hits classes for  ITS version 5
////////////////////////////////////////////////////////////////////////


#include "AliITS.h"

class TBRIK;
class AliITSv5 : public AliITS {

 public:
    AliITSv5();
    AliITSv5(const char *name, const char *title);
    AliITSv5(const AliITSv5 &source); // copy constructor
    AliITSv5& operator=(const AliITSv5 &source); // assignment operator	 
    virtual       ~AliITSv5() ;
    virtual void   BuildGeometry();
    virtual void  CreateGeometry();
    virtual void  CreateMaterials();
    virtual void  Init();   
    virtual Int_t IsVersion() const {
      // returns the ITS version number
      return 5;
    }
    virtual void  StepManager();


    ClassDef(AliITSv5,1)//Hits manager for ITS version 5 Official detailed geometry
};
 
#endif
