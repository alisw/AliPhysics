#ifndef ALIITSV5SYMM_H
#define ALIITSV5SYMM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//     Manager and hits classes for  ITS version 5 with symm services
////////////////////////////////////////////////////////////////////////


#include "AliITS.h"

class TBRIK;
class AliITSv5symm : public AliITS {

 public:
    AliITSv5symm();
    AliITSv5symm(const char *name, const char *title);
    AliITSv5symm(const AliITSv5symm &source); // copy constructor
    AliITSv5symm& operator=(const AliITSv5symm &source); // assignment operator	 
    virtual       ~AliITSv5symm() ;
    virtual void   BuildGeometry();
    virtual void  CreateGeometry();
    virtual void  CreateMaterials();
    virtual void  Init();   
    virtual Int_t IsVersion() const {
      // returns the ITS version number
      return 5;
    }
    virtual void  StepManager();


    ClassDef(AliITSv5symm,1)//Hits manager for ITS version 5 Official detailed 
	                         //geometry with symmetric services
};
 
#endif
