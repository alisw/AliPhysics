#ifndef ALIITSV5ASYMM_H
#define ALIITSV5ASYMM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//     Manager and hits classes for  ITS version 5 with symm services
////////////////////////////////////////////////////////////////////////


#include "AliITS.h"

class TBRIK;
class AliITSv5asymm : public AliITS {

 public:
    AliITSv5asymm();
    AliITSv5asymm(const char *name, const char *title);
    AliITSv5asymm(const AliITSv5asymm &source); // copy constructor
    AliITSv5asymm& operator=(const AliITSv5asymm &source); // assignment operator	 
    virtual       ~AliITSv5asymm() ;
    virtual void   BuildGeometry();
    virtual void  CreateGeometry();
    virtual void  CreateMaterials();
    virtual void  Init();   
    virtual Int_t IsVersion() const {
      // returns the ITS version number
      return 5;
    }
    virtual void  StepManager();


    ClassDef(AliITSv5asymm,1)//Hits manager for ITS version 5 Official detailed 
	                         //geometry with asymmetric services
};
 
#endif
