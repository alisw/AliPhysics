#ifndef ALIITSV1_H
#define ALIITSV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 1    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSv1 : public AliITS {

 public:
    AliITSv1();
    AliITSv1(const char *name, const char *title);
    AliITSv1(const AliITSv1 &source); // copy constructor
    AliITSv1& operator=(const AliITSv1 &source); // assignment operator
    virtual       ~AliITSv1() ;
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init(); 
    virtual Int_t  IsVersion() const {
      // returns the ITS version number 
      return 1;
    } 
    virtual void   DrawModule();
    virtual void   StepManager();

    ClassDef(AliITSv1,1)  //Hits manager for set:ITS version 1 cource Geometry
};
 
#endif
