#ifndef ALIITSVPPRASYMM_H
#define ALIITSVPPRASYMM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 6    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSvPPRasymm : public AliITS {

 public:
    AliITSvPPRasymm();
    AliITSvPPRasymm(const char *name, const char *title);
    AliITSvPPRasymm(const AliITSvPPRasymm &source); // copy constructor
    AliITSvPPRasymm& operator=(const AliITSvPPRasymm &source); // assignment operator
    virtual       ~AliITSvPPRasymm() ;
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init(); 
    virtual Int_t  IsVersion() const {
      // returns the ITS version number 
      return 8;
    } 
    virtual void   DrawModule();
    virtual void   StepManager();

    ClassDef(AliITSvPPRasymm,1)  //Hits manager for set:ITS version 8 
                                 // PPR detailed Geometry asymmetric
};
 
#endif
