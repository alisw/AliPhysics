#ifndef ALIITSVPPRSYMM_H
#define ALIITSVPPRSYMM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 7    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSvPPRsymm : public AliITS {

 public:
    AliITSvPPRsymm();
    AliITSvPPRsymm(const char *name, const char *title);
    AliITSvPPRsymm(const AliITSvPPRsymm &source); // copy constructor
    AliITSvPPRsymm& operator=(const AliITSvPPRsymm &source); // assignment operator
    virtual       ~AliITSvPPRsymm() ;
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init(); 
    virtual Int_t  IsVersion() const {
      // returns the ITS version number 
      return 9;
    } 
    virtual void   DrawModule();
    virtual void   StepManager();

    ClassDef(AliITSvPPRsymm,1)  //Hits manager for set:ITS version 9 
                                // PPR detailed Geometry symmetric
};
 
#endif
