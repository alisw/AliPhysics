#ifndef ALIITSVPPRSYMM_H
#define ALIITSVPPRSYMM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 7    //
/////////////////////////////////////////////////////////
 
#include "AliITSvPPRasymm.h"
 
class AliITSvPPRsymm : public AliITSvPPRasymm {

 public:
    AliITSvPPRsymm();
    AliITSvPPRsymm(const char *name, const char *title);
    virtual       ~AliITSvPPRsymm() ;
    virtual void   CreateGeometry();
    virtual Int_t  IsVersion() const {// returns the ITS version number 
	return 9;} 

 private:

    ClassDef(AliITSvPPRsymm,1)  //Hits manager for set:ITS version 9 
                                // PPR detailed Geometry symmetric
};
 
#endif
