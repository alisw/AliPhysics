#ifndef ALIMUONVTEMP_H
#define ALIMUONVTEMP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
 
#include "AliMUONv1.h"

class AliMUONvTemp : public AliMUONv1 {
public:
    AliMUONvTemp(){;}
    AliMUONvTemp(const char *name, const char *title);
    
    virtual  ~AliMUONvTemp() {}
    virtual void   CreateGeometry();
private:
   ClassDef(AliMUONvTemp,1)  // MUON Detector class Version Temporary
};
#endif







