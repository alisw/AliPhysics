#ifndef ALIPYTHIA_H
#define ALIPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TPythia6.h>
#include "GenTypeDefs.h"

class AliPythia:public TPythia6
{
 public:
    virtual ~AliPythia(){;}
    // convert to compressed code and print result (for debugging only)
    virtual Int_t CheckedLuComp(Int_t kf);
    // Pythia initialisation for selected processes
    virtual void ProcInit
	(Process_t process, Float_t energy, StrucFunc_t strucfunc);
    static  AliPythia* Instance();

 protected:
    Process_t     fProcess;           // Process type
    Float_t       fEcms;              // Centre of mass energy
    StrucFunc_t   fStrucFunc;         // Structure function
    static AliPythia*    fgAliPythia; // Pointer to single instance
 private: 
    AliPythia();

    ClassDef(AliPythia,1) //ALICE UI to PYTHIA
};

#endif



