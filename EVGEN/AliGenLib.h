#ifndef ALIGENLIB_H
#define ALIGENLIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id */

#include <TObject.h>
#include "GenTypeDefs.h"

class AliGenLib :
public TObject
{
 public:
//
    typedef Double_t (*GenFunc)  (Double_t *, Double_t *);
    typedef Int_t    (*GenFuncIp)();    
    virtual GenFunc   GetPt(Param_t param, const char *tname)   = 0;
    virtual GenFunc   GetY (Param_t param, const char *tname)  = 0;
    virtual GenFuncIp GetIp(Param_t param, const char *tname)  = 0;    
    ClassDef(AliGenLib,1) // Library providing y and pT parameterisations
};
#endif







