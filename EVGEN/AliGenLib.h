#ifndef ALIGENLIB_H
#define ALIGENLIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

class TRandom;

class AliGenLib :
  public TObject
{
 public:
//
    virtual ~AliGenLib(){}
    typedef Double_t (*GenFunc)  (const Double_t *, const Double_t *);
    typedef Int_t    (*GenFuncIp)(TRandom *);    
    virtual GenFunc   GetPt(Int_t param, const char *tname) const   = 0;
    virtual GenFunc   GetY (Int_t param, const char *tname) const   = 0;
    virtual GenFuncIp GetIp(Int_t param, const char *tname) const   = 0;    
    ClassDef(AliGenLib,0) // Library providing y and pT parameterisations
};
#endif







