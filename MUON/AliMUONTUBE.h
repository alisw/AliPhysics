#ifndef AliMUONTube_H
#define AliMUONTube_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:MUON     //
////////////////////////////////////////////////
#include "TTUBE.h"

class AliMUONTUBE :
public TTUBE {
    
 public:
    AliMUONTUBE();
    AliMUONTUBE(Text_t *name, Text_t *title, Text_t *material, Float_t rmin, Float_t rmax, Float_t dz,Float_t aspect);
    AliMUONTUBE(Text_t *name, Text_t *title, Text_t *material, Float_t rmax, Float_t dz);
    virtual ~AliMUONTUBE(){};
    virtual void Paint(Option_t* option);
    ClassDef(AliMUONTUBE,1)
};
#endif
