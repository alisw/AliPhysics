#ifndef RICHSegResCkv_H
#define RICHSegResCkv_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRICHSegResV0.h"

class AliRICHResponseCkv : public AliRICHResponseV0 {
    
 public:
    AliRICHResponseCkv(){}
    virtual ~AliRICHResponseCkv(){}
    
    virtual Float_t IntPH(Float_t =0);
    
    ClassDef(AliRICHResponseCkv,1)
	
	
	};
	
#endif
