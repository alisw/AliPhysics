#ifndef RICHSegResCkv_H
#define RICHSegResCkv_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRICHSegResV0.h"

class AliRICHresponseCkv : public AliRICHresponseV0 {
    
 public:
    AliRICHresponseCkv(){}
    virtual ~AliRICHresponseCkv(){}
    
    virtual Float_t IntPH(Float_t =0);
    
    ClassDef(AliRICHresponseCkv,1)
	
	
	};
	
#endif
