#ifndef ALIMUONRESPONSETRIGGER_H
#define ALIMUONRESPONSETRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMUONResponseV0.h"

class AliMUONResponseTrigger : 
public AliMUONResponseV0 {
 public:
  AliMUONResponseTrigger(){};
  virtual ~AliMUONResponseTrigger(){} 

  virtual Int_t    DigitResponse(Int_t digit);    

  // Set the GenerCluster parameter       
  virtual Int_t SetGenerCluster();

  ClassDef(AliMUONResponseTrigger,1) // Implementation of RPC response
    
};
#endif













