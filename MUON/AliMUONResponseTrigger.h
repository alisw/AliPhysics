#ifndef ALIMUONRESPONSETRIGGER_H
#define ALIMUONRESPONSETRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#include "AliMUONResponse.h"

class AliMUONResponseTrigger : public AliMUONResponse
{
 public:
  AliMUONResponseTrigger();
  virtual ~AliMUONResponseTrigger(){} 

  virtual Int_t DigitResponse(Int_t digit, AliMUONTransientDigit* where);


  // Set the GenerCluster parameter       
  virtual Int_t SetGenerCluster(){return 0;}

  ClassDef(AliMUONResponseTrigger,1) // Implementation of RPC response
    
};
#endif













