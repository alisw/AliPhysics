#ifndef ALIMUONRESPONSETRIGGER_H
#define ALIMUONRESPONSETRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup sim
/// \class AliMUONResponseTrigger
/// \brief Implementation of RPC response

#include "AliMUONResponse.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerChamberEfficiency.h"

class AliMUONResponseTrigger : public AliMUONResponse
{
 public:
  AliMUONResponseTrigger();
  virtual ~AliMUONResponseTrigger(); 

  /// Set the GenerCluster parameter       
  virtual Int_t SetGenerCluster(){return 0;}

  virtual void DisIntegrate(const AliMUONHit& hit, TList& digits, Float_t timeDif);
    
 private:
  /// Not implemented
  AliMUONResponseTrigger(const AliMUONResponseTrigger& other);
  /// Not implemented
  AliMUONResponseTrigger& operator=(const AliMUONResponseTrigger& other); // assignment operator
  

  ClassDef(AliMUONResponseTrigger,3) // Implementation of RPC response
    
};
#endif
