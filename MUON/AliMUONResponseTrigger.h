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

class AliMUONResponseTrigger : public AliMUONResponse
{
 public:
  AliMUONResponseTrigger();
  virtual ~AliMUONResponseTrigger(){} 

  virtual Int_t DigitResponse(Int_t digit, AliMUONTransientDigit* where) const;


  // Set the GenerCluster parameter       
  virtual Int_t SetGenerCluster(){return 0;}

  virtual void DisIntegrate(const AliMUONHit& hit, TList& digits);
  
 private:
  static const Float_t fgkTofLimit;  ///< Particle above this threshold are discarded
  ClassDef(AliMUONResponseTrigger,1) ///< Implementation of RPC response
    
};
#endif













