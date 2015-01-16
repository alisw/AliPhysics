#ifndef ALITENDERSUPPLY_H
#define ALITENDERSUPPLY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 01/09/2009

//==============================================================================
//   AliTenderSupply - Base class for user-defined ESD additions and corrections.
//==============================================================================

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class AliTender;

class AliTenderSupply : public TNamed {

protected:
  const AliTender          *fTender;         // Tender car
  
public:  
  AliTenderSupply();
  AliTenderSupply(const char *name, const AliTender *tender=NULL);
  AliTenderSupply(const AliTenderSupply &other);
  virtual ~AliTenderSupply();
  AliTenderSupply& operator=(const AliTenderSupply &other);

  // Run control
  virtual void              Init() = 0;
  virtual void              ProcessEvent() = 0;
  
  void                      SetTender(const AliTender *tender) {fTender = tender;}
    
  ClassDef(AliTenderSupply,1)  // Base class for tender user algorithms
};
#endif
