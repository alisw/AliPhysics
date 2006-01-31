/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONSDigitizerV2
/// \brief New sdigitizer, not deriving from MUONDigitizer, and using
/// new Response:DisIntegrate method
///
/// Note also that this one does *not* merge sdigits at all 
/// (this is deferred to the digitizer, which anyway has to do it), 
/// thus speeding a little bit this step.
///
/// \author Laurent Aphecetche

#ifndef ALIMUONSDIGITIZERV2_H
#define ALIMUONSDIGITIZERV2_H

#ifndef ROOT_TTask
#  include "TTask.h"
#endif

class AliMUONSDigitizerV2 : public TTask
{
public:
  AliMUONSDigitizerV2();
  virtual ~AliMUONSDigitizerV2();
  
  virtual void Exec(Option_t* opt="");
    
  ClassDef(AliMUONSDigitizerV2,1) // 
};

#endif
