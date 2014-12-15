/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONSDigitizerV2
/// \brief MUON SDigitizer (from Hits to SDigits).
///
/// New sdigitizer, not deriving from MUONDigitizer, and using
/// new Response:DisIntegrate method
/// Note also that this one does *not* merge sdigits at all 
/// (this is deferred to the digitizer, which anyway has to do it), 
/// thus speeding a little bit this step.
///
/// \author Laurent Aphecetche

#ifndef ALIMUONSDIGITIZERV2_H
#define ALIMUONSDIGITIZERV2_H

#include "TNamed.h"

class AliMUONSDigitizerV2 : public TNamed
{
public:
  AliMUONSDigitizerV2();
  virtual ~AliMUONSDigitizerV2();
  
  virtual void Digitize(Option_t* opt="");
  
private:
  static Float_t  fgkMaxIntTime; ///< maximum time of interaction
  static Float_t  fgkMaxPosTimeDif; ///< maximum event time after the triggered event for a hit to be digitized 
  static Float_t  fgkMaxNegTimeDif; ///< maximum event time before the triggered event for a hit to be digitized 
  static Float_t  fgkMinTimeDif; ///< minimum time difference for the reduction factor to be applied  
    
  ClassDef(AliMUONSDigitizerV2,2) // MUON SDigitizer V2-1
};

#endif
