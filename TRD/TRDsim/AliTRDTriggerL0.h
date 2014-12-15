#ifndef ALITRDTRIGGERL0_H
#define ALITRDTRIGGERL0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDTriggerL0.h 31443 2009-03-12 14:56:21Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD trigger implementation for L0 (pretrigger) simulation              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTriggerDetector.h"

class TObjArray;

class AliTRDTriggerL0 : public AliTriggerDetector {

 public:
  AliTRDTriggerL0();
  ~AliTRDTriggerL0();

  virtual void CreateInputs();
  virtual void Trigger();

 private:

  ClassDef(AliTRDTriggerL0, 1);

};

#endif
