#ifndef ALITRDTRIGGERL1_H
#define ALITRDTRIGGERL1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDTriggerL1.h 31443 2009-03-12 14:56:21Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD trigger implementation for L1 (GTU) simulation steering            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTriggerDetector.h"

class TObjArray;

class AliTRDTriggerL1 : public AliTriggerDetector {

 public:
  AliTRDTriggerL1();
  ~AliTRDTriggerL1();

  virtual void CreateInputs();
  virtual void Trigger();

 private:

  ClassDef(AliTRDTriggerL1, 1);

};

#endif
