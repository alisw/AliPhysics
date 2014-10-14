#ifndef ALITRDTRIGGER_H
#define ALITRDTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDTrigger.h 31443 2009-03-12 14:56:21Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD trigger interface class to CTP                                     //
// from this class the two classes for L0 (pretrigger) and                //
// L1 (GTU) are called
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTriggerDetector.h"

class AliTRDTrigger : public AliTriggerDetector {

 public:
  AliTRDTrigger();
  ~AliTRDTrigger();

  virtual void    AssignInputs(const TObjArray& inputs);
  virtual void CreateInputs();
  virtual void Trigger();

 private:
  TObjArray fTriggers; // array of all contributing triggers

  ClassDef(AliTRDTrigger, 1);

};

#endif
