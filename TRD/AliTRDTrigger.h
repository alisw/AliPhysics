#ifndef ALITRDTRIGGER_H
#define ALITRDTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDTrigger.h 31443 2009-03-12 14:56:21Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD trigger interface class to CTP                                     //
// currently the Trigger() method calls the GTU tracking simulation and   //
// runs two example triggers, namely on a single high pt particle and     //
// on a jet.                                                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTriggerDetector.h"

class AliTRDTrigger : public AliTriggerDetector {

 public:
  AliTRDTrigger();
  ~AliTRDTrigger();

  virtual void CreateInputs();
  virtual void Trigger();

 private:

  ClassDef(AliTRDTrigger, 1);

};

#endif
