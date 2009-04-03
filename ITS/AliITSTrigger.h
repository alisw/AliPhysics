#ifndef ALIITSTRIGGER_H
#define ALIITSTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//                                                                    //
// Simulates generation of Fast-OR signals from SPD (if needed).      //
// Processes the Fast-OR signals generated in AliITSsimulationSPD.    //
// Provides inputs for AliCentralTrigger.                             //
//                                                                    //
// Version 2, Henrik Tydesjo, Feb 2009                                //
// Version 1, D. Elia, C. Jorgensen, Mar 2006                         //
// Version 0, J. Conrad, E. Lopez Torres, Oct 2005                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "AliTriggerDetector.h"
#include "AliITSTriggerFOProcessor.h"

class AliITSTriggerConditions;

class AliITSTrigger : public AliTriggerDetector {

 public:
  AliITSTrigger();
  AliITSTrigger(AliITSTriggerConditions* cond);
  virtual ~AliITSTrigger() {}

  virtual void    SetTriggerConditions(AliITSTriggerConditions* cond);
  virtual void    CreateInputs();
  virtual void    Trigger();

 private:
  AliITSTriggerFOProcessor fPITprocessor; //! used for processing of FO signals

  ClassDef( AliITSTrigger, 2 )  // ITS SPD Trigger Detector class

};

#endif
