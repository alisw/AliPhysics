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
  Float_t fPtThresholdA;	// pt threshold A
  Float_t fPtThresholdB;	// pt threshold B
  Int_t   fPidThresholdA;	// PID threshold A
  Int_t   fPidThresholdB;	// PID threshold B
  Int_t   fNoThreshold;		// number threshold for all tracks
  Int_t   fNoThresholdA;        // number threshold for tracks above pt A
  Int_t   fNoThresholdB;        // number threshold for tracks above pt B
  Int_t   fNoThresholdJetA;     // number threshold for tracks above pt A (jets)
  Int_t   fNoThresholdJetB;     // number threshold for tracks above pt B (jets)
  Int_t   fNoThresholdElA;      // number threshold for tracks above pt A and PID A (electron)
  Int_t   fNoThresholdElB;      // number threshold for tracks above pt B and PID B (electron)

  ClassDef(AliTRDTriggerL1, 2);

};

#endif
