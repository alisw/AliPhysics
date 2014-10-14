/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliTRDTrigger.cxx 31904 2009-04-08 16:42:03Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TRD trigger interface class to CTP                                        //
// from here the trigger simulation for L0 (pretrigger) and L1 (GTU) are     //
// called                                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TClonesArray.h"

#include "AliLog.h"
#include "AliTriggerInput.h"
#include "AliTriggerDetector.h"

#include "AliTRDTrigger.h"
#include "AliTRDTriggerL0.h"
#include "AliTRDTriggerL1.h"

AliTRDTrigger::AliTRDTrigger() :
  AliTriggerDetector(),
  fTriggers()
{
  // default constructor

  fTriggers.AddLast(new AliTRDTriggerL0());
  fTriggers.AddLast(new AliTRDTriggerL1());

  SetName("TRD");
}

AliTRDTrigger::~AliTRDTrigger()
{
  // destructor
  TIter trigger(&fTriggers);
  while (AliTriggerDetector *trgDet = (AliTriggerDetector*) trigger())
    delete trgDet;

  fInputs.Clear(); // inputs are deleted either by CTP or submodule
}

void AliTRDTrigger::AssignInputs(const TObjArray& inputs)
{
  // Create inputs for all registered trigger modules.
  if( fInputs.GetEntriesFast() > 0 ) return;

  TIter trigger(&fTriggers);
  while (AliTriggerDetector *trgDet = (AliTriggerDetector*) trigger()) {
    trgDet->AssignInputs(inputs);
    fInputs.AddAll(trgDet->GetInputs());
  }
}

void AliTRDTrigger::CreateInputs()
{

}

void AliTRDTrigger::Trigger()
{
  // TRD trigger steering
  // all registered TRD trigger mechanism are
  // run from here

  TIter trigger(&fTriggers);
  while (AliTriggerDetector *trgDet = (AliTriggerDetector*) trigger()) {
    trgDet->Trigger();
  }
}
