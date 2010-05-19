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

/* $Id: AliTRDTriggerL0.cxx 31904 2009-04-08 16:42:03Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TRD trigger L0 (pretrigger) simulation                                    //
// So far no real trigger decision is done.                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"

#include "AliLog.h"
#include "AliTriggerInput.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliTRDptrgParam.h"
#include "AliTRDptrgCBB.h"

#include "AliTRDTriggerL0.h"

AliTRDTriggerL0::AliTRDTriggerL0()
{
  SetName("TRD");
}

AliTRDTriggerL0::~AliTRDTriggerL0()
{

}

void AliTRDTriggerL0::CreateInputs()
{
  if (fInputs.GetEntriesFast() > 0)
    return;

  fInputs.AddLast(new AliTriggerInput("0HWU", "TRD", 1)); // TRD wake up
  fInputs.AddLast(new AliTriggerInput("0HSG", "TRD", 1)); // single gap
  fInputs.AddLast(new AliTriggerInput("0HDG", "TRD", 1)); // double gap
}

void AliTRDTriggerL0::Trigger()
{

  AliRunLoader *runLoader = AliRunLoader::Instance();
  if (!runLoader)
    return;
  AliLoader *trdLoader = runLoader->GetLoader("TRDLoader");
  if (!trdLoader)
    return;

  AliTRDptrgParam* param = AliTRDptrgParam::Instance();

  AliTRDptrgCBB* ptrgCBB = new AliTRDptrgCBB(runLoader, param, kDigits);

  Int_t* simulationResult;
  simulationResult = ptrgCBB->Simulate();
  for (Int_t iResult = 1; iResult <= simulationResult[0]; iResult++) {
    AliDebug(5, Form("Result[%d]=0x%x\n",iResult,simulationResult[iResult]));
  }
  if ((simulationResult[0] > 0) || (simulationResult[1] > 0)) { 
    AliInfo("Fired single gap trigger");
    SetInput("0HSG");
  }

  if (simulationResult[2] > 0) {
    AliInfo("Fired  double gap trigger");
    SetInput("0HDG");
  }

  if (simulationResult[3] > 0) {
    AliInfo("Fired TRD wake up call trigger");
    SetInput("0HWU");
  }

  delete ptrgCBB;
  if (simulationResult != 0x0)
    delete[] simulationResult;
  simulationResult = 0x0;

  AliDebug(5, Form("memory state: %d", param->CheckVariables()));
}
