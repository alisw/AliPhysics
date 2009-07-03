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

#include "AliTRDTriggerL0.h"
#include "AliTRDgtuSim.h"
#include "AliTRDtrackGTU.h"

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

  fInputs.AddLast(new AliTriggerInput("0HMB", "TRD", 1)); // whatever should be there
}

void AliTRDTriggerL0::Trigger()
{
  // just an example:
  AliRunLoader *runLoader = AliRunLoader::Instance();
  if (!runLoader)
    return;
  AliLoader *trdLoader = runLoader->GetLoader("TRDLoader");
  if (!trdLoader)
    return;

  // here comes the actual pretrigger simulation

}
