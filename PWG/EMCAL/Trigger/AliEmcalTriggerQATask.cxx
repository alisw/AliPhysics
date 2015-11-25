/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <TClonesArray.h>
#include <THashList.h>

#include "AliEmcalTriggerPatchInfoAPV1.h"
#include "AliEmcalTriggerQAAP.h"

#include "AliEmcalTriggerQATask.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerQATask)
/// \endcond

/**
 * Dummy constructor
 */
AliEmcalTriggerQATask::AliEmcalTriggerQATask() : 
  AliAnalysisTaskEmcal("AliEmcalTriggerQATask",kTRUE),
  fCaloTriggersInName("EmcalTriggers"),
  fEMCALTriggerQA(0),
  fCaloTriggersIn(0)
{
}

/**
 * Named constructor.
 * \param name Name of the trigger QA task
 */
AliEmcalTriggerQATask::AliEmcalTriggerQATask(const char *name) :
  AliAnalysisTaskEmcal(name,kTRUE),
  fCaloTriggersInName("EmcalTriggers"),
  fEMCALTriggerQA(0),
  fCaloTriggersIn(0)
{
  // Constructor.
  SetMakeGeneralHistograms(kTRUE);

  TString qaName(Form("%s_AliEmcalTriggerQAAP",name));
  fEMCALTriggerQA = new AliEmcalTriggerQAAP(qaName);
}

/**
 * Destructor
 */
AliEmcalTriggerQATask::~AliEmcalTriggerQATask()
{
  delete fEMCALTriggerQA;
}

/**
 * Init the analysis.
 */
void AliEmcalTriggerQATask::ExecOnce()
{
  AliAnalysisTaskEmcal::ExecOnce();

  if (!fInitialized) return;

  fCaloTriggersIn = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloTriggersInName));

  TString objname(fCaloTriggersIn->GetClass()->GetName());
  TClass cls(objname);
  if (!cls.InheritsFrom("AliEmcalTriggerPatchInfoAPV1")) {
    AliError(Form("%s: Objects of type %s in %s are not inherited from AliEmcalTriggerPatchInfoAPV1!",
		  GetName(), cls.GetName(), fCaloTriggersInName.Data())); 
    fCaloTriggersIn = 0;
  }

  if (!fCaloTriggersIn) {
    fInitialized = kFALSE;
    AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), fCaloTriggersInName.Data()));
    return;
  }
}

/**
 * Create objects, histograms
 */
void AliEmcalTriggerQATask::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fEMCALTriggerQA->SetDebugLevel(DebugLevel());
  fEMCALTriggerQA->Init();
  
  if (fOutput) {  
    fOutput->Add(fEMCALTriggerQA->GetListOfHistograms());
    PostData(1, fOutput);
  }
}

/**
 * Run analysis.
 * \return Always true.
 */
Bool_t AliEmcalTriggerQATask::Run() 
{
  return kTRUE;
}


/**
 * Fill QA histograms
 * \return Always true.
 */
Bool_t AliEmcalTriggerQATask::FillHistograms() 
{
  Int_t nPatches = fCaloTriggersIn->GetEntriesFast();

  AliDebug(2, Form("nPatches = %d", nPatches));

  Int_t type = 0;
  
  for (Int_t i = 0; i < nPatches; i++) {
    AliDebug(2, Form("Processing patch %d", i));

    AliEmcalTriggerPatchInfoAPV1* patch = static_cast<AliEmcalTriggerPatchInfoAPV1*>(fCaloTriggersIn->At(i));
    if (!patch) continue;

    fEMCALTriggerQA->ProcessPatch(patch);
  }

  fEMCALTriggerQA->EventCompleted();

  return kTRUE;
}
