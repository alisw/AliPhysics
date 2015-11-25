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
#include "AliEmcalTriggerFastORAP.h"

#include "AliEmcalTriggerQATask.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerQATask)
/// \endcond

/**
 * Dummy constructor
 */
AliEmcalTriggerQATask::AliEmcalTriggerQATask() : 
  AliAnalysisTaskEmcal("AliEmcalTriggerQATask",kTRUE),
  fTriggerPatchesName("EmcalTriggers"),
  fEMCALTriggerQA(0),
  fBadChannels(),
  fTriggerPatches(0)
{
}

/**
 * Named constructor.
 * \param name Name of the trigger QA task
 */
AliEmcalTriggerQATask::AliEmcalTriggerQATask(const char *name) :
  AliAnalysisTaskEmcal(name,kTRUE),
  fTriggerPatchesName("EmcalTriggers"),
  fEMCALTriggerQA(0),
  fBadChannels(),
  fTriggerPatches(0)
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

  fTriggerPatches = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTriggerPatchesName));

  TString objname(fTriggerPatches->GetClass()->GetName());
  TClass cls(objname);
  if (!cls.InheritsFrom("AliEmcalTriggerPatchInfoAPV1")) {
    AliError(Form("%s: Objects of type %s in %s are not inherited from AliEmcalTriggerPatchInfoAPV1!",
		  GetName(), cls.GetName(), fTriggerPatchesName.Data()));
    fTriggerPatches = 0;
  }

  if (!fTriggerPatches) {
    fInitialized = kFALSE;
    AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), fTriggerPatchesName.Data()));
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
  if (fTriggerPatches) {
    Int_t nPatches = fTriggerPatches->GetEntriesFast();

    AliDebug(2, Form("nPatches = %d", nPatches));

    Int_t type = 0;

    for (Int_t i = 0; i < nPatches; i++) {
      AliDebug(2, Form("Processing patch %d", i));

      AliEmcalTriggerPatchInfoAPV1* patch = static_cast<AliEmcalTriggerPatchInfoAPV1*>(fTriggerPatches->At(i));
      if (!patch) continue;

      fEMCALTriggerQA->ProcessPatch(patch);
    }
  }

  if (fCaloTriggers) {
    AliEmcalTriggerFastORAP fastor;
    fCaloTriggers->Reset();
    Int_t globCol = -1, globRow = -1;
    Float_t L0amp = -1;
    Int_t L1amp = -1;
    while (fCaloTriggers->Next()) {
      // get position in global 2x2 tower coordinates
      // A0 left bottom (0,0)
      fCaloTriggers->GetPosition(globCol, globRow);
      // exclude channel completely if it is masked as hot channel
      if (fBadChannels.HasChannel(globCol, globRow)) continue;
      // for some strange reason some ADC amps are initialized in reconstruction
      // as -1, neglect those
      fCaloTriggers->GetL1TimeSum(L1amp);
      if (L1amp < 0) L1amp = 0;
      fCaloTriggers->GetAmplitude(L0amp);
      if (L0amp < 0) L0amp = 0;

      fastor.Initialize(L0amp, L1amp, globRow, globCol, fGeom);

      fEMCALTriggerQA->ProcessFastor(&fastor);
    }
  }
  fEMCALTriggerQA->EventCompleted();

  return kTRUE;
}
