/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <THnSparse.h>

#include <AliESDEvent.h>
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerQAPP.h"
#include "AliEMCALTriggerFastOR.h"
#include "AliEMCALTriggerConstants.h"

#include "AliEmcalTriggerQATaskPP.h"

using namespace EMCALTrigger;

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerQATaskPP)
/// \endcond

/**
 * Dummy constructor
 */
AliEmcalTriggerQATaskPP::AliEmcalTriggerQATaskPP() :
  AliAnalysisTaskEmcal(),
  fTriggerPatchesName("EmcalTriggers"),
  fEMCALTriggerQA(0),
  fADCperBin(20),
  fMinAmplitude(0),
  fDCalPlots(kTRUE),
  fMinTimeStamp(0),
  fMaxTimeStamp(0),
  fTimeStampBinWidth(0),
  fTriggerPatches(0),
  fESDEvent(0)
{
}

/**
 * Named constructor.
 * \param name Name of the trigger QA task
 */
AliEmcalTriggerQATaskPP::AliEmcalTriggerQATaskPP(const char *name) :
  AliAnalysisTaskEmcal(name,kTRUE),
  fTriggerPatchesName("EmcalTriggers"),
  fEMCALTriggerQA(0),
  fADCperBin(20),
  fMinAmplitude(0),
  fDCalPlots(kTRUE),
  fMinTimeStamp(0),
  fMaxTimeStamp(0),
  fTimeStampBinWidth(0),
  fTriggerPatches(0),
  fESDEvent(0)
{
  // Constructor.
  SetMakeGeneralHistograms(kTRUE);

  SetForceBeamType(AliAnalysisTaskEmcal::kpp);
  SetNCentBins(1);

  fEMCALTriggerQA = new TObjArray((fNcentBins+1)*2);
  fEMCALTriggerQA->SetOwner(kTRUE);
  fEMCALTriggerQA->AddAt(new AliEmcalTriggerQAPP(name), 0);
}

/**
 * Destructor
 */
AliEmcalTriggerQATaskPP::~AliEmcalTriggerQATaskPP()
{
  delete fEMCALTriggerQA;
}

/**
 * Init the analysis.
 */
void AliEmcalTriggerQATaskPP::ExecOnce()
{
  AliAnalysisTaskEmcal::ExecOnce();

  fESDEvent = dynamic_cast<AliESDEvent*>(InputEvent());

  if (!fESDEvent){
    fMinTimeStamp = 0;
    fMaxTimeStamp = 0;
    fTimeStampBinWidth = 0;
  }

  if (!fInitialized) return;

  fTriggerPatches = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTriggerPatchesName));

  if (fTriggerPatches) {
    TString objname(fTriggerPatches->GetClass()->GetName());
    TClass cls(objname);
    if (!cls.InheritsFrom("AliEMCALTriggerPatchInfo")) {
      AliError(Form("%s: Objects of type %s in %s are not inherited from AliEMCALTriggerPatchInfo!",
          GetName(), cls.GetName(), fTriggerPatchesName.Data()));
      fTriggerPatches = 0;
    }
  }

  if (!fTriggerPatches) {
    fInitialized = kFALSE;
    AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), fTriggerPatchesName.Data()));
    return;
  }

  for (Int_t i = 0; i < fNcentBins; i++) {
    GetTriggerQA(i)->ExecOnce();
  }
}

/**
 * Create objects, histograms
 */
void AliEmcalTriggerQATaskPP::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if (fOutput) {  
    for (Int_t i = 0; i < fNcentBins; i++) {
      GetTriggerQA(i)->EnableDCal(fDCalPlots);
      GetTriggerQA(i)->EnableHistogramsByTimeStamp(fTimeStampBinWidth);
      GetTriggerQA(i)->SetDebugLevel(DebugLevel());
      GetTriggerQA(i)->Init();
      fOutput->Add(GetTriggerQA(i)->GetListOfHistograms());
    }

    PostData(1, fOutput);
  }
}

/**
 * Run analysis.
 * \return Always true.
 */
Bool_t AliEmcalTriggerQATaskPP::Run()
{
  return kTRUE;
}


/**
 * Fill QA histograms
 * \return Always true.
 */
Bool_t AliEmcalTriggerQATaskPP::FillHistograms()
{
  if (fESDEvent) {
    if (fESDEvent->GetTimeStamp() < fMinTimeStamp) return kFALSE;
    if (fMaxTimeStamp > fMinTimeStamp && fESDEvent->GetTimeStamp() > fMaxTimeStamp) return kFALSE;
    GetTriggerQA(fCentBin)->EventTimeStamp(fESDEvent->GetTimeStamp());
  }

  if (fTriggerPatches) {
    Int_t nPatches = fTriggerPatches->GetEntriesFast();

    AliDebug(2, Form("nPatches = %d", nPatches));

    for (Int_t i = 0; i < nPatches; i++) {
      AliDebug(2, Form("Processing patch %d", i));

      AliEMCALTriggerPatchInfo* patch = static_cast<AliEMCALTriggerPatchInfo*>(fTriggerPatches->At(i));
      if (!patch) continue;
      if (patch->GetADCAmp() < fMinAmplitude) continue;

      GetTriggerQA(fCentBin)->ProcessPatch(patch);
    }
  }

  if (fCaloTriggers) {
    AliEMCALTriggerFastOR fastor;
    fCaloTriggers->Reset();
    Int_t globCol = -1, globRow = -1;
    Float_t L0amp = -1;
    Int_t L1amp = -1;
    while (fCaloTriggers->Next()) {
      // get position in global 2x2 tower coordinates
      // A0 left bottom (0,0)
      fCaloTriggers->GetPosition(globCol, globRow);
      // for some strange reason some ADC amps are initialized in reconstruction
      // as -1, neglect those
      fCaloTriggers->GetL1TimeSum(L1amp);
      if (L1amp < 0) L1amp = 0;
      fCaloTriggers->GetAmplitude(L0amp);
      L0amp *= 4;
      if (L0amp < 0) L0amp = 0;

      Int_t time = -1;
      Int_t nl0times(0);
      fCaloTriggers->GetNL0Times(nl0times);
       if(nl0times) {
         TArrayI l0times(nl0times);
         fCaloTriggers->GetL0Times(l0times.GetArray());
         for(int itime = 0; itime < nl0times; itime++){
           time = l0times[itime];
           break;
         }
       }

      fastor.Initialize(L0amp, L1amp, globRow, globCol, time, fGeom);

      GetTriggerQA(fCentBin)->ProcessFastor(&fastor, fCaloCells);
    }
  }

  if (fCaloCells) {
    const Int_t ncells = fCaloCells->GetNumberOfCells();
    AliEmcalTriggerQAPP::AliEmcalCellInfo cellInfo;
    for (Int_t pos = 0; pos < ncells; pos++) {
      Double_t amp = fCaloCells->GetAmplitude(pos);
      Int_t absId = fCaloCells->GetCellNumber(pos);
      cellInfo.Set(absId, amp);
      GetTriggerQA(fCentBin)->ProcessCell(cellInfo);
    }
  }

  GetTriggerQA(fCentBin)->EventCompleted();

  return kTRUE;
}

/**
 * Set number of ADC per bin in all the trigger QA
 * \param i number of ADC per bin.
 */
void AliEmcalTriggerQATaskPP::SetADCperBin(Int_t n)
{
  fADCperBin = n;

  for (Int_t i = 0; i < fNcentBins; i++) {
    GetTriggerQA(i)->SetADCperBin(n);
  }
}
