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
#include <AliEMCALTriggerOnlineQAPbPb.h>
#include <AliEMCALTriggerQA.h>
#include <AliEMCALTriggerPatchInfo.h>
#include <AliEMCALTriggerFastOR.h>
#include <AliEMCALTriggerConstants.h>
#include <AliEMCALTriggerOnlineQAPP.h>
#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliEMCALTriggerOfflineQAPP.h"
#include "AliEmcalTriggerQATask.h"

using namespace EMCALTrigger;

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerQATask);
/// \endcond

/**
 * Dummy constructor
 */
AliEmcalTriggerQATask::AliEmcalTriggerQATask() :
  AliAnalysisTaskEmcalLight(),
  fTriggerPatchesName("EmcalTriggers"),
  fEMCALTriggerQA(),
  fADCperBin(16),
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
AliEmcalTriggerQATask::AliEmcalTriggerQATask(const char *name, UInt_t nCentBins, Bool_t online) :
  AliAnalysisTaskEmcalLight(name,kTRUE),
  fTriggerPatchesName("EmcalTriggers"),
  fEMCALTriggerQA(),
  fADCperBin(16),
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

  SetNCentBins(nCentBins);

  if (nCentBins == 0) {
    AliInfo("Setting up the task for pp collisions.");
    SetForceBeamType(AliAnalysisTaskEmcalLight::kpp);
    fEMCALTriggerQA.SetLast(0);
    if (online){
      fEMCALTriggerQA.AddAt(new AliEMCALTriggerOnlineQAPP(name), 0);
    }
    else {
      fEMCALTriggerQA.AddAt(new AliEMCALTriggerOfflineQAPP(name), 0);
    }
  }
  else {
    AliInfo("Setting up the task for PbPb collisions.");
    // No offline class for PbPb... yet
    SetForceBeamType(AliAnalysisTaskEmcalLight::kAA);
    fEMCALTriggerQA.SetLast(nCentBins-1);
    for (Int_t i = 0; i < nCentBins; i++) {
      fEMCALTriggerQA.AddAt(new AliEMCALTriggerOnlineQAPbPb(name), i);
    }
  }

  fEMCALTriggerQA.SetOwner(kTRUE);
}

/**
 * Destructor
 */
AliEmcalTriggerQATask::~AliEmcalTriggerQATask()
{
}

/**
 * Init the analysis.
 */
void AliEmcalTriggerQATask::ExecOnce()
{
  AliAnalysisTaskEmcalLight::ExecOnce();

  fESDEvent = dynamic_cast<AliESDEvent*>(InputEvent());

  if (!fESDEvent){
    fMinTimeStamp = 0;
    fMaxTimeStamp = 0;
    fTimeStampBinWidth = 0;
  }

  if (!fLocalInitialized) return;

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
    fLocalInitialized = kFALSE;
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
void AliEmcalTriggerQATask::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalLight::UserCreateOutputObjects();

  if (fOutput) {  
    for (Int_t i = 0; i < fNcentBins; i++) {
      GetTriggerQA(i)->EnableHistogramsByTimeStamp(fTimeStampBinWidth);
      GetTriggerQA(i)->SetDebugLevel(DebugLevel());
      GetTriggerQA(i)->Init();
      fOutput->Add(GetTriggerQA(i)->GetListOfHistograms());

      AliEMCALTriggerOfflineQAPP* triggerOffline = dynamic_cast<AliEMCALTriggerOfflineQAPP*>(GetTriggerQA(i));
      if (triggerOffline) {
        triggerOffline->EnableDCal(fDCalPlots);
      }
    }

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
  if (fESDEvent) {
    if (fESDEvent->GetTimeStamp() < fMinTimeStamp) return kFALSE;
    if (fMaxTimeStamp > fMinTimeStamp && fESDEvent->GetTimeStamp() > fMaxTimeStamp) return kFALSE;
    GetTriggerQA(fCentBin)->EventTimeStamp(fESDEvent->GetTimeStamp());
  }

  if (fTriggerPatches) {
    Int_t nPatches = fTriggerPatches->GetEntriesFast();

    AliDebug(2, Form("nPatches = %d", nPatches));

    for (Int_t i = 0; i < nPatches; i++) {
      AliDebug(2, Form("Processing bkg patch %d", i));

      AliEMCALTriggerPatchInfo* patch = static_cast<AliEMCALTriggerPatchInfo*>(fTriggerPatches->At(i));
      if (!patch) continue;
      if (patch->GetADCAmp() < fMinAmplitude) continue;

      GetTriggerQA(fCentBin)->ProcessBkgPatch(patch);
    }

    GetTriggerQA(fCentBin)->ComputeBackground();

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
    AliEMCALTriggerQA::AliEMCALCellInfo cellInfo;
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
void AliEmcalTriggerQATask::SetADCperBin(Int_t n)
{
  fADCperBin = n;

  TIter next(&fEMCALTriggerQA);
  AliEMCALTriggerQA* qa = 0;
  while ((qa = static_cast<AliEMCALTriggerQA*>(next()))) {
    qa->SetADCperBin(n);
  }
}

/**
 * Create a new instance of the AliEmcalTriggerQATask and adds it to the analysis manager.
 * \param triggerPatchesName name of the trigger patch collection
 * \param cellsName name of the EMCal cell collection
 * \param triggersName name of the primitive trigger objects (FastORs)
 * \param nCentBins number of centrality bins
 * \param online switch to use the online (HLT) or offline components
 * \param suffix to be added at the end of the task name
 * \param sudbir directory inside of the root file where the output objects will be stored
 * \return a pointer to the new instance of the class
 */
AliEmcalTriggerQATask* AliEmcalTriggerQATask::AddTaskEmcalTriggerQA(TString triggerPatchesName, TString cellsName, TString triggersName, Int_t nCentBins, Bool_t online, TString subdir, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliEmcalTriggerQATask", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AliEmcalTriggerQATask", "This task requires an input event handler");
    return 0;
  }

  // Init the task and do settings
  TString taskName("AliEmcalTriggerQATask");
  if (!suffix.IsNull()) {
    taskName += "_";
    taskName += suffix;
  }
  AliEmcalTriggerQATask* eTask = new AliEmcalTriggerQATask(taskName, nCentBins, online);
  eTask->SetTriggerPatchesName(triggerPatchesName);
  if(triggersName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      triggersName = "EMCALTrigger";
      ::Info("AddTaskEmcalTriggerQA", "ESD analysis, triggersName = \"%s\"", triggersName.Data());
    }
    else {
      triggersName = "emcalTrigger";
      ::Info("AddTaskEmcalTriggerQA", "AOD analysis, triggersName = \"%s\"", triggersName.Data());
    }
  }
  if(cellsName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      cellsName = "EMCALCells";
      ::Info("AddTaskEmcalTriggerQA", "ESD analysis, cellsName = \"%s\"", cellsName.Data());
    }
    else {
      cellsName = "emcalCells";
      ::Info("AddTaskEmcalTriggerQA", "AOD analysis, cellsName = \"%s\"", cellsName.Data());
    }
  }
  eTask->SetCaloTriggersName(triggersName);
  eTask->SetCaloCellsName(cellsName);

  // Final settings, pass to manager and set the containers
  mgr->AddTask(eTask);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );

  TString commonoutput;
  if (subdir.IsNull()) {
    commonoutput = mgr->GetCommonFileName();
  }
  else {
    commonoutput = TString::Format("%s:%s", mgr->GetCommonFileName(), subdir.Data());
  }
  TString contOutName(Form("%s_histos", taskName.Data()));
  mgr->ConnectOutput(eTask, 1, mgr->CreateContainer(contOutName, TList::Class(), AliAnalysisManager::kOutputContainer, commonoutput.Data()));

  return eTask;
}

/**
 * Add this task to the QA train
 * \param runnumber Run number
 */
void AliEmcalTriggerQATask::AddTaskEmcalTriggerQA_QAtrain(Int_t runnumber)
{
  EBeamType_t beam = BeamTypeFromRunNumber(runnumber);
  Int_t nCentBins = 0;
  if (beam == kpA || beam == kAA) nCentBins = 4;
  std::vector<std::string> triggerClasses = {"CINT7", "CEMC7", "CDMC7", "EG1", "EG2", "EJ1", "EJ2", "DG1", "DG2", "DJ1", "DJ2" };
  for (auto triggerClass : triggerClasses) {
    TString suffix(triggerClass.c_str());
    suffix.ReplaceAll("-", "_");
    AliEmcalTriggerQATask* task = AddTaskEmcalTriggerQA("EmcalTriggers", "", "", nCentBins, kFALSE, "CaloQA_default", suffix);
    task->AddAcceptedTriggerClass(triggerClass.c_str());
    task->SetForceBeamType(beam);
  }
}
