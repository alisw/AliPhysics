/*
 ***********************************************************
 Manager for event plane corrections framework
Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
Instructions in AddTask_EPcorrectionsExample.C
2014/12/10
 *********************************************************
 */

#include <Riostream.h>

#include <TGrid.h>
#include <TROOT.h>
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <TChain.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliESDEvent.h>
#include "AliQnCorrectionsCutsSet.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsHistos.h"
#include "AliLog.h"

#include "AliAnalysisTaskFlowVectorCorrections.h"

ClassImp(AliAnalysisTaskFlowVectorCorrections)


AliAnalysisTaskFlowVectorCorrections::AliAnalysisTaskFlowVectorCorrections() :
AliQnCorrectionsFillEventTask(),
fCalibrateByRun(kTRUE),
fCalibrationFile(""),
fCalibrationFileSource(CALIBSRC_local),
fTriggerMask(0),
fEventQAList(0x0),
fEventCuts(NULL),
fLabel(""),
fQAhistograms(""),
fFillEventQA(kFALSE),
fProvideQnVectorsList(kFALSE),
fOutputSlotEventQA(-1),
fOutputSlotHistQA(-1),
fOutputSlotHistNveQA(-1),
fOutputSlotHistQn(-1),
fOutputSlotQnVectorsList(-1),
fOutputSlotTree(-1)
{
  //
  // Default constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskFlowVectorCorrections::AliAnalysisTaskFlowVectorCorrections(const char* name) :
    AliQnCorrectionsFillEventTask(name),
fCalibrateByRun(kTRUE),
fCalibrationFile(""),
fCalibrationFileSource(CALIBSRC_local),
fTriggerMask(0),
fEventQAList(0x0),
fEventCuts(NULL),
fLabel(""),
fQAhistograms(""),
fFillEventQA(kFALSE),
fProvideQnVectorsList(kFALSE),
fOutputSlotEventQA(-1),
fOutputSlotHistNveQA(-1),
fOutputSlotHistQA(-1),
fOutputSlotHistQn(-1),
fOutputSlotQnVectorsList(-1),
fOutputSlotTree(-1)
{
  //
  // Constructor
  //

  fEventQAList = new TList();
  fEventQAList->SetName("EventQA");
  fEventQAList->SetOwner(kTRUE);

  fEventHistos = new AliQnCorrectionsHistos();
}

//_________________________________________________________________________________
void AliAnalysisTaskFlowVectorCorrections::DefineInOutput(){

  if(!fAliQnCorrectionsManager) {
    AliFatal("First configure QnCorrecionsManager!!\n");
    return;
  }

  DefineInput(0,TChain::Class());
  Int_t outputSlot = 1;
  // Calibration histograms
  if (fAliQnCorrectionsManager->GetShouldFillOutputHistograms()) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotHistQn = outputSlot++;
  }
  // Calibrated qvector tree
  if (fAliQnCorrectionsManager->GetShouldFillQnVectorTree()) {
    DefineOutput(outputSlot, TTree::Class());
    fOutputSlotTree = outputSlot++;
  }
  // Qvector QA histograms
  if (fAliQnCorrectionsManager->GetShouldFillQAHistograms()) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotHistQA = outputSlot++;
  }
  // Qvector non validated entries QA histograms
  if (fAliQnCorrectionsManager->GetShouldFillNveQAHistograms()) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotHistNveQA = outputSlot++;
  }
  // Calibrated qvector list
  if (fProvideQnVectorsList) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotQnVectorsList=outputSlot++;
  }
  // Event QA histograms
  if (fFillEventQA) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotEventQA=outputSlot++;
  }
}

void AliAnalysisTaskFlowVectorCorrections::SetCalibrationHistogramsFile(CalibrationFileSource source, const char *filename) {
  TFile *calibfile = NULL;

  fCalibrationFile = filename;
  fCalibrationFileSource = source;

  switch (fCalibrationFileSource) {
  case CALIBSRC_local:
    if (fCalibrationFile.Length() != 0) {
      if(fCalibrationFile.Contains("alien"))
        TGrid::Connect("alien://");
      calibfile = TFile::Open(fCalibrationFile);
    }
    if (calibfile != NULL && calibfile->IsOpen()) {
      AliInfo(Form("\t Calibration file %s open", fCalibrationFile.Data()));
      fAliQnCorrectionsManager->SetCalibrationHistogramsList(calibfile);
      calibfile->Close();
    }
    break;
  case CALIBSRC_alien:
    /* sanity check before we go to the grid */
    if (!fCalibrationFile.Contains("alien"))
      AliFatal(Form("\t alien was selected as source but %s filename does not contain \"alien\". Aborting!!!", fCalibrationFile.Data()));
    break;
  default:
    AliFatal("Calibration file source not supported. Aborting!!!");
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskFlowVectorCorrections::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //
  this->SetDefaultVarNames();
  this->SetDetectors();

  TFile *calibfile = NULL;

  /* get the calibration file if needed */
  switch (fCalibrationFileSource) {
  case CALIBSRC_local:
    break;
  case CALIBSRC_alien:
    if (fCalibrationFile.Length() != 0) {
      TGrid::Connect("alien://");
      calibfile = TFile::Open(fCalibrationFile);
    }
    if (calibfile != NULL && calibfile->IsOpen()) {
      AliInfo(Form("\t Calibration file %s open", fCalibrationFile.Data()));
      fAliQnCorrectionsManager->SetCalibrationHistogramsList(calibfile);
      calibfile->Close();
    }
    break;
  default:
    break;
  }

  fAliQnCorrectionsManager->InitializeQnCorrectionsFramework();

  if (fAliQnCorrectionsManager->GetShouldFillOutputHistograms())
    PostData(fOutputSlotHistQn, fAliQnCorrectionsManager->GetOutputHistogramsList());
  if (fAliQnCorrectionsManager->GetShouldFillQnVectorTree())
    PostData(fOutputSlotTree, fAliQnCorrectionsManager->GetQnVectorTree());
  if (fAliQnCorrectionsManager->GetShouldFillQAHistograms())
    PostData(fOutputSlotHistQA, fAliQnCorrectionsManager->GetQAHistogramsList());
  if (fAliQnCorrectionsManager->GetShouldFillNveQAHistograms())
    PostData(fOutputSlotHistNveQA, fAliQnCorrectionsManager->GetNveQAHistogramsList());
  if (fFillEventQA)
    PostData(fOutputSlotEventQA, fEventQAList);
}

/// The current run has changed. Usually it is only sent before
/// the first event is handled.
/// Notify the framework manager that the current label has changed.
void AliAnalysisTaskFlowVectorCorrections::NotifyRun() {

  if (fCalibrateByRun) fAliQnCorrectionsManager->SetCurrentProcessListName(Form("%d", this->fCurrentRunNumber));
}

void AliAnalysisTaskFlowVectorCorrections::UserExec(Option_t *){
  //
  // Main loop. Called for every event
  //

  fEvent = InputEvent();
  fAliQnCorrectionsManager->ClearEvent();

  fDataBank = fAliQnCorrectionsManager->GetDataContainer();

  FillEventData();

  fEventHistos->FillHistClass("Event_NoCuts", fDataBank);

  if (IsEventSelected(fDataBank)) {
    fEventHistos->FillHistClass("Event_Analysis", fDataBank);

    fAliQnCorrectionsManager->ProcessEvent();
  }  // end if event selection

  if(fProvideQnVectorsList)
    PostData(fOutputSlotQnVectorsList, fAliQnCorrectionsManager->GetQnVectorList());
}  // end loop over events


void AliAnalysisTaskFlowVectorCorrections::FinishTaskOutput()
{
  //
  // Finish Task
  //
  fAliQnCorrectionsManager->FinalizeQnCorrectionsFramework();

  THashList* hList = (THashList*) fEventHistos->HistList();
  for(Int_t i=0; i<hList->GetEntries(); ++i) {
    THashList* list = (THashList*)hList->At(i);
    fEventQAList->Add(list);
  }
}

Bool_t AliAnalysisTaskFlowVectorCorrections::IsEventSelected(Float_t* values) {

  if(!fEventCuts) return kTRUE;
  return fEventCuts->IsSelected(values);
}


