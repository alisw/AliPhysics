/*
 ***********************************************************
 Manager for event plane corrections framework
Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
Instructions in AddTask_EPcorrectionsExample.C
2014/12/10
 *********************************************************
 */

//#include "AliSysInfo.h"
#include <iostream>

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
#include "AliQnCorrectionsFillEvent.h"
#include "AliQnCorrectionsCuts.h"
#include "AliQnCorrectionsVarManager.h"
#include "AliQnCorrectionsManager.h"
//#include "QnCorrectionsHelper.h"
#include "AliQnCorrectionsHistos.h"
#include "TSystem.h"
#include "TROOT.h"
//#include "TObjectTable.h"
//#include "QnCorrectionsReducedVarManager.h"
//#include "AliReducedEvent.h"
//#include "trainsimulator/localHelper.h"
// make a change
#include "AliAnalysisTaskFlowVectorCorrections.h"

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskFlowVectorCorrections)


//_________________________________________________________________________________
AliAnalysisTaskFlowVectorCorrections::AliAnalysisTaskFlowVectorCorrections() :
    AliAnalysisTaskSE(),
    fRunLightWeight(kFALSE),
    fCalibrateByRun(kTRUE),
    fUseFriendEvent(kFALSE),
    fTriggerMask(0),
    fInitialized(kFALSE),
    fListInputHistogramsQnCorrections(0x0),
    fEventQAList(0x0),
    fEventPlaneManager(0x0),
    fEventCuts(0x0),
    fFillEvent(0x0),
    fEventPlaneHistos(0x0),
    fLabel(""),
    fQAhistograms(""),
    fFillEventQA(kTRUE),
    fProvideQnVectorsList(kTRUE),
    fOutputSlotEventQA(-1),
    fOutputSlotHistQA(-1),
    fOutputSlotHistQn(-1),
    fOutputSlotQnVectorsList(-1),
    fOutputSlotTree(-1),
    fRunListPath(""),
    fCalibrationFilePath("")
{
  //
  // Default constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskFlowVectorCorrections::AliAnalysisTaskFlowVectorCorrections(const char* name) :
  AliAnalysisTaskSE(name),
  fRunLightWeight(kFALSE),
  fCalibrateByRun(kTRUE),
  fUseFriendEvent(kFALSE),
  fTriggerMask(0),
  fInitialized(kFALSE),
  fListInputHistogramsQnCorrections(0x0),
  fEventQAList(0x0),
  fEventPlaneManager(0x0),
  fEventCuts(0x0),
  fFillEvent(0x0),
  fEventPlaneHistos(0x0),
  fLabel(""),
  fQAhistograms(""),
  fFillEventQA(kTRUE),
  fProvideQnVectorsList(kTRUE),
  fOutputSlotEventQA(-1),
  fOutputSlotHistQA(-1),
  fOutputSlotHistQn(-1),
  fOutputSlotQnVectorsList(-1),
  fOutputSlotTree(-1),
  fRunListPath(""),
  fCalibrationFilePath("")
{
  //
  // Constructor
  //

  //fFillEvent = new QnCorrectionsReducedVarManager();
  fFillEvent = new AliQnCorrectionsFillEvent();

  fEventQAList = new TList();
  fEventQAList->SetName("EventQA");
  fEventQAList->SetOwner(kTRUE);

  fEventPlaneHistos = new AliQnCorrectionsHistos();
  fFillEvent->SetEventPlaneHistos(fEventPlaneHistos);


}

//_________________________________________________________________________________
void AliAnalysisTaskFlowVectorCorrections::DefineInOutput(){

  if(!fEventPlaneManager) {std::cout<<"First configure EventPlaneManager!!"<<std::endl; return;}

  DefineInput(0,TChain::Class());
  Int_t outputSlot=1;
  if(fEventPlaneManager->ShouldFillHistogramsQnCorrections())  {DefineOutput(outputSlot, TList::Class()); fOutputSlotHistQn=outputSlot++;}   // Calibration histograms
  if(fEventPlaneManager->ShouldFillTreeQnVectors())            {DefineOutput(outputSlot, TTree::Class()); fOutputSlotTree=outputSlot++;} // Calibrated qvector tree
  if(fEventPlaneManager->ShouldFillHistogramsQA())             {DefineOutput(outputSlot, TList::Class()); fOutputSlotHistQA=outputSlot++;}  // Qvector QA histograms
  if(fProvideQnVectorsList)                                    {DefineOutput(outputSlot, TList::Class()); fOutputSlotQnVectorsList=outputSlot++;}   // Calibrated qvector list
  if(fFillEventQA)                                             {DefineOutput(outputSlot, TList::Class()); fOutputSlotEventQA=outputSlot++;}   // Event QA histograms


  TString paths(gROOT->GetMacroPath());
  TObjArray* patharray = fRunListPath.Tokenize("/");
  TString appendpath="";
  for(Int_t i=0; i<patharray->GetEntries()-1; i++) appendpath=appendpath+patharray->At(i)->GetName()+"/";
  paths.Append(appendpath);

  gROOT->SetMacroPath(paths);
  const char* runlistfile = gSystem->Which(gROOT->GetMacroPath(), fRunListPath.Data());
  fEventPlaneManager->SetFileRunLabels(runlistfile);

}


//_________________________________________________________________________________
void AliAnalysisTaskFlowVectorCorrections::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //

  fFillEvent->SetEventPlaneManager(fEventPlaneManager);

  TFile* inputfile = 0x0;
  if(!(fCalibrationFilePath.EqualTo(""))) {
    if(fCalibrationFilePath.Contains("alien")) TGrid::Connect("alien://");
    inputfile = TFile::Open(fCalibrationFilePath);
  }
  if(inputfile){
    fEventPlaneManager->SetCalibrationFile(inputfile);
    inputfile->Close();
  }

  fEventPlaneManager->Initialize();

  if(fEventPlaneManager->ShouldFillHistogramsQnCorrections()) PostData(fOutputSlotHistQn, fEventPlaneManager->GetListOutputHistogramsQnCorrections());
  if(fEventPlaneManager->ShouldFillTreeQnVectors())   PostData(fOutputSlotTree, fEventPlaneManager->GetTreeQnVectors());
  if(fEventPlaneManager->ShouldFillHistogramsQA()) PostData(fOutputSlotHistQA, fEventPlaneManager->GetListHistogramsQA());
  if(fFillEventQA)   PostData(fOutputSlotEventQA, fEventQAList);

}



//________________________________________________________________________________________________________
void AliAnalysisTaskFlowVectorCorrections::UserExec(Option_t *){
  //
  // Main loop. Called for every event
  //

  AliVEvent* event = InputEvent();
  //AliReducedEvent* event = dynamic_cast<AliReducedEvent*>(GetInputData(0));

  fEventPlaneManager->ClearEvent();

  Float_t* values = fEventPlaneManager->GetDataContainer();


  fFillEvent->Process((AliAnalysisTaskSE*) this, event, values);
  //fFillEvent->Process(event, values);

  if(fCalibrateByRun) fEventPlaneManager->SetCalibrationFileDirectoryName(Form("%d",(Int_t) (values[AliQnCorrectionsVarManager::kRunNo]+10E-6)));

  //if(!fInitialized){ fEventPlaneManager->InitializeCalibrationHistograms(); fInitialized = kTRUE;}



  fEventPlaneHistos->FillHistClass("Event_NoCuts", values);
  //for(UShort_t ibit=0; ibit<64; ++ibit) {
  //HIST::Instance()->FillHistClass("OfflineTriggers_NoCuts", values);
  //}


  // use only selected triggers for event plane calibration averages
  if(IsEventSelected(values)){//&&TriggerSelected(event)) {
    fEventPlaneHistos->FillHistClass("Event_Analysis", values);
    //for(UShort_t ibit=0; ibit<64; ++ibit) {
    //VAR::FillEventOfflineTriggers(ibit, event, values);
    //HIST::Instance()->FillHistClass("OfflineTriggers_WithCuts", values);
    //}

    fEventPlaneManager->Process();



  }  // end if event selection

  //fQvectorList.Print();
  //gObjectTable->Print();
  if(fProvideQnVectorsList) PostData(fOutputSlotQnVectorsList, fEventPlaneManager->GetListQnVectors());

  }  // end loop over events


  //__________________________________________________________________
  void AliAnalysisTaskFlowVectorCorrections::FinishTaskOutput()
  {
    //
    // Finish Task 
    //

    fEventPlaneManager->Finalize();
    //fEventPlaneManager->WriteCalibrationHistogramsToList();
    //fEventPlaneManager->WriteQaHistogramsToList();// fEventPlaneHistos->HistList());


    THashList* hList = (THashList*) fEventPlaneHistos->HistList();
    for(Int_t i=0; i<hList->GetEntries(); ++i) {
      THashList* list = (THashList*)hList->At(i);
      fEventQAList->Add(list);
    }

  }



  //__________________________________________________________________
  Bool_t AliAnalysisTaskFlowVectorCorrections::IsEventSelected(Float_t* values) {
    if(!fEventCuts) return kTRUE;
    return fEventCuts->IsSelected(values);
  }


