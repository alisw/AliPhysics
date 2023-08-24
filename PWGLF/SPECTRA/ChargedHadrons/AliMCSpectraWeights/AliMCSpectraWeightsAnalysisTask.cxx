#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEventCuts.h"
#include "AliMCEvent.h"
#include "AliMCSpectraWeights.h"
#include "AliMultSelection.h"
#include "AliPhysicsSelection.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVMultiplicity.h"
#include "TArrayD.h"
#include "TChain.h"
#include "TList.h"
#include "TString.h"
#include <iostream>
#include <string>

#include "AliMCSpectraWeightsAnalysisTask.h"

class AliMCSpectraWeightsAnalysisTask;

ClassImp(AliMCSpectraWeightsAnalysisTask);

AliMCSpectraWeightsAnalysisTask::AliMCSpectraWeightsAnalysisTask()
    : AliAnalysisTaskSE(), // default constructor
      fDebugLevel(0), fOutputList(0), fEvent(0), fMCEvent(0), fMCStack(0),
      fEventCuts(0), fIsMC(1), fTriggerMask(AliVEvent::kMB),
      fMCSpectraWeights(0) {
    // this is used by root for IO purposes, it needs to remain empty
}

AliMCSpectraWeightsAnalysisTask::AliMCSpectraWeightsAnalysisTask(
    const char* name)
    : AliAnalysisTaskSE(name), // default constructor
    fDebugLevel(0), fOutputList(0), fEvent(0), fMCEvent(0), fMCStack(0),
    fEventCuts(0), fIsMC(1), fTriggerMask(AliVEvent::kMB),
    fMCSpectraWeights(0) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliMCSpectraWeightsAnalysisTask::~AliMCSpectraWeightsAnalysisTask() {
    if (fOutputList)
        delete fOutputList;
}

void AliMCSpectraWeightsAnalysisTask::UserCreateOutputObjects() {
    OpenFile(1, "recreate");
    fOutputList = new TList();
    fOutputList->SetOwner();
    if (fIsMC && fMCSpectraWeights) {
        fOutputList->Add(
            (TObject*)fMCSpectraWeights->GetHistMCGenPrimTrackParticles());
        fOutputList->Add((TObject*)fMCSpectraWeights->GetHistDataFraction());
        fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCFraction());
        fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCWeights());
    } else
        printf("AliMCSpectraWeightsAnalysisTask:: Either running not MC or "
               "object of AliMCSpectraWeights is null pointer\n");

    PostData(1, fOutputList);
}

void AliMCSpectraWeightsAnalysisTask::UserExec(Option_t* option) {
    if (!fIsMC) {
        printf("AliMCSpectraWeightsAnalysisTask:: WARNING task only works for "
               "MC!!\n");
        return;
    }
    AliInputEventHandler* inputHandler =
        (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()
            ->GetInputEventHandler();
    if (!inputHandler) {
        Printf("ERROR: Could not receive inputHandler");
        return;
    }
    // check if trigger mask is fine
    if (!(inputHandler->IsEventSelected() & fTriggerMask))
        return;
    AliPhysicsSelection* physicsSelection =
        static_cast<AliPhysicsSelection*>(inputHandler->GetEventSelection());
    if (!physicsSelection) {
        Printf("ERROR: Could not receive physicsSelection");
        return;
    }
    fEvent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fEvent) {
        printf("ERROR: fEvent not available\n");
        return;
    }

    fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!fMCEvent) {
        printf("ERROR: fMCEvent not available\n");
        return;
    }
    fMCStack = fMCEvent->Stack();
    if (!fMCStack) {
        printf("ERROR: fMCStack not available\n");
        return;
    }

    if (!fEventCuts.AcceptEvent(fEvent)) {
        if(fDebugLevel > 0)
        {
            printf("LOG: event not passed selection\n");
        }
        return;
    }
    // fill MC spectra
    //==============================================================================
    if (fMCSpectraWeights->GetTaskStatus() <
        AliMCSpectraWeights::TaskState::kMCSpectraObtained)
        fMCSpectraWeights->FillMCSpectra(fMCEvent);

    PostData(1, fOutputList);
}

AliMCSpectraWeightsAnalysisTask*
AliMCSpectraWeightsAnalysisTask::AddTaskWeights(const char* collisionSystem,
                                                const char* previousTrain,
                                                const char* name,
                                                const char* outfile) {

    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        printf("AddTaskWeights::No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the
    // analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        printf("AddTaskWeights::This task requires an input event handler");
        return NULL;
    }

    // Setup output file
    //===========================================================================
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":";
    fileName += name; // create a subfolder in the file
    if (outfile) {    // if a finename is given, use that one
        fileName = TString(outfile);
    }

    // create AliMCSpectraWeights
    //===========================================================================
    std::cout << "create AliMCSpectraWeights\n";
    std::string stTrainOutputPath{previousTrain};
    AliMCSpectraWeights* fMCSpectraWeights =
        new AliMCSpectraWeights(collisionSystem, "fMCSpectraWeights",
                                AliMCSpectraWeights::SysFlag::kNominal);
    if (stTrainOutputPath.size() > 5)
        fMCSpectraWeights->SetMCSpectraFile(previousTrain);
    fMCSpectraWeights->Init();

    // Configure task
    //===========================================================================
    std::cout << "Configure task\n";
    AliMCSpectraWeightsAnalysisTask* task =
        new AliMCSpectraWeightsAnalysisTask(name);
    task->SetUseMC(
        (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() !=
         0x0));
    task->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kMB); // kINT7
    // For particle composition
    task->SetMCSpectraWeightObject(fMCSpectraWeights);
    // Debug info
    task->SetDebugLevel(0);

    // attach the task to the manager and configure in and ouput
    //===========================================================================
    std::cout << "attach the task to the manager and configure in and ouput\n";
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(
        task, 1,
        mgr->CreateContainer(name, TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             fileName.Data()));

    return task;
}
