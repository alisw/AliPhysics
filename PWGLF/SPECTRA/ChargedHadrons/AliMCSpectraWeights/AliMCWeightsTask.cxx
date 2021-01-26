//
//  AliMCWeightsTask.cpp
//  ParticleComposition
//
//  Created by Patrick Huhn on 05.10.20.
//  Copyright © 2020 Patrick Huhn. All rights reserved.
//
#include "AliMCWeightsTask.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliPhysicsSelection.h"
#include "AliVEvent.h"
#include "AliMultSelection.h"
#include "AliPhysicsSelection.h"
#include "TChain.h"
#include "TList.h"
#include "TString.h"
#include <iostream>

class AliMCWeightsTask;

/// \cond CLASSIMP
//ClassImp(AliMCWeightsTask);
/// \endcond

AliMCWeightsTask::AliMCWeightsTask()
: AliAnalysisTaskSE(), fOutputList(nullptr),
fEvent(nullptr), fMCEvent(nullptr), fMCSpectraWeights(nullptr) {
    // default constructor for root
}

AliMCWeightsTask::AliMCWeightsTask(const char* name)
: AliAnalysisTaskSE(name), fOutputList(nullptr),
fEvent(nullptr), fMCEvent(nullptr), fMCSpectraWeights(nullptr) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliMCWeightsTask::~AliMCWeightsTask() {
    // destructor
    if (fOutputList) delete fOutputList;
}

void AliMCWeightsTask::UserCreateOutputObjects() {
    OpenFile(1, "recreate");
    fOutputList = new TList();
    fOutputList->SetOwner();

    fOutputList->Add(
                     (TObject*)fMCSpectraWeights->GetHistMCGenPrimTrackParticles());
    fOutputList->Add((TObject*)fMCSpectraWeights->GetHistDataFraction());
    fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCFraction());
    fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCWeights());
    //    fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCWeightsSys());
    //    //TODO: check if this works

    PostData(1, fOutputList);
}

void AliMCWeightsTask::UserExec(Option_t* option) {
    AliInputEventHandler* inputHandler =
    (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()
    ->GetInputEventHandler();
    if (!inputHandler) {
        Printf("ERROR: Could not receive inputHandler");
        return;
    }

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

    if (fMCSpectraWeights->GetTaskStatus() <
        AliMCSpectraWeights::TaskState::kMCSpectraObtained)
        fMCSpectraWeights->FillMCSpectra(fMCEvent);
    else {
        fMCSpectraWeights->SetCurrentEvent(fMCEvent);
        fMCSpectraWeights->StartNewEvent();

        TString fStoredObjectName = "fMCSpectraWeights";
        // Add to AliVEvent
        if ((!(fEvent->FindListObject(fStoredObjectName.Data())))) {
            fEvent->AddObject(fMCSpectraWeights);
        }
    }

    PostData(1, fOutputList);
}

AliMCWeightsTask*
AliMCWeightsTask::AddTaskAliMCWeightsTask() {

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
    fileName += "AliMCWeightsTask"; // create a subfolder in the file

    // create AliMCSpectraWeights
    //===========================================================================
    std::string stTrainOutputPath;
    std::string collisionSystem;
//    switch (genType) { // TODO: fill path here
//        case MCGeneratorType::PP_PYTHIA:
            stTrainOutputPath = "~/particle-composition-correction/data/train/"
            "pp_PCC_LHC17l3b_CENT_wSDD.root";
            collisionSystem = "pp";
//            break;
//        case MCGeneratorType::PPB_EPOS:
//            stTrainOutputPath = "";
//            collisionSystem = "ppb";
//            break;
//        case MCGeneratorType::PBPB_HIJING:
//            stTrainOutputPath = "";
//            collisionSystem = "pbpb";
//            break;
//        default:
//            break;
//    }

    AliMCSpectraWeights* fMCSpectraWeights =
    new AliMCSpectraWeights(collisionSystem, "fMCSpectraWeights",
                            AliMCSpectraWeights::SysFlag::kNominal);
    if (stTrainOutputPath.size() > 5) {
        fMCSpectraWeights->SetMCSpectraFile(stTrainOutputPath);
        fMCSpectraWeights->SetDoSystematics(true);
    }
    fMCSpectraWeights->Init();

    // Configure task
    //===========================================================================
    std::cout << "Configure task\n";
    AliMCWeightsTask* task = new AliMCWeightsTask("AliMCWeightsTask");
    task->SetMCSpectraWeightObject(fMCSpectraWeights);

    // attach the task to the manager and configure in and ouput
    //===========================================================================
    std::cout << "attach the task to the manager and configure in and ouput\n";
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(
                       task, 1,
                       mgr->CreateContainer("AliMCWeightsTask", TList::Class(),
                                            AliAnalysisManager::kOutputContainer,
                                            fileName.Data()));

    return task;
}
