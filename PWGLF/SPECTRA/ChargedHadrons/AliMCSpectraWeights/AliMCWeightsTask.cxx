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
#include <chrono>

#ifdef __AliMCWeightsTask_DebugPCC__
#define DebugPCC(x) std::cout << x
#else
#define DebugPCC(x)
#endif

#ifdef __AliMCWeightsTask_DebugTiming__
#define DebugChrono(x) std::cout << x
#else
#define DebugChrono(x)
#endif

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
#ifdef __AliMCWeightsTask_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    // destructor
    if (fOutputList) delete fOutputList;

#ifdef __AliMCWeightsTask_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    DebugChrono("deletion took " << duration << " microseconds\n");
#endif
}

void AliMCWeightsTask::UserCreateOutputObjects() {
#ifdef __AliMCWeightsTask_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    OpenFile(1, "recreate");
    fOutputList = new TList();
    fOutputList->SetOwner();

    fOutputList->Add(
                     (TObject*)fMCSpectraWeights->GetHistMCGenPrimTrackParticles());
    fOutputList->Add((TObject*)fMCSpectraWeights->GetHistDataFraction());
    fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCFraction());
    fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCWeights());
    //    fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCWeightsSys());

    PostData(1, fOutputList);
#ifdef __AliMCWeightsTask_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    DebugChrono("UserCreateOutputObjects took " << duration << " microseconds\n");
#endif
}

void AliMCWeightsTask::UserExec(Option_t* option) {
#ifdef __AliMCWeightsTask_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
    DebugChrono("Start new AliMCWeightsTask::UserExec\n");
#endif

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

    fEvent = InputEvent();
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
        auto tmpObject = static_cast<AliMCSpectraWeightsHandler*>(fEvent->FindListObject(fStoredObjectName.Data()));
        if (!tmpObject) {
            AliMCSpectraWeightsHandler* handler = new AliMCSpectraWeightsHandler(fMCSpectraWeights, fStoredObjectName.Data());
            fEvent->AddObject(handler);
            DebugPCC("Added fMCSpectraWeights to event\n");
        }
        else{
            DebugPCC("fMCSpectraWeights already in event\n");
        }

    }

    PostData(1, fOutputList);
#ifdef __AliMCWeightsTask_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    DebugChrono("UserExec took " << duration << " microseconds\n");
#endif
}

AliMCWeightsTask*
AliMCWeightsTask::AddTaskAliMCWeightsTask() {
#ifdef __AliMCWeightsTask_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif

    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        printf("AliMCWeightsTask::No analysis manager to connect to.");
        return nullptr;
    }

    // Check the analysis type using the event handlers connected to the
    // analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        printf("AliMCWeightsTask::This task requires an input event handler");
        return nullptr;
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
//            stTrainOutputPath = "~/particle-composition-correction/data/train/"
//            "pp_PCC_LHC17l3b_CENT_wSDD.root";
            stTrainOutputPath = "alien:///alice/cern.ch/user/p/phuhn/"
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

#ifdef __AliMCWeightsTask_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    DebugChrono("AliMCWeightsTask took " << duration << " microseconds\n");
#endif
    
    return task;
}
