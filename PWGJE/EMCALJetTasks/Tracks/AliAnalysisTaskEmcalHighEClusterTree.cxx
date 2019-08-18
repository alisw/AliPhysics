/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <TLorentzVector.h>
#include <TNtuple.h>

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalHighEClusterTree.h"
#include "AliClusterContainer.h"
#include "AliLog.h"
#include "AliVCaloCells.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalHighEClusterTree)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalHighEClusterTree::AliAnalysisTaskEmcalHighEClusterTree():
    AliAnalysisTaskEmcal(),
    fOutputTree(nullptr),
    fMinClusterE(150.)
{

}

AliAnalysisTaskEmcalHighEClusterTree::AliAnalysisTaskEmcalHighEClusterTree(const char *name):
    AliAnalysisTaskEmcal(name, kTRUE),
    fOutputTree(nullptr),
    fMinClusterE(150.)
{
    DefineOutput(2, TNtuple::Class());
    SetMakeGeneralHistograms(kTRUE);
}

AliAnalysisTaskEmcalHighEClusterTree::~AliAnalysisTaskEmcalHighEClusterTree() {
    
}

void AliAnalysisTaskEmcalHighEClusterTree::UserCreateOutputObjects(){
    AliAnalysisTaskEmcal::UserCreateOutputObjects();
    fOutputTree = new TNtuple("HighEClusters", "Tree for high-energy clusters", "Eraw:Enonlin:Eleading:Exotic:Lambda02:Lambda20:Ncell:Time:Eta:Phi:IsMinBias:IsEJ1:IsEJ2:IsDJ1:IsDJ2");   
    PostData(1, fOutput);
    PostData(2, fOutputTree);
}

bool AliAnalysisTaskEmcalHighEClusterTree::Run(){
    auto clustercont = GetClusterContainer("clustercontainer");
    float IsMinBias = (fInputHandler->IsEventSelected() & AliVEvent::kINT7) ? 1. : 0.,
          IsEJ1 = ((fInputHandler->IsEventSelected() & AliVEvent::kINT7) && fInputEvent->GetFiredTriggerClasses().Contains("EJ1")) ? 1. : 0.,
          IsEJ2 = ((fInputHandler->IsEventSelected() & AliVEvent::kINT7) && fInputEvent->GetFiredTriggerClasses().Contains("EJ2")) ? 1. : 0.,
          IsDJ1 = ((fInputHandler->IsEventSelected() & AliVEvent::kINT7) && fInputEvent->GetFiredTriggerClasses().Contains("DJ1")) ? 1. : 0.,
          IsDJ2 = ((fInputHandler->IsEventSelected() & AliVEvent::kINT7) && fInputEvent->GetFiredTriggerClasses().Contains("DJ2")) ? 1. : 0.;
    for(auto c : clustercont->all()) {
        if(c->E() < fMinClusterE) continue;

        TLorentzVector pvec;
        c->GetMomentum(pvec, fVertex);

        float datapoint[15] = {static_cast<float>(c->E()), static_cast<float>(c->GetNonLinCorrEnergy()), 0., static_cast<float>(c->GetIsExotic() ? 1. : 0.), 
                               static_cast<float>(c->GetM02()), static_cast<float>(c->GetM20()), static_cast<float>(c->GetNCells()), static_cast<float>(c->GetTOF()*1e9),
                               static_cast<float>(pvec.Eta()), static_cast<float>(pvec.Phi()), IsMinBias, IsEJ1, IsEJ2, IsDJ1, IsDJ2 };
        float emaxcell(0.);
        for(auto icell = 0; icell < c->GetNCells(); icell++){
            float celltmp = fInputEvent->GetEMCALCells()->GetCellAmplitude(c->GetCellAbsId(icell));
            if(celltmp > emaxcell) emaxcell = celltmp;
        }
        datapoint[2] = emaxcell; 
        fOutputTree->Fill(datapoint);
    }   
}

bool AliAnalysisTaskEmcalHighEClusterTree::IsTriggerSelected() {
    return ((fInputHandler->IsEventSelected() && AliVEvent::kINT7) || (fInputHandler->IsEventSelected() && AliVEvent::kEMCEJE));
}

AliAnalysisTaskEmcalHighEClusterTree *AliAnalysisTaskEmcalHighEClusterTree::AddTaskEmcalHighClusterE(const char *name) {
    auto mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr) {
        AliErrorGeneralStream("AliAnalysisTaskEmcalHighEClusterTree::AddTaskEmcalHighClusterE") << "No analysis manager available" << std::endl;
        return nullptr; 
    }    

    bool isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();

    auto task = new AliAnalysisTaskEmcalHighEClusterTree(name);
    mgr->AddTask(task);

    auto clustercont = task->AddClusterContainer("usedefault");
    clustercont->SetExoticCut(false);
    clustercont->SetName("clustercontainer");

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer("highEClusterHistograms", TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));
    mgr->ConnectOutput(task, 2, mgr->CreateContainer("highEClusterTree", TNtuple::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));

    return task;
}