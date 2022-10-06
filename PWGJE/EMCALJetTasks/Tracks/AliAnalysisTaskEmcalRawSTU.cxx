/************************************************************************************
 * Copyright (C) 2021, Copyright Holders of the ALICE Collaboration                 *
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
#include <bitset>
#include <map>

#include <THistManager.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalRawSTU.h"
#include "AliEMCALTriggerTypes.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliVCaloTrigger.h"
#include "AliVEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalRawSTU)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalRawSTU::AliAnalysisTaskEmcalRawSTU():
    AliAnalysisTaskEmcal(),
    fHistos(nullptr)
{

}

AliAnalysisTaskEmcalRawSTU::AliAnalysisTaskEmcalRawSTU(const char *name):
    AliAnalysisTaskEmcal(name, true),
    fHistos(nullptr)
{
    SetCaloTriggersName("usedefault");
}

bool AliAnalysisTaskEmcalRawSTU::IsTriggerSelected() {
    return (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA);
}

void AliAnalysisTaskEmcalRawSTU::UserCreateOutputObjects(){
    AliAnalysisTaskEmcal::UserCreateOutputObjects();

    fHistos = new THistManager("STUhistos");
    std::vector<std::string> triggers = {"EG1", "EG2", "DG1", "DG2", "EDG1", "EDG2", "EG1noDG1", "EG1andDG1", "EG2noDG2", "EG2andDG2", "DG1noEG1", "DG2noEG2"};
    for(const auto &trg : triggers) {
        fHistos->CreateTH2(Form("hRowCol%s", trg.data()), Form("Patch row-col %s", trg.data()), 48, -0.5, 47.5, 104, -0.5, 103.5);
    }

    for(auto hist : *fHistos->GetListOfHistograms()) fOutput->Add(hist);
    PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalRawSTU::Run(){
    AliDebugStream(1) << "Start new event" << std::endl;
    auto triggers = fInputEvent->GetFiredTriggerClasses();
    bool isEG1 = (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && triggers.Contains("EG1"),
         isEG2 = (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && triggers.Contains("EG2"),
         isDG1 = (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && triggers.Contains("DG1"),
         isDG2 = (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) && triggers.Contains("DG2");
    if(!(isEG1 || isEG2 || isDG1 || isDG2)) return false;
    AliDebugStream(1) << "Event was selected" << std::endl;

    fCaloTriggers->Reset();
    int col, row, triggerbitsRaw;
    int nall(0), nsel(0);
    while(fCaloTriggers->Next()) {
        nall++;
        fCaloTriggers->GetTriggerBits(triggerbitsRaw);
        if(!triggerbitsRaw) continue;
        nsel++;
        std::bitset<sizeof(triggerbitsRaw) * 8> triggerbits(triggerbitsRaw);
        fCaloTriggers->GetPosition(col, row);

        if(triggerbits.test(kL1GammaHigh) || (triggerbits.test(kL1GammaHigh + kTriggerTypeEnd))) {
            AliDebugStream(2) << "Found patch High" << std::endl;
            // Gamma high
            bool filledCombined = false;
            if(isEG1) {
                fHistos->FillTH1("hRowColEG1", col, row);
                fHistos->FillTH1("hRowColEDG1", col, row);
                filledCombined = true;
                if(isDG1) fHistos->FillTH1("hRowColEG1andDG1", col, row);
                else fHistos->FillTH1("hRowColEG1noDG1", col, row);
            }
            if(isDG1) {
                fHistos->FillTH1("hRowColDG1", col, row);
                if(!filledCombined) fHistos->FillTH1("hRowColEDG1", col, row);
                if(!isEG1) fHistos->FillTH1("hRowColDG1noEG1", col, row);
            }
        }
        if(triggerbits.test(kL1GammaLow) || (triggerbits.test(kL1GammaLow + kTriggerTypeEnd))) {
            AliDebugStream(2) << "Found patch Low" << std::endl;
            // Gamma low
            bool filledCombined = false;
            if(isEG2) {
                fHistos->FillTH1("hRowColEG2", col, row);
                fHistos->FillTH1("hRowColEDG2", col, row);
                filledCombined = true;
                if(isDG2) fHistos->FillTH1("hRowColEG2andDG2", col, row);
                else fHistos->FillTH1("hRowColEG2noDG2", col, row);
            }
            if(isDG2) {
                fHistos->FillTH1("hRowColDG2", col, row);
                if(!filledCombined) fHistos->FillTH1("hRowColEDG2", col, row);
                if(!isEG2) fHistos->FillTH1("hRowColDG2noEG2", col, row);
            }
        }
    }
    AliDebugStream(1) << "All: " << nall << ", sel " << nsel << std::endl;
    
    return true;
}

AliAnalysisTaskEmcalRawSTU *AliAnalysisTaskEmcalRawSTU::AddTaskEmcalRawSTU(const char *name) {
   auto mgr = AliAnalysisManager::GetAnalysisManager(); 

   if(!mgr) {
        AliFatalGeneral("AliAnalysisTaskEmcalRawSTU::AddTaskEmcalRawSTU", "No analysis manager created");
   }

   auto task = new AliAnalysisTaskEmcalRawSTU(name);
   mgr->AddTask(task);

   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, mgr->CreateContainer("RawSTUHistos", TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));
   return task;
}