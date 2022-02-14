#include <algorithm>
#include <iostream>
#include <vector>
#include "THistManager.h"
#include <TList.h>
#include <TMath.h>
#include <TString.h>
#include "AliEMCALTriggerPatchInfo.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalTriggerBackground.h"
#include "AliVCaloTrigger.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerBackground)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalTriggerBackground::AliAnalysisTaskEmcalTriggerBackground():
    AliAnalysisTaskEmcal(),
    fHistos(nullptr)
{
    SetCaloTriggerPatchInfoName("EmcalTriggers");
    SetCaloTriggersName("EMCALTrigger");
}

AliAnalysisTaskEmcalTriggerBackground::AliAnalysisTaskEmcalTriggerBackground(const char *name) :
    AliAnalysisTaskEmcal(name, true),
    fHistos(nullptr)
{
    SetCaloTriggerPatchInfoName("EmcalTriggers");
    SetCaloTriggersName("EMCALTrigger");
    SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalTriggerBackground::~AliAnalysisTaskEmcalTriggerBackground() {
    if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalTriggerBackground::UserCreateOutputObjects() {
    AliAnalysisTaskEmcal::UserCreateOutputObjects();

    fHistos = new THistManager("backgroundtaskhists");
    fHistos->CreateTH2("hRhoEMCAL", "Rho values in EMCAL", 1000, 0., 1000, 1000, 0., 1000.);
    fHistos->CreateTH2("hRhoEMCALMedianDCAL", "Rho value in EMCAL vs Median DCAL", 1000, 0., 1000, 1000, 0., 1000.);
    fHistos->CreateTH2("hRhoDCAL", "Rho values in DCAL", 1000, 0., 1000, 1000, 0., 1000.);
    fHistos->CreateTH2("hRhoDCALMedianEMCAL", "Rho value in DCAL vs Median EMCAL", 1000, 0., 1000, 1000, 0., 1000.);

    for(auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
}

bool AliAnalysisTaskEmcalTriggerBackground::Run() {
    std::vector<double> bkgpatchesEMCAL, bkgpatchesDCAL;
    
    for(auto p : *fTriggerPatchInfo) {
        AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(p);
        if(!patch->IsBkgRecalc()) continue;
        if(patch->IsEMCal()) bkgpatchesEMCAL.push_back(patch->GetADCAmp());
        else bkgpatchesDCAL.push_back(patch->GetADCAmp());
    }

    std::sort(bkgpatchesEMCAL.begin(), bkgpatchesEMCAL.end(), std::greater<double>());
    std::sort(bkgpatchesDCAL.begin(), bkgpatchesDCAL.end(), std::greater<double>());

    auto medEMCAL = TMath::Median<double>(bkgpatchesEMCAL.size(), bkgpatchesEMCAL.data()),
         medDCAL = TMath::Median<double>(bkgpatchesDCAL.size(), bkgpatchesDCAL.data());

    fHistos->FillTH2("hRhoEMCAL", fCaloTriggers->GetMedian(0), medEMCAL);
    fHistos->FillTH2("hRhoDCAL", fCaloTriggers->GetMedian(1), medDCAL);
    fHistos->FillTH2("hRhoEMCALMedianDCAL", fCaloTriggers->GetMedian(0), medDCAL);
    fHistos->FillTH2("hRhoDCALMedianEMCAL", fCaloTriggers->GetMedian(1), medEMCAL);
    return true;
}

AliAnalysisTaskEmcalTriggerBackground * AliAnalysisTaskEmcalTriggerBackground::AddTaskEmcalTriggerBackground(const char *name){
    auto mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr) {
        std::cerr << "No Analysis manager available." << std::endl;
        return nullptr;
    }

    AliAnalysisTaskEmcalTriggerBackground *task = new AliAnalysisTaskEmcalTriggerBackground(name);
    mgr->AddTask(task);

    TString outfilename(mgr->GetCommonFileName());
    outfilename += Form(":EmcalTriggerBackground%s", name);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("EmcalTriggerBackground%s", name), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename.Data()));
    std::cout << "Done" << std::endl;
    return task;
}