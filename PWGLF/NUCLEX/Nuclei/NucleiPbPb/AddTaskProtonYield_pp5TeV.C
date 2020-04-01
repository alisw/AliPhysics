#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskNucleiYield.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliPID.h"
#endif

AliAnalysisTaskNucleiYield* AddTaskProtonYield_pp5TeV(Bool_t isMC = kFALSE,
    AliPID::EParticleType part = AliPID::kProton,
    Int_t pdgCode = 2212,
    TString tskname = "proton",
    bool saveTrees = true,
    TString suffix = "") {

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskProtonYield_pp5TeV", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskProtonYield_pp5TeV", "This task requires an input event handler");
    return 0x0;
  }

  tskname.Append(Form("%s",suffix.Data()));

  AliAnalysisTaskNucleiYield *task = new AliAnalysisTaskNucleiYield(tskname);

  task->SetParticleType(part);
  task->SetPDG(pdgCode);
  task->SetIsMC(isMC);
  task->SetDCABins(80,-0.5,0.5);
  if (saveTrees)
	  task->SaveTrees();

  task->SetRequireTPCpidSigmas(3.f);
  float cent[14] = {-5.f,0.f,1.f,5.f,10.f,20.f,30.f,40.f,50.f,60.f,70.f,80.f,90.f,100.f};
  task->SetCentBins(13, cent);
  task->SetUseFlattening(false);
  float pt[20] = {
    0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,1.6f,
    1.8f,2.0f,2.2f,2.6f,3.0f,3.4f,3.8f,4.4f,5.0f,6.0f
  };
  task->SetPtBins(19,pt);

  float dcabins[53] = {
    -1.30,-1.20,-1.10,-1.00,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,
    -0.35,-0.30,-0.25,-0.20,-0.15,-0.12,-0.10,-0.09,-0.08,-0.07,
    -0.06,-0.05,-0.04,-0.03,-0.02,-0.01, 0.00, 0.01, 0.02, 0.03,
     0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.15, 0.20,
     0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00,
     1.10, 1.20, 1.30
  };
  task->SetDCABins(52,dcabins);

  mgr->AddTask(task);

  int slot = 0;
  TString output = "AnalysisResults.root";
  AliAnalysisDataContainer *cont = mgr->CreateContainer(Form("nuclei_%s",tskname.Data()),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      output.Data());
  mgr->ConnectInput  (task,  slot++, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  slot++, cont);
  if (saveTrees) {
    cont = mgr->CreateContainer(Form("%s_rectree",tskname.Data()),
        TTree::Class(),
        AliAnalysisManager::kOutputContainer,
        output.Data());
    mgr->ConnectOutput (task,  slot++, cont);
    if (isMC) {
      cont = mgr->CreateContainer(Form("%s_simtree",tskname.Data()),
        TTree::Class(),
        AliAnalysisManager::kOutputContainer,
        output.Data());
      mgr->ConnectOutput (task,  slot++, cont);
    }
  }

  return task;
}
