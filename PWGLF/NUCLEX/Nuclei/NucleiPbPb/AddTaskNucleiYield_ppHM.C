/// \file AddTaskNucleiYield.C
/// \brief Simple macro to add the task to a grid job
///
/// The task is here added several times to analyse different particle species and to investigate
/// different set of cuts in only one run.
///
/// \author Maximiliano Puccio <maximiliano.puccio@cern.ch>, University and INFN Torino
/// \date Feb 19th, 2015

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskNucleiYield.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliPID.h"
#endif

AliAnalysisTaskNucleiYield* AddTaskNucleiYield_ppHM(Bool_t isMC = kFALSE,
    AliPID::EParticleType part = AliPID::kDeuteron,
    Int_t pdgCode = 1000010020,
    TString tskname = "deuteron",
    TString suffix = "") {

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNucleiYield", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskNucleiYield", "This task requires an input event handler");
    return 0x0;
  }

  tskname.Append(Form("%s",suffix.Data()));

  AliAnalysisTaskNucleiYield *deu = new AliAnalysisTaskNucleiYield(tskname);
  deu->fEventCut.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0 | AliVEvent::kINT7);

  deu->SetParticleType(part);
  deu->SetPDG(pdgCode);
  deu->SetIsMC(isMC);
  deu->SetDCABins(80,-0.5,0.5);

  deu->SetRequireTPCpidSigmas(3.f);

  Float_t lDesiredBoundaries[1000];
  Long_t   lNDesiredBoundaries=0;
  lDesiredBoundaries[0] = 0.0;
  //From High To Low Multiplicity
  for( Int_t ib = 1; ib < 101; ib++) { // 100 bins  ] 0.0 , 0.1 ]
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.001;
  }
  for( Int_t ib = 1; ib < 91; ib++) { // 90 bins  ] 0.1 , 1.0 ]
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.01;
  }
  for( Int_t ib = 1; ib < 91; ib++) { // 90 bins ] 1.0 , 10. ]
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.1;
  }
  for( Int_t ib = 1; ib < 96; ib++) { // 95 bins ] 10.0 , 105.0 ]
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 1.0;
  }
  deu->SetCentBins(lNDesiredBoundaries, lDesiredBoundaries);
  deu->SetUseFlattening(false);
  float pt[20] = {
    0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,1.6f,
    1.8f,2.0f,2.2f,2.6f,3.0f,3.4f,3.8f,4.4f,5.0f,6.0f
  };
  deu->SetPtBins(19,pt);

  float dcabins[53] = {
    -1.30,-1.20,-1.10,-1.00,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,
    -0.35,-0.30,-0.25,-0.20,-0.15,-0.12,-0.10,-0.09,-0.08,-0.07,
    -0.06,-0.05,-0.04,-0.03,-0.02,-0.01, 0.00, 0.01, 0.02, 0.03,
     0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.15, 0.20,
     0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00,
     1.10, 1.20, 1.30
  };
  deu->SetDCABins(52,dcabins);

  mgr->AddTask(deu);

  TString output = "AnalysisResults.root";
  AliAnalysisDataContainer *deuCont = mgr->CreateContainer(Form("nuclei_%s",tskname.Data()),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      output.Data());
  AliAnalysisDataContainer *deuCont1 = mgr->CreateContainer(Form("RTree%s",tskname.Data()),
      TTree::Class(),
      AliAnalysisManager::kOutputContainer,
      output.Data());
  AliAnalysisDataContainer *deuCont2 = mgr->CreateContainer(Form("STree%s",tskname.Data()),
      TTree::Class(),
      AliAnalysisManager::kOutputContainer,
      output.Data());
  mgr->ConnectInput  (deu,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu,  1, deuCont);
  mgr->ConnectOutput (deu,  2, deuCont1);
  mgr->ConnectOutput (deu,  3, deuCont2);
  return deu;
}
