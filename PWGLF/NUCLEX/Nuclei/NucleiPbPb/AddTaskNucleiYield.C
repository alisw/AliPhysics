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

AliAnalysisTaskNucleiYield* AddTaskNucleiYield(Bool_t isMC = kFALSE,
                                               AliPID::EParticleType part = AliPID::kDeuteron,
                                               Int_t pdgCode = 1000010020,
                                               TString tskname = "",
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

  // Common variables
  Float_t tritonCentBins[7] = {0.f,5.f,10.f,20.f,30.f,40.f,50.f};
  Float_t deuteronPtBins[27] = {
    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
    3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f
  };
  Float_t deutBBpar[5] = {4.69637f,7.51827f,0.0183746f,2.60f,2.7f};
  Float_t sigmaBBpar = 0.1f;
  Float_t dcabins[32] = {
    -0.50,-0.40,-0.30,-0.20,-0.15,-0.12,-0.10,-0.09,-0.08,-0.07,
    -0.06,-0.05,-0.04,-0.03,-0.02,-0.01, 0.00, 0.01, 0.02, 0.03,
    0.04, 0.05, 0.06, 0.07, 0.09, 0.10, 0.12, 0.15, 0.20, 0.30,
    0.40, 0.50
  };

  tskname.Append(Form("%s",suffix.Data()));
  AliAnalysisTaskNucleiYield *deu = new AliAnalysisTaskNucleiYield(tskname);
  deu->SetParticleType(part);
  deu->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu->SetPDG(pdgCode);
  deu->SetIsMC(isMC);
  deu->SetCentBins(6, tritonCentBins);
  deu->SetPtBins(26,deuteronPtBins);
  deu->SetDCABins(80,-0.5,0.5);
  deu->SetDCABins(31,dcabins);
  mgr->AddTask(deu);
  TString output = "AnalysisResults.root";
  AliAnalysisDataContainer *deuCont = mgr->CreateContainer(Form("mpuccio_%s",tskname.Data()),
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           output.Data());
  mgr->ConnectInput  (deu,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu,  1, deuCont);
  return deu;
}

