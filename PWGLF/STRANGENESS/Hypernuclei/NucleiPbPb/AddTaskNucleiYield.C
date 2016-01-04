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
  Float_t tritonCentBins[4] = {0.,10.,30.,50.};
  Float_t deuteronPtBins[27] = {
    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
    3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f
  };
  Float_t deutBBpar[5] = {4.69637f,7.51827f,0.0183746f,2.60f,2.7f};
  Float_t sigmaBBpar = 0.1f;

  TString tskname = "deuterons";
  tskname.Append(Form("%s",suffix.Data()));
  AliAnalysisTaskNucleiYield *deu = new AliAnalysisTaskNucleiYield(tskname);
  deu->SetParticleType(AliPID::kDeuteron);
  deu->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu->SetPDG(pdgCode);
  deu->SetIsMC(isMC);
  deu->SetCentBins(3, tritonCentBins);
  deu->SetPtBins(26,deuteronPtBins);
  deu->SetDCABins(80,-0.5,0.5);
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
