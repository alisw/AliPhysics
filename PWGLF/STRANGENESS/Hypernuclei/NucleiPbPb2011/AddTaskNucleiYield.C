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

void AddTaskNucleiYield(Bool_t isMC = kFALSE, TString output = "nuclei"){
  
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNucleiYield", "No analysis manager found.");
    return;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskNucleiYield", "This task requires an input event handler");
    return;
  }
  
  // Common variables
  if (isMC) output.Append("MC.root");
  else output.Append(".root");
  Float_t deuteronCentBins[5] = {0.,10.,20.,40.,60.};
  Float_t tritonCentBins[4] = {0.,10.,30.,50.};
  Float_t deuteronPtBins[27] = {
    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
    3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f
  };
  Float_t deutBBpar[5] = {4.69637f,7.51827f,0.0183746f,2.60f,2.7f};
  Float_t sigmaBBpar = 0.1f;
  Float_t tritonPtBins[17] = {
    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f
  };
  TString tskname;
  
  /// ### Deuterons
  tskname = "deu4cent";
  AliAnalysisTaskNucleiYield *deu = new AliAnalysisTaskNucleiYield(tskname);
  deu->SetParticleType(AliPID::kDeuteron);
  deu->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu->SetPDG(1000010020);
  deu->SetIsMC(isMC);
  deu->SetCentBins(4, deuteronCentBins);
  deu->SetPtBins(26,deuteronPtBins);
  deu->SetDCABins(80,-0.5,0.5);
  mgr->AddTask(deu);
  AliAnalysisDataContainer *deuCont = mgr->CreateContainer(Form("mp_%s",tskname.Data()),
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           output.Data());
  mgr->ConnectInput  (deu,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu,  1, deuCont);
  
  /// ### Deuterons with other centralities classes
  tskname = "deu3cent";
  AliAnalysisTaskNucleiYield *deu2 = new AliAnalysisTaskNucleiYield(tskname);
  deu2->SetParticleType(AliPID::kDeuteron);
  deu2->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu2->SetPDG(1000010020);
  deu2->SetIsMC(isMC);
  deu2->SetCentBins(3, tritonCentBins);
  deu2->SetPtBins(26,deuteronPtBins);
  deu2->SetDCABins(80,-0.5,0.5);
  mgr->AddTask(deu2);
  AliAnalysisDataContainer *deuCont2 = mgr->CreateContainer(Form("mp_%s",tskname.Data()),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            output.Data());
  mgr->ConnectInput  (deu2,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu2,  1, deuCont2);
  
  /// \f$\chi^{2}\f$ studies
  tskname = "deu3cent_chi2low";
  AliAnalysisTaskNucleiYield *deu3 = new AliAnalysisTaskNucleiYield(tskname);
  deu3->SetParticleType(AliPID::kDeuteron);
  deu3->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu3->SetPDG(1000010020);
  deu3->SetIsMC(isMC);
  deu3->SetCentBins(3, tritonCentBins);
  deu3->SetPtBins(26,deuteronPtBins);
  deu3->SetDCABins(80,-0.5,0.5);
  deu3->SetRequireMaxChi2(3.5);
  mgr->AddTask(deu3);
  AliAnalysisDataContainer *deuCont3 = mgr->CreateContainer(Form("mp_%s",tskname.Data()),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            output.Data());
  mgr->ConnectInput  (deu3,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu3,  1, deuCont3);

  tskname = "deu3cent_chi2high";
  AliAnalysisTaskNucleiYield *deu4 = new AliAnalysisTaskNucleiYield(tskname);
  deu4->SetParticleType(AliPID::kDeuteron);
  deu4->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu4->SetPDG(1000010020);
  deu4->SetIsMC(isMC);
  deu4->SetCentBins(3, tritonCentBins);
  deu4->SetPtBins(26,deuteronPtBins);
  deu4->SetDCABins(80,-0.5,0.5);
  deu4->SetRequireMaxChi2(6.);
  mgr->AddTask(deu4);
  AliAnalysisDataContainer *deuCont4 = mgr->CreateContainer(Form("mp_%s",tskname.Data()),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            output.Data());
  mgr->ConnectInput  (deu4,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu4,  1, deuCont4);

  /// ### Studies using different TPC minimum number of clusters
  tskname = "deu3cent_tpclow";
  AliAnalysisTaskNucleiYield *deu5 = new AliAnalysisTaskNucleiYield(tskname);
  deu5->SetParticleType(AliPID::kDeuteron);
  deu5->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu5->SetPDG(1000010020);
  deu5->SetIsMC(isMC);
  deu5->SetCentBins(3, tritonCentBins);
  deu5->SetPtBins(26,deuteronPtBins);
  deu5->SetDCABins(80,-0.5,0.5);
  deu5->SetRequireTPCsignal(60);
  mgr->AddTask(deu5);
  AliAnalysisDataContainer *deuCont5 = mgr->CreateContainer(Form("mp_%s",tskname.Data()),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            output.Data());
  mgr->ConnectInput  (deu5,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu5,  1, deuCont5);
  
  tskname = "deu3cent_tpchigh";
  AliAnalysisTaskNucleiYield *deu6 = new AliAnalysisTaskNucleiYield(tskname);
  deu6->SetParticleType(AliPID::kDeuteron);
  deu6->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu6->SetPDG(1000010020);
  deu6->SetIsMC(isMC);
  deu6->SetCentBins(3, tritonCentBins);
  deu6->SetPtBins(26,deuteronPtBins);
  deu6->SetDCABins(80,-0.5,0.5);
  deu6->SetRequireTPCsignal(80);
  mgr->AddTask(deu6);
  AliAnalysisDataContainer *deuCont6 = mgr->CreateContainer(Form("mp_%s",tskname.Data()),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            output.Data());
  mgr->ConnectInput  (deu6,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu6,  1, deuCont6);
  
  /// ### Studies using different cuts on DCA\f$_{z}\f$
  tskname = "deu3cent_dcazhigh";
  AliAnalysisTaskNucleiYield *deu7 = new AliAnalysisTaskNucleiYield(tskname);
  deu7->SetParticleType(AliPID::kDeuteron);
  deu7->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu7->SetPDG(1000010020);
  deu7->SetIsMC(isMC);
  deu7->SetCentBins(3, tritonCentBins);
  deu7->SetPtBins(26,deuteronPtBins);
  deu7->SetDCABins(80,-0.5,0.5);
  deu7->SetRequireMaxDCAz(2.f);
  mgr->AddTask(deu7);
  AliAnalysisDataContainer *deuCont7 = mgr->CreateContainer(Form("mp_%s",tskname.Data()),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            output.Data());
  mgr->ConnectInput  (deu7,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu7,  1, deuCont7);

  tskname = "deu3cent_dcazlow";
  AliAnalysisTaskNucleiYield *deu8 = new AliAnalysisTaskNucleiYield(tskname);
  deu8->SetParticleType(AliPID::kDeuteron);
  deu8->SetCustomTPCpid(deutBBpar, sigmaBBpar);
  deu8->SetPDG(1000010020);
  deu8->SetIsMC(isMC);
  deu8->SetCentBins(3, tritonCentBins);
  deu8->SetPtBins(26,deuteronPtBins);
  deu8->SetDCABins(80,-0.5,0.5);
  deu8->SetRequireMaxDCAz(0.5);
  mgr->AddTask(deu8);
  AliAnalysisDataContainer *deuCont8 = mgr->CreateContainer(Form("mp_%s",tskname.Data()),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            output.Data());
  mgr->ConnectInput  (deu8,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu8,  1, deuCont8);
  
  /// ### Tritons
  tskname = "tritons";
  AliAnalysisTaskNucleiYield *tri = new AliAnalysisTaskNucleiYield(tskname);
  tri->SetParticleType(AliPID::kTriton);
  tri->SetPDG(1000010030);
  tri->SetIsMC(isMC);
  tri->SetCentBins(3, tritonCentBins);
  tri->SetPtBins(16,tritonPtBins);
  tri->SetDCABins(80,-0.5,0.5);
  mgr->AddTask(tri);
  AliAnalysisDataContainer *triCont = mgr->CreateContainer("mp_tri", TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           output.Data());
  mgr->ConnectInput  (tri,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (tri,  1, triCont);
  
  /// ### \f$^{3}\f$He
  tskname = "he3";
  AliAnalysisTaskNucleiYield *he3 = new AliAnalysisTaskNucleiYield(tskname);
  he3->SetParticleType(AliPID::kHe3);
  he3->SetPDG(1000020030);
  he3->SetIsMC(isMC);
  he3->SetCentBins(3, tritonCentBins);
  he3->SetPtBins(16,tritonPtBins);
  he3->SetDCABins(80,-0.5,0.5);
  mgr->AddTask(he3);
  AliAnalysisDataContainer *he3Cont = mgr->CreateContainer("mp_he3", TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           output.Data());
  mgr->ConnectInput  (he3,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (he3,  1, he3Cont);
}