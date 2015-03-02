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
//  Float_t deuteronCentBins[5] = {0.,10.,20.,40.,60.};
  Float_t tritonCentBins[4] = {0.,10.,30.,50.};
  Float_t deuteronPtBins[27] = {
    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
    3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f
  };
  Float_t deutBBpar[5] = {4.69637f,7.51827f,0.0183746f,2.60f,2.7f};
  Float_t sigmaBBpar = 0.1f;
//  Float_t tritonPtBins[17] = {
//    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
//    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f
//  };
//  
  TString tskname = "deuterons";
//  if (isMC) tskname.Append("_MC");
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
  AliAnalysisDataContainer *deuCont = mgr->CreateContainer(Form("mpuccio:%s",tskname.Data()),
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           output.Data());
  mgr->ConnectInput  (deu,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (deu,  1, deuCont);
  return deu;
  
//  tskname = "deuteron4cent_chi2low";
//  tskname.Append(Form("_%s",suffix.Data()));
//  AliAnalysisTaskNucleiYield *deu41 = new AliAnalysisTaskNucleiYield(tskname);
//  deu41->SetParticleType(AliPID::kDeuteron);
//  deu41->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu41->SetPDG(pdgCode);
//  deu41->SetIsMC(isMC);
//  deu41->SetCentBins(4, deuteronCentBins);
//  deu41->SetPtBins(26,deuteronPtBins);
//  deu41->SetDCABins(80,-0.5,0.5);
//  deu41->SetRequireMaxChi2(3.5);
//  mgr->AddTask(deu41);
//  AliAnalysisDataContainer *deuCont41 = mgr->CreateContainer(tskname.Data(),
//                                                           TList::Class(),
//                                                           AliAnalysisManager::kOutputContainer,
//                                                           output.Data());
//  mgr->ConnectInput  (deu41,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu41,  1, deuCont41);
//  
//  tskname = "deuteron4cent_chi2high";
//  tskname.Append(Form("_%s",suffix.Data()));
//  AliAnalysisTaskNucleiYield *deu42 = new AliAnalysisTaskNucleiYield(tskname);
//  deu42->SetParticleType(AliPID::kDeuteron);
//  deu42->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu42->SetPDG(pdgCode);
//  deu42->SetIsMC(isMC);
//  deu42->SetCentBins(4, deuteronCentBins);
//  deu42->SetPtBins(26,deuteronPtBins);
//  deu42->SetDCABins(80,-0.5,0.5);
//  deu42->SetRequireMaxChi2(6.);
//  mgr->AddTask(deu42);
//  AliAnalysisDataContainer *deuCont42 = mgr->CreateContainer(tskname.Data(),
//                                                             TList::Class(),
//                                                             AliAnalysisManager::kOutputContainer,
//                                                             output.Data());
//  mgr->ConnectInput  (deu42,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu42,  1, deuCont42);
//  
//  tskname = "deuteron4cent_tpclow";
//  tskname.Append(Form("_%s",suffix.Data()));
//  AliAnalysisTaskNucleiYield *deu43 = new AliAnalysisTaskNucleiYield(tskname);
//  deu43->SetParticleType(AliPID::kDeuteron);
//  deu43->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu43->SetPDG(pdgCode);
//  deu43->SetIsMC(isMC);
//  deu43->SetCentBins(4, deuteronCentBins);
//  deu43->SetPtBins(26,deuteronPtBins);
//  deu43->SetDCABins(80,-0.5,0.5);
//  deu43->SetRequireTPCsignal(60);
//  mgr->AddTask(deu43);
//  AliAnalysisDataContainer *deuCont43 = mgr->CreateContainer(tskname.Data(),
//                                                             TList::Class(),
//                                                             AliAnalysisManager::kOutputContainer,
//                                                             output.Data());
//  mgr->ConnectInput  (deu43,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu43,  1, deuCont43);
//  
//  tskname = "deuteron4cent_tpchigh";
//  tskname.Append(Form("_%s",suffix.Data()));
//  AliAnalysisTaskNucleiYield *deu44 = new AliAnalysisTaskNucleiYield(tskname);
//  deu44->SetParticleType(AliPID::kDeuteron);
//  deu44->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu44->SetPDG(pdgCode);
//  deu44->SetIsMC(isMC);
//  deu44->SetCentBins(4, deuteronCentBins);
//  deu44->SetPtBins(26,deuteronPtBins);
//  deu44->SetDCABins(80,-0.5,0.5);
//  deu44->SetRequireTPCsignal(80);
//  mgr->AddTask(deu44);
//  AliAnalysisDataContainer *deuCont44 = mgr->CreateContainer(tskname.Data(),
//                                                             TList::Class(),
//                                                             AliAnalysisManager::kOutputContainer,
//                                                             output.Data());
//  mgr->ConnectInput  (deu44,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu44,  1, deuCont44);
//  
//  tskname = "deuteron4cent_dcazlow";
//  tskname.Append(Form("_%s",suffix.Data()));
//  AliAnalysisTaskNucleiYield *deu45 = new AliAnalysisTaskNucleiYield(tskname);
//  deu45->SetParticleType(AliPID::kDeuteron);
//  deu45->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu45->SetPDG(pdgCode);
//  deu45->SetIsMC(isMC);
//  deu45->SetCentBins(4, deuteronCentBins);
//  deu45->SetPtBins(26,deuteronPtBins);
//  deu45->SetDCABins(80,-0.5,0.5);
//  deu45->SetRequireMaxDCAz(0.5);
//  mgr->AddTask(deu45);
//  AliAnalysisDataContainer *deuCont45 = mgr->CreateContainer(tskname.Data(),
//                                                             TList::Class(),
//                                                             AliAnalysisManager::kOutputContainer,
//                                                             output.Data());
//  mgr->ConnectInput  (deu45,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu45,  1, deuCont45);
//
//  tskname = "deuteron4cent_dcazhigh";
//  tskname.Append(Form("_%s",suffix.Data()));
//  AliAnalysisTaskNucleiYield *deu46 = new AliAnalysisTaskNucleiYield(tskname);
//  deu46->SetParticleType(AliPID::kDeuteron);
//  deu46->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu46->SetPDG(pdgCode);
//  deu46->SetIsMC(isMC);
//  deu46->SetCentBins(4, deuteronCentBins);
//  deu46->SetPtBins(26,deuteronPtBins);
//  deu46->SetDCABins(80,-0.5,0.5);
//  deu46->SetRequireMaxDCAz(2.);
//  mgr->AddTask(deu46);
//  AliAnalysisDataContainer *deuCont46 = mgr->CreateContainer(tskname.Data(),
//                                                             TList::Class(),
//                                                             AliAnalysisManager::kOutputContainer,
//                                                             output.Data());
//  mgr->ConnectInput  (deu46,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu46,  1, deuCont46);
//
////
////  return deu;
//  
//  /// ### Deuterons with other centralities classes
//  tskname = "deuteron3cent";
//  AliAnalysisTaskNucleiYield *deu2 = new AliAnalysisTaskNucleiYield(tskname);
//  deu2->SetParticleType(AliPID::kDeuteron);
//  deu2->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu2->SetPDG(1000010020);
//  deu2->SetIsMC(isMC);
//  deu2->SetCentBins(3, tritonCentBins);
//  deu2->SetPtBins(26,deuteronPtBins);
//  deu2->SetDCABins(80,-0.5,0.5);
//  mgr->AddTask(deu2);
//  AliAnalysisDataContainer *deuCont2 = mgr->CreateContainer(tskname.Data(),
//                                                            TList::Class(),
//                                                            AliAnalysisManager::kOutputContainer,
//                                                            output.Data());
//  mgr->ConnectInput  (deu2,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu2,  1, deuCont2);
//  
//  /// \f$\chi^{2}\f$ studies
//  tskname = "deuteron3cent_chi2low";
//  AliAnalysisTaskNucleiYield *deu3 = new AliAnalysisTaskNucleiYield(tskname);
//  deu3->SetParticleType(AliPID::kDeuteron);
//  deu3->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu3->SetPDG(1000010020);
//  deu3->SetIsMC(isMC);
//  deu3->SetCentBins(3, tritonCentBins);
//  deu3->SetPtBins(26,deuteronPtBins);
//  deu3->SetDCABins(80,-0.5,0.5);
//  deu3->SetRequireMaxChi2(3.5);
//  mgr->AddTask(deu3);
//  AliAnalysisDataContainer *deuCont3 = mgr->CreateContainer(tskname.Data(),
//                                                            TList::Class(),
//                                                            AliAnalysisManager::kOutputContainer,
//                                                            output.Data());
//  mgr->ConnectInput  (deu3,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu3,  1, deuCont3);
//
//  tskname = "deuteron3cent_chi2high";
//  AliAnalysisTaskNucleiYield *deu4 = new AliAnalysisTaskNucleiYield(tskname);
//  deu4->SetParticleType(AliPID::kDeuteron);
//  deu4->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu4->SetPDG(1000010020);
//  deu4->SetIsMC(isMC);
//  deu4->SetCentBins(3, tritonCentBins);
//  deu4->SetPtBins(26,deuteronPtBins);
//  deu4->SetDCABins(80,-0.5,0.5);
//  deu4->SetRequireMaxChi2(6.);
//  mgr->AddTask(deu4);
//  AliAnalysisDataContainer *deuCont4 = mgr->CreateContainer(tskname.Data(),
//                                                            TList::Class(),
//                                                            AliAnalysisManager::kOutputContainer,
//                                                            output.Data());
//  mgr->ConnectInput  (deu4,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu4,  1, deuCont4);
//
//  /// ### Studies using different TPC minimum number of clusters
//  tskname = "deuteron3cent_tpclow";
//  AliAnalysisTaskNucleiYield *deu5 = new AliAnalysisTaskNucleiYield(tskname);
//  deu5->SetParticleType(AliPID::kDeuteron);
//  deu5->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu5->SetPDG(1000010020);
//  deu5->SetIsMC(isMC);
//  deu5->SetCentBins(3, tritonCentBins);
//  deu5->SetPtBins(26,deuteronPtBins);
//  deu5->SetDCABins(80,-0.5,0.5);
//  deu5->SetRequireTPCsignal(60);
//  mgr->AddTask(deu5);
//  AliAnalysisDataContainer *deuCont5 = mgr->CreateContainer(tskname.Data(),
//                                                            TList::Class(),
//                                                            AliAnalysisManager::kOutputContainer,
//                                                            output.Data());
//  mgr->ConnectInput  (deu5,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu5,  1, deuCont5);
//  
//  tskname = "deuteron3cent_tpchigh";
//  AliAnalysisTaskNucleiYield *deu6 = new AliAnalysisTaskNucleiYield(tskname);
//  deu6->SetParticleType(AliPID::kDeuteron);
//  deu6->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu6->SetPDG(1000010020);
//  deu6->SetIsMC(isMC);
//  deu6->SetCentBins(3, tritonCentBins);
//  deu6->SetPtBins(26,deuteronPtBins);
//  deu6->SetDCABins(80,-0.5,0.5);
//  deu6->SetRequireTPCsignal(80);
//  mgr->AddTask(deu6);
//  AliAnalysisDataContainer *deuCont6 = mgr->CreateContainer(tskname.Data(),
//                                                            TList::Class(),
//                                                            AliAnalysisManager::kOutputContainer,
//                                                            output.Data());
//  mgr->ConnectInput  (deu6,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu6,  1, deuCont6);
//  
//  /// ### Studies using different cuts on DCA\f$_{z}\f$
//  tskname = "deuteron3cent_dcazhigh";
//  AliAnalysisTaskNucleiYield *deu7 = new AliAnalysisTaskNucleiYield(tskname);
//  deu7->SetParticleType(AliPID::kDeuteron);
//  deu7->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu7->SetPDG(1000010020);
//  deu7->SetIsMC(isMC);
//  deu7->SetCentBins(3, tritonCentBins);
//  deu7->SetPtBins(26,deuteronPtBins);
//  deu7->SetDCABins(80,-0.5,0.5);
//  deu7->SetRequireMaxDCAz(2.f);
//  mgr->AddTask(deu7);
//  AliAnalysisDataContainer *deuCont7 = mgr->CreateContainer(tskname.Data(),
//                                                            TList::Class(),
//                                                            AliAnalysisManager::kOutputContainer,
//                                                            output.Data());
//  mgr->ConnectInput  (deu7,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu7,  1, deuCont7);
//
//  tskname = "deuteron3cent_dcazlow";
//  AliAnalysisTaskNucleiYield *deu8 = new AliAnalysisTaskNucleiYield(tskname);
//  deu8->SetParticleType(AliPID::kDeuteron);
//  deu8->SetCustomTPCpid(deutBBpar, sigmaBBpar);
//  deu8->SetPDG(1000010020);
//  deu8->SetIsMC(isMC);
//  deu8->SetCentBins(3, tritonCentBins);
//  deu8->SetPtBins(26,deuteronPtBins);
//  deu8->SetDCABins(80,-0.5,0.5);
//  deu8->SetRequireMaxDCAz(0.5);
//  mgr->AddTask(deu8);
//  AliAnalysisDataContainer *deuCont8 = mgr->CreateContainer(tskname.Data(),
//                                                            TList::Class(),
//                                                            AliAnalysisManager::kOutputContainer,
//                                                            output.Data());
//  mgr->ConnectInput  (deu8,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (deu8,  1, deuCont8);
  
//  /// ### Tritons
//  tskname = "tritons";
//  AliAnalysisTaskNucleiYield *tri = new AliAnalysisTaskNucleiYield(tskname);
//  tri->SetParticleType(AliPID::kTriton);
//  tri->SetPDG(1000010030);
//  tri->SetIsMC(isMC);
//  tri->SetCentBins(3, tritonCentBins);
//  tri->SetPtBins(16,tritonPtBins);
//  tri->SetDCABins(80,-0.5,0.5);
//  mgr->AddTask(tri);
//  AliAnalysisDataContainer *triCont = mgr->CreateContainer("mp_tri", TList::Class(),
//                                                           AliAnalysisManager::kOutputContainer,
//                                                           output.Data());
//  mgr->ConnectInput  (tri,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (tri,  1, triCont);
//  
//  /// ### \f$^{3}\f$He
//  tskname = "he3";
//  AliAnalysisTaskNucleiYield *he3 = new AliAnalysisTaskNucleiYield(tskname);
//  he3->SetParticleType(AliPID::kHe3);
//  he3->SetPDG(1000020030);
//  he3->SetIsMC(isMC);
//  he3->SetCentBins(3, tritonCentBins);
//  he3->SetPtBins(16,tritonPtBins);
//  he3->SetDCABins(80,-0.5,0.5);
//  mgr->AddTask(he3);
//  AliAnalysisDataContainer *he3Cont = mgr->CreateContainer("mp_he3", TList::Class(),
//                                                           AliAnalysisManager::kOutputContainer,
//                                                           output.Data());
//  mgr->ConnectInput  (he3,  0, mgr->GetCommonInputContainer());
//  mgr->ConnectOutput (he3,  1, he3Cont);

}
