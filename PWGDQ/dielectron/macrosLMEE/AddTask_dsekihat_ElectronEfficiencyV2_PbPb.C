//// ROOT6 modifications
//#ifdef __CLING__
//#include <AliAnalysisManager.h>
//#include <AliAODInputHandler.h>
//#include <AliDielectronVarCuts.h>
//
//// Tell ROOT where to find AliPhysics headers
//R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
//#endif

void AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task);
void AddPairMCSignal(AliAnalysisTaskElectronEfficiencyV2* task);

AliAnalysisTaskElectronEfficiencyV2* AddTask_dsekihat_ElectronEfficiencyV2_PbPb(
    Bool_t getFromAlien = kFALSE,
    TString configFile="Config_dsekihat_ElectronEfficiencyV2_PbPb.C",
    UInt_t trigger = AliVEvent::kINT7,
    const Int_t CenMin =  0,
    const Int_t CenMax = 10,
    const Float_t PtMin =  0.2,
    const Float_t PtMax = 10.0,
    const Float_t EtaMin = -0.8,
    const Float_t EtaMax = +0.8,
    const TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;Pythia CC_0;Pythia BB_0;Pythia B_0;",
    const Bool_t isLHC19f2 = kTRUE,
    const std::string resolutionAlien ="",
    const std::string cocktailAlien   ="",
    const std::string centralityAlien =""
    ){

  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  Bool_t isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();
  printf("isAOD = %d\n",isAOD);

  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(Form("TaskElectronEfficiencyV2_Cen%d_%d_kINT7",CenMin,CenMax));

  TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");
  //TString configBasePath("./");
  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data())))){
    configBasePath=Form("%s/",gSystem->pwd());
  }
  TString configFilePath(configBasePath + configFile);

  //load cut library first
  TString libFilePath(configBasePath + "LMEECutLib_dsekihat.C");
  std::cout << "Configpath:  " << configFilePath << std::endl;

#if defined(__CLING__)
  printf("ROOT6\n");
  gROOT->LoadMacro(libFilePath.Data());//library first
  gROOT->LoadMacro(configFilePath.Data());
#elif defined(__CINT__)
  printf("ROOT5\n");
  gROOT->LoadMacro(libFilePath.Data());//library first
  gROOT->LoadMacro(configFilePath.Data());
#endif

  // Adding cutsettings
  Int_t nCut = gROOT->ProcessLine("GetN()");
  printf("Add %d cuts\n",nCut);
  for (int iCut = 0; iCut < nCut; ++iCut){
    AliAnalysisFilter *filter = reinterpret_cast<AliAnalysisFilter*>(gROOT->ProcessLine(Form("Config_dsekihat_ElectronEfficiencyV2_PbPb(%d,%d,%f,%f,%f,%f)",iCut,isAOD,PtMin,PtMax,EtaMin,EtaMax)));
    task->AddTrackCuts(filter);
  }

  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);//always ON in Run2 analyses for both data and MC.
  task->SetTriggerMask(trigger);
  task->SetEventFilter(reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("LMEECutLib::SetupEventCuts(%f,%f,%d,\"%s\")",(Float_t)CenMin,(Float_t)CenMax,kTRUE,"V0M"))));//kTRUE is for Run2
  //task->SetCentralityEstimator("V0M");
  //task->SetCentrality(CenMin, CenMax);

  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(0.1);
  task->SetMaxPtGen(1e+10);
  task->SetMinEtaGen(-1.5);
  task->SetMaxEtaGen(+1.5);

  // Set minimum and maximum values for pairing
  task->SetKinematicCuts(PtMin, PtMax, EtaMin, EtaMax);

  // Set Binning
  task->SetPtBinsLinear    (0, 10, 1000);
  task->SetEtaBinsLinear   (-1, +1, 20);
  task->SetPhiBinsLinear   (0, TMath::TwoPi(), 90);
  task->SetThetaBinsLinear (0, TMath::TwoPi(), 60);
  task->SetMassBinsLinear  (0, 5, 500);
  task->SetPairPtBinsLinear(0, 20, 200);

  task->SetSmearGenerated(kFALSE);
  task->SetResolutionDeltaPtBinsLinear( -10,  +10, 2000);
  task->SetResolutionRelPtBinsLinear  (   0,    2,  400);
  task->SetResolutionEtaBinsLinear    (-0.4, +0.4,  200);
  task->SetResolutionPhiBinsLinear    (-0.4, +0.4,  200);
  task->SetResolutionThetaBinsLinear  (-0.4, +0.4,  200);

  // Set MCSignal and Cutsetting to fill the support histograms
  task->SetSupportHistoMCSignalAndCutsetting(0,0);//fill support histograms for first MCsignal and first cutsetting

  // Pairing related config
  task->SetDoPairing(kTRUE);
  //task->SetULSandLS(kTRUE);
  //task->SetDeactivateLS(kTRUE);

  //TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;";
  //TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;Pythia CC_0;Pythia BB_0;Pythia B_0;";
  //TString generators = "Pythia CC_0;Pythia BB_0;Pythia B_0;";
  //TString generators = "Pythia CC_0;";
  //TString generators = "Pythia BB_0;Pythia B_0;";
  //TString generators = "Pythia BB_0;";
  //TString generators = "Pythia B_0";

  //TString generators = reinterpret_cast<TString>(gROOT->ProcessLine("LMEECutLib::GetGeneratorMCSignalName()"));
  //TString generators = TString(gROOT->ProcessLine("LMEECutLib::GetGeneratorMCSignalName()"));

  cout<<"Efficiency based on MC generators: " << generators <<endl;
  TString generatorsPair=generators;
  task->SetGeneratorMCSignalName(generatorsPair);
  task->SetGeneratorULSSignalName(generators);
  task->SetLHC19f2MC(isLHC19f2);

  // Resolution File, If resoFilename = "" no correction is applied
  task->SetCentralityFile(centralityAlien);
  task->SetResolutionFileFromAlien(resolutionAlien);
  task->SetCocktailWeightingFromAlien(cocktailAlien);

  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  AddSingleLegMCSignal(task);
  AddPairMCSignal(task);

  const TString fileName = AliAnalysisManager::GetCommonFileName();
	const TString dirname = Form("PWGDQ_LMEE_ElectronEfficiencyV2_Cen%d_%d_kINT7",CenMin,CenMax);
  mgr->AddTask(task);
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("Efficiency_Cen%d_%d_kINT7",CenMin,CenMax), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",fileName.Data(),dirname.Data())));
  return task;
}
//___________________________________________________________________
void AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){
  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalState.SetMotherPDGs(22, 22, kTRUE, kTRUE); // Exclude conversion electrons
  task->AddSingleLegMCSignal(eleFinalState);

  //AliDielectronSignalMC eleFinalStateFromLF("eleFinalStateFromLF","eleFinalStateFromLF");
  //eleFinalStateFromLF.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  //eleFinalStateFromLF.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //eleFinalStateFromLF.SetMotherPDGs(600, 600);
  //eleFinalStateFromLF.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //task->AddSingleLegMCSignal(eleFinalStateFromLF);

//  AliDielectronSignalMC eleFinalStateFromPion("eleFinalStateFromPion","eleFinalStateFromPion");
//  eleFinalStateFromPion.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
//  eleFinalStateFromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalStateFromPion.SetMotherPDGs(111, 111);
//  eleFinalStateFromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  task->AddSingleLegMCSignal(eleFinalStateFromPion);
//
//  AliDielectronSignalMC eleFinalStateFromEta("eleFinalStateFromEta","eleFinalStateFromEta");
//  eleFinalStateFromEta.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
//  eleFinalStateFromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalStateFromEta.SetMotherPDGs(221, 221);
//  eleFinalStateFromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  task->AddSingleLegMCSignal(eleFinalStateFromEta);
//
//  AliDielectronSignalMC eleFinalStateFromEtaPrime("eleFinalStateFromEtaPrime","eleFinalStateFromEtaPrime");
//  eleFinalStateFromEtaPrime.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
//  eleFinalStateFromEtaPrime.SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalStateFromEtaPrime.SetMotherPDGs(331, 331);
//  eleFinalStateFromEtaPrime.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  task->AddSingleLegMCSignal(eleFinalStateFromEtaPrime);
//
//  AliDielectronSignalMC eleFinalStateFromRho("eleFinalStateFromRho","eleFinalStateFromRho");
//  eleFinalStateFromRho.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
//  eleFinalStateFromRho.SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalStateFromRho.SetMotherPDGs(113, 113);
//  eleFinalStateFromRho.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  task->AddSingleLegMCSignal(eleFinalStateFromRho);
//
//  AliDielectronSignalMC eleFinalStateFromOmega("eleFinalStateFromOmega","eleFinalStateFromOmega");
//  eleFinalStateFromOmega.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
//  eleFinalStateFromOmega.SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalStateFromOmega.SetMotherPDGs(223, 223);
//  eleFinalStateFromOmega.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  task->AddSingleLegMCSignal(eleFinalStateFromOmega);
//
//  AliDielectronSignalMC eleFinalStateFromPhi("eleFinalStateFromPhi","eleFinalStateFromPhi");
//  eleFinalStateFromPhi.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
//  eleFinalStateFromPhi.SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalStateFromPhi.SetMotherPDGs(333, 333);
//  eleFinalStateFromPhi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  task->AddSingleLegMCSignal(eleFinalStateFromPhi);
//
//  AliDielectronSignalMC eleFinalStateFromJpsi("eleFinalStateFromJpsi","eleFinalStateFromJpsi");
//  eleFinalStateFromJpsi.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
//  eleFinalStateFromJpsi.SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalStateFromJpsi.SetMotherPDGs(443, 443); // open charm mesons and baryons together
//  eleFinalStateFromJpsi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  task->AddSingleLegMCSignal(eleFinalStateFromJpsi);

  //
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddSingleLegMCSignal(eleFinalStateFromD);
  //
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502); // open beauty mesons and baryons together
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddSingleLegMCSignal(eleFinalStateFromB);

//  AliDielectronSignalMC eleFinalStateFromGamma("eleFinalStateFromGamma","eleFinalStateFromGamma");
//  eleFinalStateFromGamma.SetLegPDGs(11,-11);
//  eleFinalStateFromGamma.SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalStateFromGamma.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // for e from gamma: kFinalState is empty.
//  eleFinalStateFromGamma.SetMotherPDGs(22,22,kFALSE,kFALSE);
//  task->AddSingleLegMCSignal(eleFinalStateFromGamma);

//  // This is used to get electrons not from same mother for pair efficiency.
//  // Needed to look at D and B meson electrons as functionality to pair those is
//  // not implemented in the framework. Instead, use all final start electrons
//  // from D or B decays for efficiency correction, for example.
//  // The ordering must match the ordering of the added signals above*.
//  std::vector<Bool_t> DielectronsPairNotFromSameMother;
//  DielectronsPairNotFromSameMother.push_back(kFALSE);//all same mother
//  DielectronsPairNotFromSameMother.push_back(kTRUE);//D
//  DielectronsPairNotFromSameMother.push_back(kTRUE);//B
//  task->AddMCSignalsWhereDielectronPairNotFromSameMother(DielectronsPairNotFromSameMother);
}
//___________________________________________________________________
void AddPairMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){
  AliDielectronSignalMC pair_sameMother("sameMother_pair","sameMother_pair");
  pair_sameMother.SetLegPDGs(11,-11);
  pair_sameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_sameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //mother
  pair_sameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
  pair_sameMother.SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. should have no effect on final state ele.
  task->AddPairMCSignal(pair_sameMother);

  //AliDielectronSignalMC pair_sameMother_LF("sameMotherLF_par","sameMotherLF_pair");
  //pair_sameMother_LF.SetLegPDGs(11,-11);
  //pair_sameMother_LF.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //pair_sameMother_LF.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  ////mother
  //pair_sameMother_LF.SetMothersRelation(AliDielectronSignalMC::kSame);
  //pair_sameMother_LF.SetMotherPDGs(600,600); //
  //task->AddPairMCSignal(pair_sameMother_LF);


  //    AliDielectronSignalMC pair_sameMother_pion("sameMother_pion","sameMother_pion");
  //    pair_sameMother_pion.SetLegPDGs(11,-11);
  //    pair_sameMother_pion.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //    pair_sameMother_pion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //    //mother
  //    pair_sameMother_pion.SetMothersRelation(AliDielectronSignalMC::kSame);
  //    pair_sameMother_pion.SetMotherPDGs(111,111); //
  //    task->AddPairMCSignal(pair_sameMother_pion);
  //
  //    AliDielectronSignalMC pair_sameMother_eta("sameMother_eta","sameMother_eta");
  //    pair_sameMother_eta.SetLegPDGs(11,-11);
  //    pair_sameMother_eta.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //    pair_sameMother_eta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //    //mother
  //    pair_sameMother_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
  //    pair_sameMother_eta.SetMotherPDGs(221,221); //
  //    task->AddPairMCSignal(pair_sameMother_eta);
  //
  //    AliDielectronSignalMC pair_sameMother_etaprime("sameMother_etaprime","sameMother_etaprime");
  //    pair_sameMother_etaprime.SetLegPDGs(11,-11);
  //    pair_sameMother_etaprime.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //    pair_sameMother_etaprime.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //    //mother
  //    pair_sameMother_etaprime.SetMothersRelation(AliDielectronSignalMC::kSame);
  //    pair_sameMother_etaprime.SetMotherPDGs(331,331); //
  //    task->AddPairMCSignal(pair_sameMother_etaprime);
  //
  //    AliDielectronSignalMC pair_sameMother_rho("sameMother_rho","sameMother_rho");
  //    pair_sameMother_rho.SetLegPDGs(11,-11);
  //    pair_sameMother_rho.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //    pair_sameMother_rho.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //    //mother
  //    pair_sameMother_rho.SetMothersRelation(AliDielectronSignalMC::kSame);
  //    pair_sameMother_rho.SetMotherPDGs(113,113); //
  //    task->AddPairMCSignal(pair_sameMother_rho);
  //
  //    AliDielectronSignalMC pair_sameMother_omega("sameMother_omega","sameMother_omega");
  //    pair_sameMother_omega.SetLegPDGs(11,-11);
  //    pair_sameMother_omega.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //    pair_sameMother_omega.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //    //mother
  //    pair_sameMother_omega.SetMothersRelation(AliDielectronSignalMC::kSame);
  //    pair_sameMother_omega.SetMotherPDGs(223,223); //
  //    task->AddPairMCSignal(pair_sameMother_omega);
  //
  //    AliDielectronSignalMC pair_sameMother_phi("sameMother_phi","sameMother_phi");
  //    pair_sameMother_phi.SetLegPDGs(11,-11);
  //    pair_sameMother_phi.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //    pair_sameMother_phi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //    //mother
  //    pair_sameMother_phi.SetMothersRelation(AliDielectronSignalMC::kSame);
  //    pair_sameMother_phi.SetMotherPDGs(333,333); //
  //    task->AddPairMCSignal(pair_sameMother_phi);
  //
  //    AliDielectronSignalMC pair_sameMother_jpsi("sameMother_jpsi","sameMother_jpsi");
  //    pair_sameMother_jpsi.SetLegPDGs(11,-11);
  //    pair_sameMother_jpsi.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //    pair_sameMother_jpsi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //    //mother
  //    pair_sameMother_jpsi.SetMothersRelation(AliDielectronSignalMC::kSame);
  //    pair_sameMother_jpsi.SetMotherPDGs(443,443); //
  //    task->AddPairMCSignal(pair_sameMother_jpsi);

  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD_pair","eleFinalStateFromD_pair");
  eleFinalStateFromD.SetLegPDGs(11,-11);
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //mother
  eleFinalStateFromD.SetMotherPDGs(402, 402); //
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddPairMCSignal(eleFinalStateFromD);

  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB_pair","eleFinalStateFromB_pair");
  eleFinalStateFromB.SetLegPDGs(11,-11);
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //mother
  eleFinalStateFromB.SetMotherPDGs(502, 502); //
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddPairMCSignal(eleFinalStateFromB);

//  AliDielectronSignalMC eleFinalStateFromGamma("eleFinalStateFromGamma_pair","eleFinalStateFromGamma_pair");
//  eleFinalStateFromGamma.SetLegPDGs(11,-11);
//  eleFinalStateFromGamma.SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalStateFromGamma.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
//  //mother
//  eleFinalStateFromGamma.SetMothersRelation(AliDielectronSignalMC::kSame);
//  eleFinalStateFromGamma.SetMotherPDGs(22, 22); //
//  task->AddPairMCSignal(eleFinalStateFromGamma);

}
//___________________________________________________________________

