#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>

TString TrackCutnames[] = {
// "ResolutionTrackCuts"
 //,"DefaultTrackCut_Nsc0"
//  "DefaultTrackCut_Nsc01"
 "DefaultTrackCut_Nsc01_woPU"
// ,"DefaultTrackCut_Nsc01_onlyPU"
//,"LooseTrackCut"
//,"TightTrackCut"
//,"PIDCalibTrackCut"
};
const Int_t nTC = sizeof(TrackCutnames)/sizeof(TrackCutnames[0]);
Int_t GetNTC(){return nTC;}

TString PIDnames[] = {
 "noPID"
// ,"DefaultPID"
// ,"ITSTPChadrejORTOFrec"
// ,"ITSTPChadrej"
// ,"ITSTOFrecover"
// "TPChadrejORTOFrec"
// ,"TPChadrej"
// ,"TOFrecover"
};
const Int_t nPID = sizeof(PIDnames)/sizeof(PIDnames[0]);
Int_t GetNPID(){return nPID;}

TString PFnames[] = {
  "woPF"
// ,"wPF"
};
const Int_t nPF = sizeof(PFnames)/sizeof(PFnames[0]);
Int_t GetNPF(){return nPF;}

const Int_t nDie = nTC * nPID * nPF;
Int_t GetN(){return nDie;}

void AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task);
void AddPairMCSignal(AliAnalysisTaskElectronEfficiencyV2* task);

AliAnalysisFilter *Config_dsekihat_ElectronEfficiencyV2_PbPb(
		const Int_t cutDefinitionTC,
		const Int_t cutDefinitionPID,
		const Int_t cutDefinitionPF,
    const Bool_t applyPairCut
    )
{

  TString name = Form("%s_%s_%s_PairCut%d",TrackCutnames[cutDefinitionTC].Data(),PIDnames[cutDefinitionPID].Data(),PFnames[cutDefinitionPF].Data(),applyPairCut);

  Bool_t isAOD = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();
  printf("isAOD = %d\n",isAOD);

  AliAnalysisFilter *anaFilter = new AliAnalysisFilter(Form("anaFilter_%s",name.Data()),Form("anaFilter_%d_%d_%d",cutDefinitionTC,cutDefinitionPID,cutDefinitionPF)); // named constructor seems mandatory!
  LMEECutLib *lib = new LMEECutLib(name); 
  anaFilter->AddCuts(lib->SetupTrackCuts());

  if(name.Contains("Resolution",TString::kIgnoreCase) || name.Contains("noPID",TString::kIgnoreCase)){
    printf("Do not add PID cut for ResolutionTrackCuts/noPID\n");
  }
  else anaFilter->AddCuts(lib->SetupPIDCuts());

  if(!isAOD)       anaFilter->AddCuts(lib->SetupESDtrackCuts());
	if(applyPairCut) anaFilter->AddCuts(lib->SetupPairCuts());

  anaFilter->Print();
  return anaFilter;

}
//___________________________________________________________________
void AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){
  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalState.SetMotherPDGs(22, 22, kTRUE, kTRUE); // Exclude conversion electrons
  task->AddSingleLegMCSignal(eleFinalState);

  AliDielectronSignalMC eleFinalStateFromLF("eleFinalStateFromLF","eleFinalStateFromLF");
  eleFinalStateFromLF.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromLF.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromLF.SetMotherPDGs(600, 600);
  eleFinalStateFromLF.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  task->AddSingleLegMCSignal(eleFinalStateFromLF);

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

  //AliDielectronSignalMC eleFinalStateFromGamma("eleFinalStateFromGamma","eleFinalStateFromGamma");
  //eleFinalStateFromGamma.SetLegPDGs(11,-11);
  //eleFinalStateFromGamma.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //eleFinalStateFromGamma.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // for e from gamma: kFinalState is empty.
  //eleFinalStateFromGamma.SetMotherPDGs(22,22,kFALSE,kFALSE);
  //task->AddSingleLegMCSignal(eleFinalStateFromGamma);

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

  AliDielectronSignalMC pair_sameMother_LF("sameMotherLF_par","sameMotherLF_pair");
  pair_sameMother_LF.SetLegPDGs(11,-11);
  pair_sameMother_LF.SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_sameMother_LF.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //mother
  pair_sameMother_LF.SetMothersRelation(AliDielectronSignalMC::kSame);
  pair_sameMother_LF.SetMotherPDGs(600,600);
  task->AddPairMCSignal(pair_sameMother_LF);


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

  AliDielectronSignalMC diEleOpenCharm("DMesons","di-electrons from open charm hadrons no B grandmother");
  diEleOpenCharm.SetLegPDGs(11,-11);
  diEleOpenCharm.SetMotherPDGs(402,402);
  diEleOpenCharm.SetGrandMotherPDGs(0,0);//daiki
  diEleOpenCharm.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenCharm.SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenCharm.SetCheckStackForPDG(kTRUE);
  diEleOpenCharm.SetPDGforStack(503);
  diEleOpenCharm.SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharm.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddPairMCSignal(diEleOpenCharm);


//B meson (3)
  AliDielectronSignalMC diEleOneOpenB("B2ee","di-electrons from one B meson");
  diEleOneOpenB.SetLegPDGs(11,-11);
  diEleOneOpenB.SetMotherPDGs(402,502);
  diEleOneOpenB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOneOpenB.SetGrandMotherPDGs(501,0);
  diEleOneOpenB.SetCheckMotherGrandmotherRelation(kTRUE,kTRUE);
  diEleOneOpenB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOneOpenB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOneOpenB.SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  task->AddPairMCSignal(diEleOneOpenB);

  //B meson (1)(1)
  AliDielectronSignalMC diEleOpenB("BMesons","di-electrons from B mesons");
  diEleOpenB.SetLegPDGs(11,-11);
  diEleOpenB.SetMotherPDGs(502,502);
  diEleOpenB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenB.SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddPairMCSignal(diEleOpenB);


  //B meson (2)(2)
  AliDielectronSignalMC diEleOpenBtoD("B2D2ee","di-electrons from B->D-> e");
  diEleOpenBtoD.SetLegPDGs(11,-11);
  diEleOpenBtoD.SetMotherPDGs(402,402);
  diEleOpenBtoD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenBtoD.SetGrandMotherPDGs(502,502);
  diEleOpenBtoD.SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenBtoD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBtoD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBtoD.SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  task->AddPairMCSignal(diEleOpenBtoD);


  //B meson (1)(2)
  AliDielectronSignalMC diEleOpenBandBtoD("B2eAndB2D2e","di-electrons from B->e and B->D->e");
  diEleOpenBandBtoD.SetLegPDGs(11,-11);
  diEleOpenBandBtoD.SetMotherPDGs(402,502);
  diEleOpenBandBtoD.SetGrandMotherPDGs(502,0);
  diEleOpenBandBtoD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenBandBtoD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBandBtoD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBandBtoD.SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  //do i need this?
  diEleOpenBandBtoD.SetCheckMotherGrandmotherRelation(kTRUE,kFALSE);
  task->AddPairMCSignal(diEleOpenBandBtoD);

  //AliDielectronSignalMC eleFinalStateFromGamma("eleFinalStateFromGamma_pair","eleFinalStateFromGamma_pair");
  //eleFinalStateFromGamma.SetLegPDGs(11,-11);
  //eleFinalStateFromGamma.SetCheckBothChargesLegs(kTRUE,kTRUE);
  //eleFinalStateFromGamma.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  //mother
  //eleFinalStateFromGamma.SetMothersRelation(AliDielectronSignalMC::kSame);
  //eleFinalStateFromGamma.SetMotherPDGs(22, 22);
  //task->AddPairMCSignal(eleFinalStateFromGamma);

}
//___________________________________________________________________
