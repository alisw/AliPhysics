#include <iostream>

AliAnalysisTaskHadronElectronAnalysis* AddHadronElectronAnalysisTask(
  TString INVESTIGATION_NAME,
  TString SUBMISSION_TIME,
  TString name = "hadronElectronAnalysis",
  float ETA_CUT=0.8,
  double INVARIANT_MASS_CUT=0.14,
  double PARTNER_ELECTRON_PT_CUT=0.3,
  double TPC_NSIGMA_ELEC_MIN=-0.7,
  double TPC_NSIGMA_ELEC_MAX=3,
  double PARTNER_TPC_NSIGMA_ELEC_MIN=-3,
  double PARTNER_TPC_NSIGMA_ELEC_MAX=3,
  float ASSOC_TRK_BIT_ELEC=AliAODTrack::kTrkGlobalNoDCA,
  float ASSOC_TRK_BIT_H=1024,
  float PARTNER_TRK_BIT = AliAODTrack::kTrkTPCOnly,
  float TRIG_TRK_BIT = AliAODTrack::kIsHybridGCG,
  bool USE_SPD_KBOTH = false,
  bool RUN_ON_MC = false
) {

  AliAnalysisManager *manage = AliAnalysisManager::GetAnalysisManager();

  if (!manage) return 0x0;

  if(!manage->GetInputEventHandler()) return 0x0;



  TString file_name = Form("%s_partner_electron_pt_cut_%f_tpc_nsigma_elec_min_%f_tpc_nsigma_elec_max_%f_assoc_trk_bit_elec_%f_time_%s.root", INVESTIGATION_NAME.Data(), PARTNER_ELECTRON_PT_CUT, TPC_NSIGMA_ELEC_MIN, TPC_NSIGMA_ELEC_MAX, ASSOC_TRK_BIT_ELEC, SUBMISSION_TIME.Data());
  
  AliAnalysisTaskHadronElectronAnalysis* task = new AliAnalysisTaskHadronElectronAnalysis(name.Data());

  if(!task) return 0x0;

  task->LoadEfficiencies();

  task->SetEtaCut(ETA_CUT);
  task->SetInvariantMassCut(INVARIANT_MASS_CUT);
  task->SetPartnerElectronPtCut(PARTNER_ELECTRON_PT_CUT);
  task->SetTPCNSigmaElecMin(TPC_NSIGMA_ELEC_MIN);
  task->SetTPCNSigmaElecMax(TPC_NSIGMA_ELEC_MAX);
  task->SetPartnerTPCNSigmaElecMin(PARTNER_TPC_NSIGMA_ELEC_MIN);
  task->SetPartnerTPCNSigmaElecMax(PARTNER_TPC_NSIGMA_ELEC_MAX);
  task->SetAssocTrkBitElec(ASSOC_TRK_BIT_ELEC);
  task->SetAssocTrkBitH(ASSOC_TRK_BIT_H);
  task->SetPartnerTrkBit(PARTNER_TRK_BIT);
  task->SetTrigTrkBit(TRIG_TRK_BIT);
  task->SetUseSPDKBoth(USE_SPD_KBOTH);
  task->SetRunOnMC(RUN_ON_MC);

  manage->AddTask(task);

  manage->ConnectInput(task, 0, manage->GetCommonInputContainer());
  manage->ConnectOutput(task, 1, manage->CreateContainer("hadron-electron", TList::Class(), AliAnalysisManager::kOutputContainer, file_name.Data()));

  return task;

}
