PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt* AddTaskJetHardestKt(
 const char* njetsBase, const char* njetsUS, const char* njetsTrue, const char* njetsPartLevel, const Double_t R,
 const char* nrhoBase, const char* ntracks, const char* ntracksUS, const char* ntracksPartLevel, const char* nclusters,
 const char* ntracksTrue, const char* type, const char* CentEst, Int_t pSel,
 PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt::JetShapeType_t jetShapeType =
  PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt::kMCTrue,
 PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt::JetShapeSub_t jetShapeSub =
  PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt::kNoSub,
 PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt::JetSelectionType_t jetSelection =
  PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt::kInclusive,
 Float_t minpTHTrigger = 0., Float_t maxpTHTrigger = 0., Float_t acut = 0.6,
 PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt::DerivSubtrOrder_t derivSubtrOrder =
  PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt::kSecondOrder,
 const std::string& suffix = "")
{
  PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt* task =
   PWGJE::EMCALJetTasks::AliAnalysisTaskJetHardestKt::AddTaskJetHardestKt(
    njetsBase, njetsUS, njetsTrue, njetsPartLevel, R, nrhoBase, ntracks, ntracksUS, ntracksPartLevel, nclusters,
    ntracksTrue, type, CentEst, pSel, jetShapeType, jetShapeSub, jetSelection, minpTHTrigger, maxpTHTrigger, acut,
    derivSubtrOrder, suffix);
  return task;
}
