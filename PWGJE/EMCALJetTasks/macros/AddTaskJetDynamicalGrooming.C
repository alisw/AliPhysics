PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming* AddTaskJetDynamicalGrooming(
 const char* njetsBase, const char* njetsUS, const char* njetsTrue, const char* njetsPartLevel, const Double_t R,
 const char* nrhoBase, const char* ntracks, const char* ntracksUS, const char* ntracksPartLevel, const char* nclusters,
 const char* ntracksTrue, const char* type, const char* CentEst, Int_t pSel,
 AliAnalysisTaskJetDynamicalGrooming::JetShapeType_t jetShapeType = AliAnalysisTaskJetDynamicalGrooming::kMCTrue,
 AliAnalysisTaskJetDynamicalGrooming::JetShapeSub_t jetShapeSub = AliAnalysisTaskJetDynamicalGrooming::kNoSub,
 AliAnalysisTaskJetDynamicalGrooming::JetSelectionType_t jetSelection =
  AliAnalysisTaskJetDynamicalGrooming::kInclusive,
 Float_t minpTHTrigger = 0., Float_t maxpTHTrigger = 0., Float_t acut = 0.6,
 AliAnalysisTaskJetDynamicalGrooming::DerivSubtrOrder_t derivSubtrOrder =
  AliAnalysisTaskJetDynamicalGrooming::kSecondOrder,
 const std::string& suffix = "")
{
  PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming* task =
   PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming::AddTaskJetDynamicalGrooming(
    njetsBase, njetsUS, njetsTrue, njetsPartLevel, R, nrhoBase, ntracks, ntracksUS, ntracksPartLevel, nclusters,
    ntracksTrue, type, CentEst, pSel, jetShapeType, jetShapeSub, jetSelection, minpTHTrigger, maxpTHTrigger, acut,
    derivSubtrOrder, suffix);
  return task;
}
