AliAnalysisTaskRecursiveSoftDrop* AddTaskRecursiveSoftDrop(
							   const char * njetsData, //data jets
							   const char * njetsTrue, //Pyhthia Particle Level
							   const char * njetsDet,
							   const char * njetsHybridUs,
							   const char * njetsHybridS,
							   const Double_t R,
							   const char * nrhoBase,
							   const char * ntracksData,
							   const char * ntracksTrue,
							   const char * ntracksDet,
							   const char * ntracksHybridUs,
							   const char * ntracksHybridS,
							   const char *type,
							   const char *CentEst,
							   Int_t       pSel,
							   TString     trigClass,
							   TString     kEmcalTriggers,
							   TString     tag,
							   AliAnalysisTaskRecursiveSoftDrop::JetShapeSub jetShapeSub = AliAnalysisTaskRecursiveSoftDrop::kNoSub,
							   AliAnalysisTaskRecursiveSoftDrop::JetType fjetType = AliAnalysisTaskRecursiveSoftDrop::kData
							   ) {

  AliAnalysisTaskRecursiveSoftDrop *task = AliAnalysisTaskRecursiveSoftDrop::AddTaskRecursiveSoftDrop(njetsData,njetsTrue,njetsDet,njetsHybridUs,njetsHybridS,R,nrhoBase,ntracksData,ntracksTrue,ntracksDet,ntracksHybridUs,ntracksHybridS,type,CentEst,pSel,trigClass,kEmcalTriggers,tag,jetShapeSub,fjetType);
  return task;

}
