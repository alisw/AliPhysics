Int_t AddPWG3MuonTrain(Int_t iESDAnalysis,
                       Int_t iAODAnalyis,
                       Int_t addMuonDistributions,
		       Int_t addSingleMuonAnalysis,
		       Int_t addMuonHFAnalysis) {

  // Analysis wagons for PWG3Muon (Roberta)

  TString taskName="",loadMacroPath="$ALICE_ROOT/PWG3/muon/";
  Int_t ntasks=0;
  
  if(addMuonDistributions) {
    taskName="AddTaskMuonDistributions.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    Bool_t doInvMassFit = kTRUE;
    if(iESDAnalysis){
      AliAnalysisTaskMuonDistributions *esdmuondistributionstask = AddTaskMuonDistributions("ESD",doInvMassFit);	
    } else if(iAODAnalysis){
      AliAnalysisTaskMuonDistributions *aodmuondistributionstask = AddTaskMuonDistributions("AOD",doInvMassFit);	
    }
    ntasks++;
  }

  if(addSingleMuonAnalysis) {
    taskName="AddTaskSingleMuonAnalysis.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSingleMu *singlemutask = AddTaskSingleMuonAnalysis();	
    ntasks++;
  }

  if(addMuonHFAnalysis) {
    taskName="AddTaskMuonsHF.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    Int_t runMode = 0;
    Bool_t isAOD;
    if(iAODAnalysis==1) isAOD=kTRUE;
    else isAOD=kFALSE;
    Bool_t isTree = kFALSE;
    Bool_t isMC   = kFALSE;
    AliAnalysisTaskSEMuonsHF *muonhftask = AddTaskMuonsHF(runMode, isMC, isTree);	
    ntasks++;
  }

  return ntasks;
}
