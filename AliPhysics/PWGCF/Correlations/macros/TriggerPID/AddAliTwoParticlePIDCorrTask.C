AliAnalysisTask*  AddAliTwoParticlePIDCorrTask(TString  SampleType="pPb",//pp,pPb,PbPb
					       TString CentralityMethod="V0A",//V0A/V0M          
					       Int_t FilterBit=768,//768,16,32
					       Float_t minPt=0.2,
					       Float_t maxPt=10.0,
					       Float_t mineta=-0.8,
					       Float_t maxeta=0.8,
					       TString  AnalysisType="AOD",//AOD/MCAOD
					       const char* outputFileName = 0,
					       const char* containerName = "TwoParticlePIDCorr",
					       const char* QAContainername = "TwoParticlePIDCorr_PIDQA",
					       const char* folderName = "PWGCF_TwoParticlePIDCorr",
					       Bool_t RequestEventPlane=kFALSE,
					       const char* EPContainername = "TwoParticlePIDCorr_EP"

)
{

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddAliTwoParticlePIDCorr", "No analysis manager to connect to.");
    return NULL;
  }  
 
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPIDCorr", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType();

   AliTwoParticlePIDCorr *PIDCorr = new AliTwoParticlePIDCorr(containerName);
   PIDCorr->SetSampleType( SampleType);
   PIDCorr->SetCentralityEstimator( CentralityMethod);
   PIDCorr->SetFilterBit(FilterBit);
   PIDCorr->SetKinematicCuts( minPt,  maxPt, mineta, maxeta);
   PIDCorr->SetAnalysisType(AnalysisType);
   
   //Trigger - Physics Selection
   // PIDCorr->SelectCollisionCandidates(AliVEvent::kINT7);

   mgr->AddTask(PIDCorr);

 // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  if (!outputFileName)
    outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName, TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(QAContainername, TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));

  if(RequestEventPlane==kTRUE) AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(EPContainername, TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));

  
  mgr->ConnectInput  (PIDCorr, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (PIDCorr, 1, coutput1 );
  mgr->ConnectOutput (PIDCorr, 2, coutput2 );
if(RequestEventPlane==kTRUE)  mgr->ConnectOutput (PIDCorr, 3, coutput3 );

  return PIDCorr;
}
