AliAnalysisTempFluc*  AddTaskTempFluc(TString  analysisType="mc",//pp,pPb,PbPb
				  TString  dataType="AOD",
					       TString CentralityMethod="V0M",//V0A/V0M          
					       Int_t FilterBit=786,//768,16,32
					       Float_t minPt=0.15,
					       Float_t maxPt=3.15,
					       Float_t mineta=-0.5,
					       Float_t maxeta=0.5,
					       Float_t minrap=-0.5,
					       Float_t maxrap=0.5,
                           Bool_t  isCorrect_efficiency=kFALSE,
					       const char* outputFileName = "Result_mc_fb786_ns3.root",
					       const char* containerName = "TempFluc",
                           const char* listname = "coutput"
                                  )


{

 // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTempFluc", "No analysis manager to connect to.");
    return NULL;
  }  
 
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskTempFluc", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType();

   AliAnalysisTempFluc *task = new AliAnalysisTempFluc(containerName);
   
   task->SetAnalysisMode(analysisType);
   task->SetDataType(dataType);
   task->SetCentralityEstimator("V0M");
   task->SetCentralityCutL(0.0);
   task->SetCentralityCutH(5.0);
   task->SetgEtaL(mineta);
   task->SetgEtaH(maxeta);
   task->SetgRapL(minrap);
   task->SetgRapH(maxrap);
   task->SetPtCutL(minPt);
   task->SetPtCutH(maxPt);
   task->SetFilterBit(FilterBit);
   task->SetCutNSigmaTPC(3.0);
   task->SetCutNSigmaTPCTOF(3.0);
   task->Correct_efficiency(isCorrect_efficiency);
   task->SetCustomBinningptL("0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.5,0.5,0.5");
   task->SetCustomBinningptH("1.0,1.2,1.5,2.0,1.0,1.2,1.5,1.0,1.2,1.5");
    
    task->SetEffcorectionfilePathName("alien:////alice/cern.ch/user/s/subasu/submit_test/efficiency_0005.root");

   mgr->AddTask(task);

 // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  if (!outputFileName)
    outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName, TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, listname));
  
  
  
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, coutput1 );
 

  return task;













}
