AliAnalysisTaskSELambdac *AddTaskLambdac(TString finname,Bool_t storeNtuple,Bool_t readMC,Bool_t MCPid,Bool_t realPid,Bool_t resPid,Bool_t useKF)
{
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLambdac", "No analysis manager to connect to.");
    return NULL;
  }


  Bool_t stdcuts=kFALSE;
  TFile* filecuts=new TFile(finname.Data());
  if(!filecuts->IsOpen()){
   cout<<"Input file not found: using std cut object"<<endl;
   stdcuts=kTRUE;
  }
  AliRDHFCutsLctopKpi* prodcuts=new AliRDHFCutsLctopKpi();
  if(stdcuts) prodcuts->SetStandardCutsPP2010();
  else   prodcuts = (AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiProdCuts");
  prodcuts->SetName("LctopKpiProdCuts");
  prodcuts->SetMinPtCandidate(-1.);
  prodcuts->SetMaxPtCandidate(10000.);

  AliRDHFCutsLctopKpi *analysiscuts = new AliRDHFCutsLctopKpi();
  if(stdcuts) analysiscuts->SetStandardCutsPP2010();
  else analysiscuts = (AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiAnalysisCuts");
  analysiscuts->SetName("LctopKpiAnalysisCuts");
  analysiscuts->SetMinPtCandidate(-1.);
  analysiscuts->SetMaxPtCandidate(10000.);

  // Aanalysis task                                                                                                                     
  AliAnalysisTaskSELambdac *lambdacTask = new AliAnalysisTaskSELambdac("LambdacAnalysis",storeNtuple,analysiscuts,prodcuts);
  lambdacTask->SetReadMC(readMC);
  if(MCPid) lambdacTask->SetMCPid();
  if(resPid) lambdacTask->SetResonantPid();
  if(realPid) lambdacTask->SetRealPid();
  lambdacTask->SetDebugLevel(0);
  if(useKF) {
    lambdacTask->SetUseKF();
     Float_t cuts[10]={0.1,0.1,1.5,0.5,0.1,1.5,0.5,0.1,1.5,0.5};
      lambdacTask->SetCutsKF(cuts);
       }
  mgr->AddTask(lambdacTask);

  //                                                                                                                                    
  // Create containers for input/output                                                                                                 
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassLambdac";
  AliAnalysisDataContainer *cinputLambdac = mgr->CreateContainer("cinputLambdac",TChain::Class(),
                                                          AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutputLambdacCuts = mgr->CreateContainer("coutputLambdacCuts",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  AliAnalysisDataContainer *coutputLambdac = mgr->CreateContainer("coutputLambdac",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  if(storeNtuple){
    AliAnalysisDataContainer *coutputLambdac2 = mgr->CreateContainer("coutputLambdac2",TNtuple::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								 "InvMassLambdac_nt1.root");

    coutputLambdac2->SetSpecialOutput();
  }

  mgr->ConnectInput(lambdacTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(lambdacTask,1,coutputLambdac);
  mgr->ConnectOutput(lambdacTask,2,coutputLambdacCuts);
  
  if(storeNtuple){
    mgr->ConnectOutput(lambdacTask,3,coutputLambdac2);
  }
  return lambdacTask;
}
