AliAnalysisTaskSELambdac *AddTaskLambdac(TString finname,Bool_t storeNtuple,Bool_t readMC,Bool_t MCPid,Bool_t realPid,Bool_t resPid,Bool_t useKF,
					 Bool_t fillVarHists=kFALSE, Bool_t priorsHists=kFALSE, Bool_t multiplicityHists=kFALSE)
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
  lambdacTask->SetFillVarHists(fillVarHists);
  lambdacTask->SetPriorsHists(priorsHists);
  lambdacTask->SetMultiplicityHists(multiplicityHists);

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
  mgr->ConnectInput(lambdacTask,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutputLambdacCuts = mgr->CreateContainer("coutputLambdacCuts",TList::Class(),
								      AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,2,coutputLambdacCuts);

  AliAnalysisDataContainer *coutputLambdac = mgr->CreateContainer("coutputLambdac",TList::Class(),
								  AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,1,coutputLambdac);

  AliAnalysisDataContainer *coutputLambdacMC = mgr->CreateContainer("coutputLambdacMC",TList::Class(),
								    AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,3,coutputLambdacMC);

  AliAnalysisDataContainer *coutputLambdacNev = mgr->CreateContainer("coutputLambdacNev",TH1F::Class(),
								     AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,4,coutputLambdacNev);

  AliAnalysisDataContainer *coutputAPriori = mgr->CreateContainer("coutputAPriori",TList::Class(),
								  AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,5,coutputAPriori);
  AliAnalysisDataContainer *coutputMultiplicity = mgr->CreateContainer("coutputMultiplicity",TList::Class(),
								       AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,6,coutputMultiplicity);

  if (storeNtuple) {
    AliAnalysisDataContainer *coutputLambdac2 = mgr->CreateContainer("coutputLambdac2",TNtuple::Class(),
								     AliAnalysisManager::kOutputContainer,"InvMassLambdac_nt1.root");
    coutputLambdac2->SetSpecialOutput();
    mgr->ConnectOutput(lambdacTask,7,coutputLambdac2);
  }


  return lambdacTask;
}
