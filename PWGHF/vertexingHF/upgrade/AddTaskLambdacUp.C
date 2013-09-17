AliAnalysisTaskSE *AddTaskLambdacUp(TString finname,Bool_t storeNtuple,Bool_t readMC,Bool_t MCPid,Bool_t realPid,Bool_t resPid,Bool_t useKF,
				  Bool_t fillVarHists=kFALSE,TString postname="")
{
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLambdac", "No analysis manager to connect to.");
    return NULL;
  }


  TString taskName=("AliAnalysisTaskSELambdacUp.cxx+");
  gROOT->LoadMacro(taskName.Data());
  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( finname.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
      filecuts=TFile::Open(finname.Data());
      if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
	AliFatal("Input file not found : check your cut object");
      }
  }
  AliRDHFCutsLctopKpi* prodcuts=new AliRDHFCutsLctopKpi();
  prodcuts = (AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiProdCuts");
  prodcuts->SetName("LctopKpiProdCuts");
  prodcuts->SetMinPtCandidate(-1.);
  prodcuts->SetMaxPtCandidate(10000.);

  AliRDHFCutsLctopKpi *analysiscuts = new AliRDHFCutsLctopKpi();
  analysiscuts = (AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiAnalysisCuts");
  analysiscuts->SetName("LctopKpiAnalysisCuts");
  analysiscuts->SetMinPtCandidate(-1.);
  analysiscuts->SetMaxPtCandidate(10000.);

  // Aanalysis task                                                                                                                     
  AliAnalysisTaskSELambdacUp *lambdacTask = new AliAnalysisTaskSELambdacUp("LambdacAnalysis",storeNtuple,analysiscuts,prodcuts);
  lambdacTask->SetReadMC(readMC);
  if(MCPid) lambdacTask->SetMCPid();
  if(resPid) lambdacTask->SetResonantPid();
  if(realPid) lambdacTask->SetRealPid();
  lambdacTask->SetFillVarHists(fillVarHists);
  
  lambdacTask->SetAnalysis(kTRUE);

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

  TString finDirname="pp";
  TString inname = "cinputLc";
  TString outname = "coutputLc";
  TString cutsname = "coutputLcCuts";
  TString normname = "coutputLcNorm";
  TString ntuplename = "coutputLc2";
  TString nev2 = "coutputNev";
  TString outname2 = "coutputLambdacMC";
  
  inname += finDirname.Data();
  outname += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  ntuplename += finDirname.Data();
  nev2 += finDirname.Data();
  outname2 += finDirname.Data();
  
  inname +=  postname.Data();
  outname +=  postname.Data();
  cutsname +=  postname.Data();
  normname += postname.Data();
  ntuplename +=  postname.Data();
  nev2 +=  postname.Data();
  outname2 +=  postname.Data();
  

  AliAnalysisDataContainer *cinputLambdac = mgr->CreateContainer(inname,TChain::Class(),
								 AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(lambdacTask,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutputLambdacCuts = mgr->CreateContainer(cutsname,TList::Class(),
								      AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,2,coutputLambdacCuts);

  AliAnalysisDataContainer *coutputLambdac = mgr->CreateContainer(outname,TList::Class(),
								  AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,1,coutputLambdac);

  AliAnalysisDataContainer *coutputLambdacMC = mgr->CreateContainer(outname2,TList::Class(),
								    AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,3,coutputLambdacMC);

  AliAnalysisDataContainer *coutputLambdacNev = mgr->CreateContainer(nev2,TH1F::Class(),
								     AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,4,coutputLambdacNev);

 

  AliAnalysisDataContainer *coutputLambdacNorm = mgr->CreateContainer(normname,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

 mgr->ConnectOutput(lambdacTask,5,coutputLambdacNorm);

  if (storeNtuple) {
    AliAnalysisDataContainer *coutputLambdac2 = mgr->CreateContainer(ntuplename,TNtuple::Class(),
								     AliAnalysisManager::kOutputContainer,"InvMassLambdac_nt1.root");
    coutputLambdac2->SetSpecialOutput();
    mgr->ConnectOutput(lambdacTask,6,coutputLambdac2);
  }


  return lambdacTask;
}
