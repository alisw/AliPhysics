AliAnalysisTaskSEDplus *AddTaskDplus(Int_t system=0/*0=pp,1=PbPb*/,
				     Float_t minC=0, Float_t maxC=100,
				     Bool_t storeNtuple=kFALSE,
				     Bool_t readMC=kFALSE,
				     TString finDirname="Loose",
				     TString filename="DplustoKpipiCuts.root",
				     TString finAnObjname="AnalysisCuts", TString finProdObjname="ProduCuts")
{
  //                                                                                                                                    
  // Test macro for the AliAnalysisTaskSE for D+ candidates 

  //Invariant mass histogram and                                                 
  // association with MC truth (using MC info in AOD)                                                                                   
  //  R. Bala, bala@to.infn.it                                                                                                                                  
  // Get the pointer to the existing analysis manager via the static access method.                                                     
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDplus", "No analysis manager to connect to.");
  }

  Bool_t stdcuts=kFALSE;
  TFile* filecuts=new TFile(filename.Data());
  if(!filecuts->IsOpen()){
    cout<<"Input file not found: using standard cuts"<<endl;
    stdcuts=kTRUE;
  }
  
  
  //Analysis Task

  
  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  if(stdcuts) {
    if(system==0) analysiscuts->SetStandardCutsPP2010();
    else if(system==1){
      analysiscuts->SetStandardCutsPbPb2010();
      analysiscuts->SetMinCentrality(minC);
      analysiscuts->SetMaxCentrality(maxC);
      //      analysiscuts->SetUseAOD049(kTRUE);
      analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M);
    }
  }
  else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(finAnObjname);

  AliRDHFCutsDplustoKpipi* prodcuts=new AliRDHFCutsDplustoKpipi();
  if(stdcuts)  {
    if(system==0) prodcuts->SetStandardCutsPP2010();
    else if(system==1) {
      prodcuts->SetStandardCutsPbPb2010();
      prodcuts->SetMinCentrality(minC);
      prodcuts->SetMaxCentrality(maxC);
      //      prodcuts->SetUseAOD049(kTRUE);
      prodcuts->SetUseCentrality(AliRDHFCuts::kCentV0M);
    }
  }
  else prodcuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(finProdObjname);
  
  //AliRDHFCutsDplustoKpipi *prodcuts = (AliRDHFCutsDplustoKpipi*)fileCuts->Get("ProdCuts");
  //AliRDHFCutsDplustoKpipi *analysiscuts = (AliRDHFCutsDplustoKpipi*)fileCuts->Get("AnalysisCuts");

  
  AliAnalysisTaskSEDplus *dplusTask = new AliAnalysisTaskSEDplus("DplusAnalysis",analysiscuts,prodcuts,storeNtuple);
  dplusTask->SetReadMC(readMC);
  dplusTask->SetDoLikeSign(kFALSE);
  //  dplusTask->SetUseTPCpid(kTRUE);
  //dplusTask->SetUseTOFpid(kTRUE);
  dplusTask->SetDebugLevel(0);
  dplusTask->SetMassLimits(0.2);
  dplusTask->SetUseBit(kTRUE);

  if (system==0) dplusTask->SetDoImpactParameterHistos(kTRUE);

  mgr->AddTask(dplusTask);
  
  // Create containers for input/output 

  TString inname = "cinputDplus";
  TString outname = "coutputDplus";
  TString cutsname = "coutputDplusCuts";
  TString normname = "coutputDplusNorm";
  TString ntuplename = "coutputDplus2";
  inname += finDirname.Data();
  outname += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  ntuplename += finDirname.Data();
  TString centr=Form("%.0f%.0f",analysiscuts->GetMinCentrality(),analysiscuts->GetMaxCentrality());
  inname += centr;
  outname += centr;
  cutsname += centr;
  normname += centr;
  ntuplename += centr;


  AliAnalysisDataContainer *cinputDplus = mgr->CreateContainer(inname,TChain::Class(),
							       AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassDplus";
  
  AliAnalysisDataContainer *coutputDplusCuts = mgr->CreateContainer(cutsname,TList::Class(),
								    AliAnalysisManager::kOutputContainer,
								    outputfile.Data());
  
  AliAnalysisDataContainer *coutputDplus = mgr->CreateContainer(outname,TList::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  AliAnalysisDataContainer *coutputDplusNorm = mgr->CreateContainer(normname,AliNormalizationCounter::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  
  if(storeNtuple){
    AliAnalysisDataContainer *coutputDplus2 = mgr->CreateContainer(ntuplename,TNtuple::Class(),
								   AliAnalysisManager::kOutputContainer,
								   "InvMassDplus_nt1.root");
    
    coutputDplus2->SetSpecialOutput();
  }
  mgr->ConnectInput(dplusTask,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(dplusTask,1,coutputDplus);
  
  mgr->ConnectOutput(dplusTask,2,coutputDplusCuts);

  mgr->ConnectOutput(dplusTask,3,coutputDplusNorm);  
  if(storeNtuple){
    mgr->ConnectOutput(dplusTask,4,coutputDplus2);
  }
  return dplusTask;
}
