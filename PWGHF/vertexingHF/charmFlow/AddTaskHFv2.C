AliAnalysisTaskSEHFv2 *AddTaskHFv2(TString filename="DplustoKpipiCutsPbPb.root",AliAnalysisTaskSEHFv2::DecChannel decCh=AliAnalysisTaskSEHFv2::kDplustoKpipi,TString cutsobjname="AnalysisCuts", Bool_t readMC=kFALSE,TString name="",Int_t flagep=1 /*0=tracks,1=V0*/)
{
  //
  // Test macro for the AliAnalysisTaskSE for  D 
  // mesons v2 analysis with event plane method
  // Authors: Chiara Bianchin, cbianchi@pd.infn.it, 
  //          Robert Grajcarek, grajcarek@physi.uni-heidelberg.de
  //          Giacomo Ortona, ortona@to.infn.it,
  //          Carlos Perez Lara, carlos.eugenio.perez.lara@cern.ch
  //          Francesco Prino, prino@to.infn.it
  // Get the pointer to the existing analysis manager via the static access method.
  //============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFv2", "No analysis manager to connect to.");
    return NULL;
  }
  Bool_t stdcuts=kFALSE;
  TFile* filecuts=new TFile(filename.Data());
  if(!filecuts->IsOpen()){
    cout<<"Input file not found:  using std cut object"<<endl;
    stdcuts=kTRUE;
  }
  
  AliRDHFCuts *analysiscuts=0x0;
  TString suffix="";

  //Analysis cuts
  if(decCh==AliAnalysisTaskSEHFv2::kDplustoKpipi){
    analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
  }else if(decCh==AliAnalysisTaskSEHFv2::kD0toKpi){
    analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
  }else if(decCh==AliAnalysisTaskSEHFv2::kDstartoKpipi){
    analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(cutsobjname);
  }

  if(!analysiscuts){
    AliFatal("Specific AliRDHFCuts not found");
  }

  analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M);
  analysiscuts->SetRemoveDaughtersFromPrim(kFALSE);
  analysiscuts->SetTriggerClass("");
  analysiscuts->ResetMaskAndEnableMBTrigger();
  analysiscuts->EnableCentralTrigger();
  analysiscuts->EnableSemiCentralTrigger();
  analysiscuts->SetMinCentrality(30.);
  analysiscuts->SetMaxCentrality(80.);
  ((AliRDHFCutsDplustoKpipi*)analysiscuts)->SetMinPtCandidate(3.);
  ((AliRDHFCutsDplustoKpipi*)analysiscuts)->SetUseImpParProdCorrCut(kFALSE);
  AliAODPidHF *pid = analysiscuts->GetPidHF();
  pid->SetOldPid(kFALSE);
  analysiscuts->SetPidHF(pid);
  // Analysis task    
  AliAnalysisTaskSEHFv2 *v2Task = new AliAnalysisTaskSEHFv2("HFv2Analysis",analysiscuts,decCh);
  v2Task->SetReadMC(kFALSE);
  v2Task->SetEtaGapFeatureForEventplaneFromTracks(kFALSE);
  v2Task->SetNMassBins(104);
  v2Task->SetMassLimits(0.2,411);
  v2Task->SetDebugLevel(0);
  v2Task->SetV0EventPlaneOrder(2);
  v2Task->SetTPCEP();//SetVZEROEPOnly();

  mgr->AddTask(v2Task);

  // Create containers for input/output

  TString contname=Form("cinputv2%s",suffix.Data());
  AliAnalysisDataContainer *cinputv2 = mgr->CreateContainer(contname.Data(),TChain::Class(),AliAnalysisManager::kInputContainer);

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  TString outputhistos = outputfile += ":PWG3_D2H_HFv2"; 

  contname=Form("hEventsInfo%s",suffix.Data());
  AliAnalysisDataContainer *coutputstat = mgr->CreateContainer(contname.Data(),TH1F::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  contname=Form("coutputv2%s",suffix.Data());
  AliAnalysisDataContainer *coutputv2 = mgr->CreateContainer(contname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  contname=Form("cutobj%s",suffix.Data());
  AliAnalysisDataContainer *cutobj = mgr->CreateContainer(contname.Data(),AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  mgr->ConnectInput(v2Task,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(v2Task,1,coutputstat);
  
  mgr->ConnectOutput(v2Task,2,coutputv2);

  mgr->ConnectOutput(v2Task,3,cutobj);
 
  return v2Task;
}
