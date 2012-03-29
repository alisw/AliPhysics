AliAnalysisTaskSEHFv2 *AddTaskHFv2(TString filename="DplustoKpipiCutsPbPb.root",AliAnalysisTaskSEHFv2::DecChannel decCh=AliAnalysisTaskSEHFv2::kDplustoKpipi,TString cutsobjname="AnalysisCuts", Bool_t readMC=kFALSE, TString suffix="", Int_t flagep=0 /*0=tracks,1=V0*/)
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
  TFile* filecuts;
  if( filename.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
    filecuts=TFile::Open(filename.Data());
    if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
      AliFatal("Input file not found : check your cut object");
    }
  }
  
  AliRDHFCuts *analysiscuts=0x0;
  Int_t pdgmes=-1;
  //Analysis cuts

  if(decCh==AliAnalysisTaskSEHFv2::kDplustoKpipi){
    if(stdcuts){
      analysiscuts = new AliRDHFCutsDplustoKpipi();
      analysiscuts->SetStandardCutsPbPb2011();
    }else  analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Dplus");
    pdgmes=411;
  }else if(decCh==AliAnalysisTaskSEHFv2::kD0toKpi){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Dzero");
    pdgmes=421;
  }else if(decCh==AliAnalysisTaskSEHFv2::kDstartoKpipi){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Dstar");
    pdgmes=413;
  }
  if(pdgmes==-1){
    AliFatal("Wrong meson setting");
  }
  if(!analysiscuts){
    AliFatal("Specific AliRDHFCuts not found");
  }

  // Analysis task    
  AliAnalysisTaskSEHFv2 *v2Task = new AliAnalysisTaskSEHFv2("HFv2Analysis",analysiscuts,decCh);
  v2Task->SetReadMC(readMC);
  v2Task->SetEtaGapFeatureForEventplaneFromTracks(kFALSE);
  v2Task->SetNMassBins(104);
  v2Task->SetMassLimits(0.2,pdgmes);
  v2Task->SetDebugLevel(0);
  v2Task->SetV0EventPlaneOrder(2);
  if(flagep==0){
    v2Task->SetTPCEP();//SetVZEROEPOnly();
    suffix+="TPC";
  }else{
    v2Task->SetVZEROEP();
    suffix+="VZERO";
  }
  mgr->AddTask(v2Task);

  // Create containers for input/output

  TString contname=Form("cinputv2%s",suffix.Data());
  AliAnalysisDataContainer *cinputv2 = mgr->CreateContainer(contname.Data(),TChain::Class(),AliAnalysisManager::kInputContainer);

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  TString outputhistos = Form("%s:PWGHF_D2H_HFv2_%s",outputfile.Data(),suffix.Data()); 

  contname=Form("hEventsInfo%s",suffix.Data());
  AliAnalysisDataContainer *coutputstat = mgr->CreateContainer(contname.Data(),TH1F::Class(),AliAnalysisManager::kOutputContainer,outputhistos.Data());

  contname=Form("coutputv2%s",suffix.Data());
  AliAnalysisDataContainer *coutputv2 = mgr->CreateContainer(contname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputhistos.Data());

  contname=Form("cutobj%s",suffix.Data());
  AliAnalysisDataContainer *cutobj = mgr->CreateContainer(contname.Data(),AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer,outputhistos.Data());

  mgr->ConnectInput(v2Task,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(v2Task,1,coutputstat);
  
  mgr->ConnectOutput(v2Task,2,coutputv2);

  mgr->ConnectOutput(v2Task,3,cutobj);
 
  return v2Task;
}
