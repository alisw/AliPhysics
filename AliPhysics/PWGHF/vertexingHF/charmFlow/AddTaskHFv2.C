AliAnalysisTaskSEHFv2 *AddTaskHFv2(TString filename="alien:///alice/cern.ch/user/a/abarbano/DstoKKpiCutsCentrality20to50_strongPID.root",AliAnalysisTaskSEHFv2::DecChannel decCh=AliAnalysisTaskSEHFv2::kDstoKKpi,TString cutsobjname="AnalysisCuts", Bool_t readMC=kFALSE, TString suffix="", AliAnalysisTaskSEHFv2::EventPlaneMeth flagep=AliAnalysisTaskSEHFv2::kVZERO/*kTPC,kTPCVZERO,kVZEROA,kVZEROC*/,Float_t minC=20.,Float_t maxC=50., Bool_t useNewQnFw=kTRUE, AliAnalysisTaskSEHFv2::FlowMethod meth=AliAnalysisTaskSEHFv2::kEP/*kSP,kEvShape*/, TString normMethod="QoverM"/*"QoverQlength","QoverSqrtM"*/,AliAnalysisTaskSEHFv2::q2Method q2meth=AliAnalysisTaskSEHFv2::kq2TPC/*kq2PosTPC,kq2NegTPC,kq2VZERO,kq2VZEROA,kq2VZEROC}*/, Int_t useAODProtection=1)
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
      Printf("FATAL: Input file not found : check your cut object");
      return NULL;
    }
  }
  
  AliRDHFCuts *analysiscuts=0x0;
  Int_t pdgmes=-1;
  //Analysis cuts
  
  if(decCh==AliAnalysisTaskSEHFv2::kDplustoKpipi) {
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDplustoKpipi();
      analysiscuts->SetStandardCutsPbPb2011();
    } else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Dplus");
    pdgmes=411;
  } else if(decCh==AliAnalysisTaskSEHFv2::kD0toKpi) {
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPbPb2011();
    } else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Dzero");
    pdgmes=421;
  } else if(decCh==AliAnalysisTaskSEHFv2::kDstartoKpipi) {
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      analysiscuts->SetStandardCutsPbPb2011();
    } else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Dstar");
    pdgmes=413;
  }
  else if(decCh==AliAnalysisTaskSEHFv2::kDstoKKpi) {
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDstoKKpi();
      //analysiscuts->SetStandardCutsPbPb2011(); //to be implemented in AliRDHFCutsDstoKKpi
    } else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Ds");
    pdgmes=431;
  }
  if(pdgmes==-1){
    Printf("FATAL: Wrong meson setting");
    return NULL;
  }
  if(!analysiscuts){
    Printf("FATAL: Specific AliRDHFCuts not found");
    return NULL;
  }
    
  // Analysis task
  AliAnalysisTaskSEHFv2 *v2Task = new AliAnalysisTaskSEHFv2("HFv2Analysis",analysiscuts,decCh);
  v2Task->SetReadMC(readMC);
  v2Task->SetEtaGapFeatureForEventplaneFromTracks(kFALSE);
  if(decCh == AliAnalysisTaskSEHFv2::kDstartoKpipi) {
    v2Task->SetNMassBins(200);
  } else if(decCh == AliAnalysisTaskSEHFv2::kDplustoKpipi || decCh == AliAnalysisTaskSEHFv2::kD0toKpi) {
    v2Task->SetNMassBins(104);
    v2Task->SetMassLimits(0.2,pdgmes);
  } else if(decCh == AliAnalysisTaskSEHFv2::kDstoKKpi) {
    v2Task->SetNMassBins(200);  // to be decided
  }
  v2Task->SetMinCentrality(minC);
  v2Task->SetMaxCentrality(maxC);
  v2Task->SetDebugLevel(0);
  v2Task->SetV0EventPlaneOrder(2);
  v2Task->SetEventPlaneMethod(flagep);
  if(flagep==AliAnalysisTaskSEHFv2::kTPC) {
    suffix+="TPC";
  } else if(flagep==AliAnalysisTaskSEHFv2::kTPCVZERO) {
    suffix+="TPCVZERO";
  } else if(flagep==AliAnalysisTaskSEHFv2::kVZERO) {
      suffix+="VZERO";
  } else if(flagep==AliAnalysisTaskSEHFv2::kVZEROA) {
    suffix+="VZEROA";
  } else if(flagep==AliAnalysisTaskSEHFv2::kVZEROC) {
    suffix+="VZEROC";
  } else if(flagep==AliAnalysisTaskSEHFv2::kPosTPCVZERO) {
    suffix+="POSTPCVZERO";
  } else if(flagep==AliAnalysisTaskSEHFv2::kNegTPCVZERO) {
    suffix+="NEGTPCVZERO";
  }
  v2Task->SetFlowMethod(meth);
  if(meth==AliAnalysisTaskSEHFv2::kEP) {
    suffix+="_EP";
  } else if(meth==AliAnalysisTaskSEHFv2::kSP) {
    suffix+="_SP";
  } else if(meth==AliAnalysisTaskSEHFv2::kEvShape) {
    suffix+="_EvShape";
  }
  v2Task->SetNormMethod(normMethod);
  v2Task->SetUseNewQnCorrFw(useNewQnFw);
  v2Task->SetAODMismatchProtection(useAODProtection);
  if(meth==AliAnalysisTaskSEHFv2::kEvShape) {v2Task->Setq2Method(q2meth);}
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
