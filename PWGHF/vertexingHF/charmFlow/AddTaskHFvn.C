AliAnalysisTaskSEHFvn *AddTaskHFvn(Int_t harm, TString filename="alien:///alice/cern.ch/user/a/abarbano/DstoKKpiCutsCentrality20to50_strongPID.root",AliAnalysisTaskSEHFvn::DecChannel decCh=AliAnalysisTaskSEHFvn::kDstoKKpi,TString cutsobjname="AnalysisCuts", Bool_t readMC=kFALSE, TString suffix="", AliAnalysisTaskSEHFvn::EventPlaneMeth flagep=AliAnalysisTaskSEHFvn::kVZERO/*kTPC,kTPCVZERO,kVZEROA,kVZEROC*/,Float_t minC=20.,Float_t maxC=50., Bool_t useNewQnFw=kTRUE, AliAnalysisTaskSEHFvn::FlowMethod meth=AliAnalysisTaskSEHFvn::kEP/*kSP,kEvShape*/, TString normMethod="QoverM"/*"QoverQlength","QoverSqrtM"*/,AliAnalysisTaskSEHFvn::q2Method q2meth=AliAnalysisTaskSEHFvn::kq2TPC/*kq2PosTPC,kq2NegTPC,kq2VZERO,kq2VZEROA,kq2VZEROC}*/, Int_t useAODProtection=1)
{
  //
  // Test macro for the AliAnalysisTaskSE for  D
  // mesons v2 analysis with event plane method

  // Get the pointer to the existing analysis manager via the static access method.
  //============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFvn", "No analysis manager to connect to.");
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
  
  if(decCh==AliAnalysisTaskSEHFvn::kDplustoKpipi) {
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDplustoKpipi();
      analysiscuts->SetStandardCutsPbPb2011();
    } else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Dplus");
    pdgmes=411;
  } else if(decCh==AliAnalysisTaskSEHFvn::kD0toKpi) {
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPbPb2011();
    } else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Dzero");
    pdgmes=421;
  } else if(decCh==AliAnalysisTaskSEHFvn::kDstartoKpipi) {
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      analysiscuts->SetStandardCutsPbPb2011();
    } else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Dstar");
    pdgmes=413;
  }
  else if(decCh==AliAnalysisTaskSEHFvn::kDstoKKpi) {
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDstoKKpi();
      //analysiscuts->SetStandardCutsPbPb2011(); //to be implemented in AliRDHFCutsDstoKKpi
    } else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(cutsobjname);
    suffix.Prepend("Ds");
    pdgmes=431;
  }
  else if(decCh==AliAnalysisTaskSEHFvn::kD0toKpiFromDstar) {
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPbPb2011();
    } else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
    suffix.Prepend("DzeroFromDstar");
    pdgmes=421;
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
  AliAnalysisTaskSEHFvn *v2Task = new AliAnalysisTaskSEHFvn("HFvnAnalysis",analysiscuts,decCh);
  v2Task->SetHarmonic(harm);
  v2Task->SetReadMC(readMC);
  v2Task->SetEtaGapFeatureForEventplaneFromTracks(kFALSE);
  if(decCh == AliAnalysisTaskSEHFvn::kDstartoKpipi) {
    v2Task->SetNMassBins(200);
  } else if(decCh == AliAnalysisTaskSEHFvn::kDplustoKpipi || decCh == AliAnalysisTaskSEHFvn::kD0toKpi || decCh == AliAnalysisTaskSEHFvn::kD0toKpiFromDstar) {
    v2Task->SetNMassBins(104);
    v2Task->SetMassLimits(0.2,pdgmes);
  } else if(decCh == AliAnalysisTaskSEHFvn::kDstoKKpi) {
    v2Task->SetNMassBins(200);  // to be decided
  }
  v2Task->SetMinCentrality(minC);
  v2Task->SetMaxCentrality(maxC);
  v2Task->SetDebugLevel(0);
  v2Task->SetEventPlaneMethod(flagep);
  if(flagep==AliAnalysisTaskSEHFvn::kTPC) {
    suffix+="TPC";
  } else if(flagep==AliAnalysisTaskSEHFvn::kTPCVZERO) {
    suffix+="TPCVZERO";
  } else if(flagep==AliAnalysisTaskSEHFvn::kVZERO) {
      suffix+="VZERO";
  } else if(flagep==AliAnalysisTaskSEHFvn::kVZEROA) {
    suffix+="VZEROA";
  } else if(flagep==AliAnalysisTaskSEHFvn::kVZEROC) {
    suffix+="VZEROC";
  } else if(flagep==AliAnalysisTaskSEHFvn::kPosTPCVZERO) {
    suffix+="POSTPCVZERO";
  } else if(flagep==AliAnalysisTaskSEHFvn::kNegTPCVZERO) {
    suffix+="NEGTPCVZERO";
  }
  v2Task->SetFlowMethod(meth);
  if(meth==AliAnalysisTaskSEHFvn::kEP) {
    suffix+="_EP";
  } else if(meth==AliAnalysisTaskSEHFvn::kSP) {
    suffix+="_SP";
  } else if(meth==AliAnalysisTaskSEHFvn::kEvShape) {
    suffix+="_EvShape";
  }
  v2Task->SetNormMethod(normMethod);
  v2Task->SetUseNewQnCorrFw(useNewQnFw);
  v2Task->SetAODMismatchProtection(useAODProtection);
  if(meth==AliAnalysisTaskSEHFvn::kEvShape) {v2Task->Setq2Method(q2meth);}
  mgr->AddTask(v2Task);
    
  // Create containers for input/output
    
  TString contname=Form("cinputv2%s",suffix.Data());
  AliAnalysisDataContainer *cinputv2 = mgr->CreateContainer(contname.Data(),TChain::Class(),AliAnalysisManager::kInputContainer);
    
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  TString outputhistos = Form("%s:PWGHF_D2H_HFvn_%s",outputfile.Data(),suffix.Data());
    
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
