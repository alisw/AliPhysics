AliAnalysisTaskSEDs *AddTaskDs(Int_t system=0/*0=pp,1=PbPb*/,
			       Int_t storeNtuple=0,Bool_t storeNsparse=kFALSE,Bool_t storeNsparseDplus=kFALSE,Bool_t readMC=kFALSE,
			       TString filename="", TString postname="", Bool_t doCutVarHistos = kFALSE, Int_t AODProtection = 1,
			       Bool_t fillNTrklAxis = kFALSE, Int_t fillCentrAxis = 0,
			       Bool_t useRotBkg=kFALSE, Bool_t useBkgFromPhiSB=kFALSE, Bool_t useCutV0multTPCout=kFALSE,
			       Bool_t storeNsparseImpPar = kFALSE)
{
  //
  // Test macro for the AliAnalysisTaskSE for Ds candidates
    
  //Invariant mass histogram and
  // association with MC truth (using MC info in AOD)
  // Origin: R. Bala, bala@to.infn.it
  // Modified for Ds meson: G.M. Innocenti innocent@to.infn.it
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDs", "No analysis manager to connect to.");
  }
    
  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( filename.EqualTo("") ) {
    stdcuts=kTRUE;
  } else {
    filecuts=TFile::Open(filename.Data());
    if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
      AliFatal("Cut object not found: analysis will not start!\n");
    }
    else printf("Cut object correctly found\n");
  }
    
  //Analysis Task
    
  AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    
  if(stdcuts) {
    if(system==0) {
      printf("Cut object not found: standard pp cut object used\n");
      analysiscuts->SetStandardCutsPP2010();
    }
    else AliFatal("Standard cut object not available for PbPb: analysis will not start!\n");
  }
  else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get("AnalysisCuts");
    
  AliAnalysisTaskSEDs *dsTask = new AliAnalysisTaskSEDs("DsAnalysis",analysiscuts,storeNtuple);
    
  dsTask->SetReadMC(readMC);
  //dsTask->SetDoLikeSign(kTRUE);
  //dsTask->SetUseTPCpid(kTRUE);
  //dsTask->SetUseTOFpid(kTRUE);
  dsTask->SetDebugLevel(0);
  dsTask->SetUseSelectionBit(kTRUE);
  //dsTask->SetMassLimits(0.2);
  dsTask->SetFillNSparse(storeNsparse);
  dsTask->SetFillNSparseDplus(storeNsparseDplus);
  dsTask->SetFillNSparseImpPar(storeNsparseImpPar);
  dsTask->SetAODMismatchProtection(AODProtection);
  dsTask->SetDoCutVarHistos(doCutVarHistos);
  dsTask->SetUseRotBkg(useRotBkg);
  dsTask->SetUseBkgFromPhiSB(useBkgFromPhiSB);
  dsTask->SetUseCutV0multVsTPCout(useCutV0multTPCout);
  dsTask->SetSystem(system);
  dsTask->SetFillTracklets(fillNTrklAxis);
  dsTask->SetFillCentralityAxis(fillCentrAxis);
  mgr->AddTask(dsTask);
    
  // Create containers for input/output
  TString name="cinputDs";
  name+=postname;
  AliAnalysisDataContainer *cinputDs = mgr->CreateContainer(name,TChain::Class(),
							    AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassDs";
  outputfile+=postname;
    
  name="coutputDsCuts"; name+=postname;
  AliAnalysisDataContainer *coutputDsCuts = mgr->CreateContainer(name,TList::Class(),
								 AliAnalysisManager::kOutputContainer,
								 outputfile.Data());
    
  name="coutputDs"; name+=postname;
  AliAnalysisDataContainer *coutputDs = mgr->CreateContainer(name,TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputfile.Data());
  name="coutputDsNorm"; name+=postname;
  AliAnalysisDataContainer *coutputDsNorm = mgr->CreateContainer(name,AliNormalizationCounter::Class(),
								 AliAnalysisManager::kOutputContainer,
								 outputfile.Data());
    
  name="coutputDs2"; name+=postname;
  if(storeNtuple){
    AliAnalysisDataContainer *coutputDs2 = mgr->CreateContainer(name,TNtuple::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
    coutputDs2->SetSpecialOutput();
  }
    
  mgr->ConnectInput(dsTask,0,mgr->GetCommonInputContainer());
    
  mgr->ConnectOutput(dsTask,1,coutputDs);
    
  mgr->ConnectOutput(dsTask,2,coutputDsCuts);
    
  mgr->ConnectOutput(dsTask,3,coutputDsNorm);  
    
  if(storeNtuple){
    mgr->ConnectOutput(dsTask,4,coutputDs2);
  }
    
  return dsTask;
}
