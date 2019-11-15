AliAnalysisTaskSEDs *AddTaskDs(Int_t system = AliAnalysisTaskSEDs::kpp, Bool_t readMC = kFALSE, Int_t AODProtection = 1, 
                               Bool_t storeNsparse = kFALSE,  TString filename = "", TString postname = "", 
                               Bool_t applyML = kFALSE, TString confFileML = "", TString cutObjName = "AnalysisCuts", 
                               Bool_t storeNsparseDplus =  kFALSE, Bool_t doCutVarHistos = kFALSE,   Bool_t storeNsparseImpPar = kFALSE)
{
  // AddTask for AliAnalysisTaskSEDs
  // Origin: R. Bala, bala@to.infn.it
  // Modified for Ds meson: G.M. Innocenti innocent@to.infn.it
  // Mantainers: F. Catalano, fabio.catalano@cern.ch
  //             F. Grosa, fabrizio.grosa@cern.ch
  //==============================================================================

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    ::Error("AddTaskDs", "No analysis manager to connect to.");

  TFile* filecuts = TFile::Open(filename.Data());
  if(!filecuts || (filecuts && !filecuts->IsOpen()))
    ::Fatal("AddTaskDs", "Cut file not found on Grid: analysis will not start!\n");
  else
    printf("Cut file correctly found\n");

  //Analysis Task
  AliRDHFCutsDstoKKpi* analysiscuts = new AliRDHFCutsDstoKKpi();
  analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(cutObjName.Data());
  if (!analysiscuts){
    ::Fatal("AddTaskDs", "Cut object not found in cutfile");
	  return NULL;
	}

  AliAnalysisTaskSEDs *dsTask = new AliAnalysisTaskSEDs("DsAnalysis",analysiscuts);
  dsTask->SetDebugLevel(0);
  dsTask->SetUseSelectionBit(kTRUE);
  dsTask->SetSystem(system);
  dsTask->SetReadMC(readMC);
  dsTask->SetAODMismatchProtection(AODProtection);
  dsTask->SetFillNSparse(storeNsparse);
  dsTask->SetDoMLApplication(applyML);
  if(applyML)
    dsTask->SetMLConfigFile(confFileML);
  dsTask->SetFillNSparseDplus(storeNsparseDplus);
  dsTask->SetDoCutVarHistos(doCutVarHistos);
  dsTask->SetFillNSparseImpPar(storeNsparseImpPar);
  if(system == AliAnalysisTaskSEDs::kPbPb || system == AliAnalysisTaskSEDs::kUpgr)
    dsTask->SetKeepOnlyBkgFromHIJING(kTRUE);
  
  mgr->AddTask(dsTask);

  // Create containers for input/output
  TString name = Form("cinputDs%s", postname);
  AliAnalysisDataContainer *cinputDs = mgr->CreateContainer(name, TChain::Class(), AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWGHF_D2H_InvMassDs";
  outputfile += postname;

  name = Form("coutputDsCuts%s", postname);
  AliAnalysisDataContainer *coutputDsCuts = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer,
								                                                 outputfile.Data());

  name =  Form("coutputDs%s", postname);
  AliAnalysisDataContainer *coutputDs = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer,
							                                               outputfile.Data());
  
  name =  Form("coutputDsNorm%s", postname);
  AliAnalysisDataContainer *coutputDsNorm = mgr->CreateContainer(name, AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer,
								                                                 outputfile.Data());

  //coutputDs2->SetSpecialOutput();

  mgr->ConnectInput(dsTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(dsTask, 1, coutputDs);
  mgr->ConnectOutput(dsTask, 2, coutputDsCuts);
  mgr->ConnectOutput(dsTask, 3, coutputDsNorm);

  return dsTask;
}
