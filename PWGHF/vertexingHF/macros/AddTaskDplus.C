AliAnalysisTaskSEDplus *AddTaskDplus(Int_t system=AliAnalysisTaskSEDplus::kpp,
                                     Float_t minC=0, Float_t maxC=100,
                                     Bool_t storeTree=0,
                                     Int_t doSparse=0,/*0=cutvar=kFALSE && imppar=kFALSE, 1=cutvar=kTRUE && imppar=kFALSE*/
                                     /*2=cutvar=kFALSE && imppar=kTRUE, 3=cutvar=kTRUE && imppar=kTRUE*/
                                     Bool_t doTrackVarSparse=kFALSE,
                                     Bool_t readMC=kFALSE,
                                     TString finDirname="Loose",
                                     TString filename="",
                                     TString finAnObjname="AnalysisCuts",
                                     Int_t etaRange=0,
                                     Bool_t cutsDistr=kFALSE,
                                     Int_t trackletsmin=-1,
                                     Int_t trackletsmax=-1)
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
  TFile* filecuts;
  if( filename.EqualTo("") ) {
    stdcuts=kTRUE;
  } else {
      filecuts=TFile::Open(filename.Data());
      if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
	::Fatal("AddTaskDplus", "Input file not found : check your cut object");
      }
  }


  //Analysis Task


  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  if(stdcuts) {
    if(system==0) analysiscuts->SetStandardCutsPP2010();
    else if(system==1){
      analysiscuts->SetStandardCutsPbPb2011();
      analysiscuts->SetMinCentrality(minC);
      analysiscuts->SetMaxCentrality(maxC);
      //      analysiscuts->SetUseAOD049(kTRUE);
      analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M);
    }
  }
  else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(finAnObjname);


  AliAnalysisTaskSEDplus *dplusTask = new AliAnalysisTaskSEDplus("DplusAnalysis",analysiscuts,storeTree);
  dplusTask->SetReadMC(readMC);
  dplusTask->SetDoLikeSign(kFALSE);
  dplusTask->SetDebugLevel(0);
  dplusTask->SetMassLimits(0.2);
  dplusTask->SetUseBit(kTRUE);
  dplusTask->SetCutsDistr(cutsDistr);
  dplusTask->SetSystem(system);
  if(system == AliAnalysisTaskSEDplus::kPbPb) {
      dplusTask->SetKeepOnlyBkgFromHIJING(kTRUE);
  }
  if(doSparse==0) {
    dplusTask->SetDoCutVarsSparses(kFALSE);
    dplusTask->SetDoImpactParameterHistos(kFALSE);
  }
  else if(doSparse==1) {
    dplusTask->SetDoCutVarsSparses(kTRUE);
    dplusTask->SetDoImpactParameterHistos(kFALSE);
  }
  else if(doSparse==2) {
    dplusTask->SetDoCutVarsSparses(kFALSE);
    dplusTask->SetDoImpactParameterHistos(kTRUE);
  }
  else if(doSparse==3){
    dplusTask->SetDoCutVarsSparses(kTRUE);
    dplusTask->SetDoImpactParameterHistos(kTRUE);
  }
  else {
    cerr << "The doSparse flag can only be 0,1,2,3!" << endl;
  }

  if((doSparse==1 || doSparse==3) && readMC)
    dplusTask->SetDoMCAcceptanceHistos(kTRUE);
  if (doTrackVarSparse)
    dplusTask->SetDoTrackVarHistos(kTRUE);
  if(etaRange==1) dplusTask->SetUseOnlyPositiveEta();
  if(etaRange==-1) dplusTask->SetUseOnlyNegativeEta();
  if(trackletsmin>=0 && trackletsmax>=0 && trackletsmax>trackletsmin) dplusTask->SetCutOnNtracklets(kTRUE,trackletsmin,trackletsmax);
  mgr->AddTask(dplusTask);

  // Create containers for input/output

  TString inname = "cinputDplus";
  TString outname = "coutputDplus";
  TString cutsname = "coutputDplusCuts";
  TString normname = "coutputDplusNorm";
  TString treename = "coutputDplus2";
  inname += finDirname.Data();
  outname += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  treename += finDirname.Data();
  TString centr=Form("%.0f%.0f",analysiscuts->GetMinCentrality(),analysiscuts->GetMaxCentrality());
  inname += centr;
  outname += centr;
  cutsname += centr;
  normname += centr;
  treename += centr;


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
  AliAnalysisDataContainer *coutputDplus2 = 0x0;
  if(storeTree){
    coutputDplus2 = mgr->CreateContainer(treename,TTree::Class(),
								   AliAnalysisManager::kOutputContainer,
								   outputfile.Data());

    coutputDplus2->SetSpecialOutput();
  }
  mgr->ConnectInput(dplusTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(dplusTask,1,coutputDplus);

  mgr->ConnectOutput(dplusTask,2,coutputDplusCuts);

  mgr->ConnectOutput(dplusTask,3,coutputDplusNorm);
  if(storeTree){
    mgr->ConnectOutput(dplusTask,4,coutputDplus2);
  }

  return dplusTask;
}
