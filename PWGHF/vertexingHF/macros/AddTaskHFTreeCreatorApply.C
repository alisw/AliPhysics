AliAnalysisTaskSEHFTreeCreatorApply *AddTaskHFTreeCreatorApply(Bool_t readMC=kFALSE,
                                                               Int_t system=0 /*0=pp,1=PbPb*/,
                                                               TString finDirname="HFTreeCreator",
                                                               TString cutsfile="",
                                                               TString confFileML = "",
                                                               Int_t AODProtection = 0,
                                                               Bool_t writeOnlySignalTree=kFALSE,
                                                               Bool_t fillMGgenTrees=kFALSE,
                                                               Int_t fillTreeDs=0,
                                                               Int_t fillTreeLc2V0bachelor=0,
                                                               Int_t fillTreeDstar=0,
                                                               Int_t fillTreeParticle=0, //0=None, 1=pKpi, 2=p, 3=K, 4=pi
                                                               Int_t pidOpt=AliHFTreeHandlerApply::kNsigmaPID,
                                                               Int_t singletrackvarsopt=AliHFTreeHandlerApply::kRedSingleTrackVars,
                                                               Bool_t useDefaultDirName=kTRUE,
                                                               Bool_t fillMixEvBkg=kFALSE)
{
  //
  //
  
  if(fillTreeDs && fillTreeLc2V0bachelor && fillTreeDstar && confFileML != ""){
    ::Error("AddTaskHFTreeCreatorApply", "Not possible to enable more than one hadron when applying.");
    return NULL;
  }
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFTreeCreatorApply", "No analysis manager to connect to.");
    return NULL;
  }
  
  //getting the cuts
  TFile* filecuts;
  if( cutsfile.EqualTo("") ) {
    ::Fatal("AddTaskHFTreeCreatorApply", "Input file not provided");
  } else {
    filecuts=TFile::Open(cutsfile.Data());
    if(!filecuts || (filecuts && !filecuts->IsOpen())){
      ::Fatal("AddTaskHFTreeCreatorApply", "Input file not found : check your cut object");
    }
  }
  
  AliRDHFCutsDstoKKpi*     looseCutsDstoKKpi       =(AliRDHFCutsDstoKKpi*)filecuts->Get("DstoKKpiFilteringCuts");
  if(!looseCutsDstoKKpi && fillTreeDs)        ::Fatal("AddTaskHFTreeCreatorApply", "looseCutsDstoKKpi : check your cut file");
  AliRDHFCutsLctoV0* looseCutsLc2V0bachelor        =(AliRDHFCutsLctoV0*)filecuts->Get("Lc2V0bachelorFilteringCuts");
  if(!looseCutsLc2V0bachelor && (fillTreeLc2V0bachelor || fillMixEvBkg))    ::Fatal("AddTaskHFTreeCreatorApply", "looseCutsLc2V0bachelor : check your cut file");
  AliRDHFCutsDStartoKpipi* looseCutsDstartoKpipi   =(AliRDHFCutsDStartoKpipi*)filecuts->Get("DstartoKpipiFilteringCuts");
  if(!looseCutsDstartoKpipi && fillTreeDstar)    ::Fatal("AddTaskHFTreeCreator", "looseCutsDstartoKpipi : check your cut file");
  AliRDHFCutsDstoKKpi*     analysisCutsDstoKKpi    =(AliRDHFCutsDstoKKpi*)filecuts->Get("DstoKKpiAnalysisCuts");
  if(!analysisCutsDstoKKpi && fillTreeDs)     ::Fatal("AddTaskHFTreeCreatorApply", "analysisCutsDstoKKpi : check your cut file");
  AliRDHFCutsLctoV0* analysisCutsLc2V0bachelor=(AliRDHFCutsLctoV0*)filecuts->Get("Lc2V0bachelorAnalysisCuts");
  if(!analysisCutsLc2V0bachelor && (fillTreeLc2V0bachelor || fillMixEvBkg)) ::Fatal("AddTaskHFTreeCreatorApply", "analysisCutsLc2V0bachelor : check your cut file");
  AliRDHFCutsDStartoKpipi* analysisCutsDstartoKpipi=(AliRDHFCutsDStartoKpipi*)filecuts->Get("DstartoKpipiAnalysisCuts");
  if(!analysisCutsDstartoKpipi && fillTreeDstar) ::Fatal("AddTaskHFTreeCreator", "analysisCutsDstartoKpipi : check your cut file");
  
  TList *cutsList=new TList();
  cutsList->SetOwner(kTRUE);
  cutsList->SetName("cut_objects");
  cutsList->Add(looseCutsDstoKKpi);
  cutsList->Add(looseCutsLc2V0bachelor);
  cutsList->Add(looseCutsDstartoKpipi);
  cutsList->Add(analysisCutsDstoKKpi);
  cutsList->Add(analysisCutsLc2V0bachelor);
  cutsList->Add(analysisCutsDstartoKpipi);
  
  AliAnalysisTaskSEHFTreeCreatorApply *task = new AliAnalysisTaskSEHFTreeCreatorApply("TreeCreatorTask",cutsList);
  
  task->SetReadMC(readMC);
  if(readMC) {
    task->SetFillMCGenTrees(fillMGgenTrees);
  }
  task->SetSystem(system);
  task->SetAODMismatchProtection(AODProtection);
  task->SetWriteOnlySignalTree(writeOnlySignalTree);
  task->SetFillDsTree(fillTreeDs);
  task->SetFillLc2V0bachelorTree(fillTreeLc2V0bachelor);
  task->SetFillDstarTree(fillTreeDstar);
  task->SetFillParticleTree(fillTreeParticle);
  task->SetPIDoptDsTree(pidOpt);
  task->SetPIDoptLc2V0bachelorTree(pidOpt);
  task->SetPIDoptDstarTree(pidOpt);
  task->SetTreeSingleTrackVarsOpt(singletrackvarsopt);
  if(fillMixEvBkg){
    task->SetSaveMixedEventBkg(kTRUE);
  }

  task->SetMLConfigFile(confFileML);
  //task->SetDebugLevel(4);
  
  mgr->AddTask(task);
  
  // Create containers for input/output
  
  TString inname = "cinput";
  TString histoname = "coutputEntries";
  TString countername = "coutputNormCounterHisto";
  TString cutsname = "coutputCuts";
  TString normname = "coutputNorm";
  TString treeevcharname = "coutputTreeEvChar";
  TString treeDsname = "coutputTreeDs";
  TString treeLc2V0bachelorname = "coutputTreeLc2V0bachelor";
  TString treeLc2V0bachelornameMixEv = "coutputTreeLc2V0bachelorMixEv";
  TString treeDstarname = "coutputTreeDstar";
  TString treeParticlename = "coutputTreeParticle";
  TString treeGenDsname = "coutputTreeGenDs";
  TString treeGenLc2V0bachelorname = "coutputTreeGenLc2V0bachelor";
  TString treeGenDstarname = "coutputTreeGenDstar";
  
  inname += finDirname.Data();
  histoname += finDirname.Data();
  countername += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  treeevcharname += finDirname.Data();
  treeDsname += finDirname.Data();
  treeLc2V0bachelorname += finDirname.Data();
  treeLc2V0bachelornameMixEv += finDirname.Data();
  treeDstarname += finDirname.Data();
  treeParticlename += finDirname.Data();
  treeGenDsname += finDirname.Data();
  treeGenLc2V0bachelorname += finDirname.Data();
  treeGenDstarname += finDirname.Data();
  
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWGHF_TreeCreatorApply";
  if(!useDefaultDirName) outputfile += finDirname;
  
  AliAnalysisDataContainer *coutputEntries = mgr->CreateContainer(histoname,TH1F::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputCounter = mgr->CreateContainer(countername,TH2F::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputCuts    = mgr->CreateContainer(cutsname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputNorm    = mgr->CreateContainer(normname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputTreeEvChar    = mgr->CreateContainer(treeevcharname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  AliAnalysisDataContainer *coutputTreeDs = 0x0;
  AliAnalysisDataContainer *coutputTreeGenDs = 0x0;
  AliAnalysisDataContainer *coutputTreeLc2V0bachelor = 0x0;
  AliAnalysisDataContainer *coutputTreeLc2V0bachelorMixEv = 0x0;
  AliAnalysisDataContainer *coutputTreeGenLc2V0bachelor = 0x0;
  AliAnalysisDataContainer *coutputTreeDstar = 0x0;
  AliAnalysisDataContainer *coutputTreeGenDstar = 0x0;
  AliAnalysisDataContainer *coutputTreeParticle = 0x0;

  if(fillTreeDs) {
    coutputTreeDs = mgr->CreateContainer(treeDsname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeDs->SetSpecialOutput();
    if(readMC && fillMGgenTrees) {
      coutputTreeGenDs = mgr->CreateContainer(treeGenDsname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeGenDs->SetSpecialOutput();
    }
  }
  
  if(fillTreeLc2V0bachelor) {
    coutputTreeLc2V0bachelor = mgr->CreateContainer(treeLc2V0bachelorname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeLc2V0bachelor->SetSpecialOutput();
    if(readMC && fillMGgenTrees) {
      coutputTreeGenLc2V0bachelor = mgr->CreateContainer(treeGenLc2V0bachelorname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeGenLc2V0bachelor->SetSpecialOutput();
    }
  }

  if(fillMixEvBkg){
    coutputTreeLc2V0bachelorMixEv = mgr->CreateContainer(treeLc2V0bachelornameMixEv,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeLc2V0bachelorMixEv->SetSpecialOutput();
  }

  if(fillTreeDstar) {
    coutputTreeDstar = mgr->CreateContainer(treeDstarname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeDstar->SetSpecialOutput();
    if(readMC && fillMGgenTrees) {
      coutputTreeGenDstar = mgr->CreateContainer(treeGenDstarname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeGenDstar->SetSpecialOutput();
    }
  }

  if(fillTreeParticle){
    coutputTreeParticle = mgr->CreateContainer(treeParticlename,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    coutputTreeParticle->SetSpecialOutput();
  }

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutputEntries);
  mgr->ConnectOutput(task,2,coutputCounter);
  mgr->ConnectOutput(task,3,coutputCuts);
  mgr->ConnectOutput(task,4,coutputNorm);
  mgr->ConnectOutput(task,5,coutputTreeEvChar);
  if(fillTreeDs) {
    mgr->ConnectOutput(task,6,coutputTreeDs);
    if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,7,coutputTreeGenDs);
  }
  if(fillTreeLc2V0bachelor) {
    mgr->ConnectOutput(task,8,coutputTreeLc2V0bachelor);
    if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,9,coutputTreeGenLc2V0bachelor);
  }
  if(fillMixEvBkg) mgr->ConnectOutput(task,10,coutputTreeLc2V0bachelorMixEv);
  if(fillTreeDstar) {
    mgr->ConnectOutput(task,11,coutputTreeDstar);
    if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,12,coutputTreeGenDstar);
  }
  if(fillTreeParticle){
    mgr->ConnectOutput(task,13,coutputTreeParticle);
  }

  return task;
}
