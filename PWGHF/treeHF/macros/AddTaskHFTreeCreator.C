AliAnalysisTaskSEHFTreeCreator *AddTaskHFTreeCreator(Bool_t readMC=kTRUE,
                                                     Int_t system=1/*0=pp,1=PbPb*/,
                                                     TString finDirname="HFTreeCreator",
                                                     TString cutsfile="",
                                                     Int_t AODProtection = 1,
                                                     Bool_t writeOnlySignalTree=kFALSE,
                                                     Bool_t fillMGgenTrees=kFALSE,
                                                     Int_t fillTreeD0=1,
                                                     Int_t fillTreeDs=1,
                                                     Int_t fillTreeDplus=1,
                                                     Int_t pidOptD0=AliHFTreeHandler::kRawAndNsigmaPID,
                                                     Int_t pidOptDs=AliHFTreeHandler::kRawAndNsigmaPID,
                                                     Int_t pidOptDplus=AliHFTreeHandler::kRawAndNsigmaPID)
{
    //
    //


    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
        return NULL;
    }

    //getting the cuts
    TFile* filecuts;
    if( cutsfile.EqualTo("") ) {
        ::Fatal("AddTaskHFTreeCreator", "Input file not provided");
    } else {
        filecuts=TFile::Open(cutsfile.Data());
        if(!filecuts || (filecuts && !filecuts->IsOpen())){
            ::Fatal("AddTaskHFTreeCreator", "Input file not found : check your cut object");
        }
    }

    AliRDHFCutsD0toKpi*      looseCutsD0toKpi        =(AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiFilteringCuts");
    if(!looseCutsD0toKpi)         ::Fatal("AddTaskHFTreeCreator", "looseCutsD0toKpi : check your cut file");
    AliRDHFCutsDstoKKpi*     looseCutsDstoKKpi       =(AliRDHFCutsDstoKKpi*)filecuts->Get("DstoKKpiFilteringCuts");
    if(!looseCutsDstoKKpi)        ::Fatal("AddTaskHFTreeCreator", "looseCutsDstoKKpi : check your cut file");
    AliRDHFCutsDplustoKpipi* looseCutsDplustoKpipi   =(AliRDHFCutsDplustoKpipi*)filecuts->Get("DplustoKpipiFilteringCuts");
    if(!looseCutsDplustoKpipi)    ::Fatal("AddTaskHFTreeCreator", "looseCutsDplustoKpipi : check your cut file");
    AliRDHFCutsD0toKpi*      analysisCutsD0toKpi     =(AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiAnalysisCuts");
    if(!analysisCutsD0toKpi)      ::Fatal("AddTaskHFTreeCreator", "analysisCutsD0toKpi : check your cut file");
    AliRDHFCutsDstoKKpi*     analysisCutsDstoKKpi    =(AliRDHFCutsDstoKKpi*)filecuts->Get("DstoKKpiAnalysisCuts");
    if(!analysisCutsDstoKKpi)     ::Fatal("AddTaskHFTreeCreator", "analysisCutsDstoKKpi : check your cut file");
    AliRDHFCutsDplustoKpipi* analysisCutsDplustoKpipi=(AliRDHFCutsDplustoKpipi*)filecuts->Get("DplustoKpipiAnalysisCuts");
    if(!analysisCutsDplustoKpipi) ::Fatal("AddTaskHFTreeCreator", "analysisCutsDplustoKpipi : check your cut file");

    TList *cutsList=new TList();
    cutsList->SetOwner(kTRUE);
    cutsList->SetName("cut_objects");
    cutsList->Add(looseCutsD0toKpi);
    cutsList->Add(looseCutsDstoKKpi);
    cutsList->Add(looseCutsDplustoKpipi);
    cutsList->Add(analysisCutsD0toKpi);
    cutsList->Add(analysisCutsDstoKKpi);
    cutsList->Add(analysisCutsDplustoKpipi);

    AliAnalysisTaskSEHFTreeCreator *task = new AliAnalysisTaskSEHFTreeCreator("TreeCreatorTask",cutsList);

    task->SetReadMC(readMC);
    if(readMC) {
      task->SetFillMCGenTrees(fillMGgenTrees);
    }
    task->SetSystem(system);
    task->SetAODMismatchProtection(AODProtection);
    task->SetWriteOnlySignalTree(writeOnlySignalTree);
    task->SetFillD0Tree(fillTreeD0);
    task->SetFillDsTree(fillTreeDs);
    task->SetFillDplusTree(fillTreeDplus);
    task->SetPIDoptD0Tree(pidOptD0);
    task->SetPIDoptDsTree(pidOptDs);
    task->SetPIDoptDplusTree(pidOptDplus);

    //task->SetDebugLevel(4);

    mgr->AddTask(task);

    // Create containers for input/output

    TString inname = "cinput";
    TString histoname = "coutputEntries";
    TString countername = "coutputNormCounterHisto";
    TString cutsname = "coutputCuts";
    TString normname = "coutputNorm";
    TString treeevcharname = "coutputTreeEvChar";
    TString treeD0name = "coutputTreeD0";
    TString treeDplusname = "coutputTreeDplus";
    TString treeDsname = "coutputTreeDs";
    TString treeGenD0name = "coutputTreeGenD0";
    TString treeGenDplusname = "coutputTreeGenDplus";
    TString treeGenDsname = "coutputTreeGenDs";

    inname += finDirname.Data();
    histoname += finDirname.Data();
    countername += finDirname.Data();
    cutsname += finDirname.Data();
    normname += finDirname.Data();
    treeevcharname += finDirname.Data();
    treeD0name += finDirname.Data();
    treeDplusname += finDirname.Data();
    treeDsname += finDirname.Data();
    treeGenD0name += finDirname.Data();
    treeGenDplusname += finDirname.Data();
    treeGenDsname += finDirname.Data();

    AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_TreeCreator";

    AliAnalysisDataContainer *coutputEntries = mgr->CreateContainer(histoname,TH1F::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputCounter = mgr->CreateContainer(countername,TH2F::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputCuts    = mgr->CreateContainer(cutsname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputNorm    = mgr->CreateContainer(normname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputTreeEvChar    = mgr->CreateContainer(treeevcharname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  
    AliAnalysisDataContainer *coutputTreeD0 = 0x0;
    AliAnalysisDataContainer *coutputTreeGenD0 = 0x0;
    AliAnalysisDataContainer *coutputTreeDplus = 0x0;
    AliAnalysisDataContainer *coutputTreeGenDplus = 0x0;
    AliAnalysisDataContainer *coutputTreeDs = 0x0;
    AliAnalysisDataContainer *coutputTreeGenDs = 0x0;

    if(fillTreeD0) {
      coutputTreeD0 = mgr->CreateContainer(treeD0name,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeD0->SetSpecialOutput();
      if(readMC && fillMGgenTrees) {
        coutputTreeGenD0 = mgr->CreateContainer(treeGenD0name,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeGenD0->SetSpecialOutput();
      }
    }
  
    if(fillTreeDplus) {
      coutputTreeDplus = mgr->CreateContainer(treeDplusname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeDplus->SetSpecialOutput();
      if(readMC && fillMGgenTrees) {
        coutputTreeGenDplus = mgr->CreateContainer(treeGenDplusname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeGenDplus->SetSpecialOutput();
      }
    }
  
    if(fillTreeDs) {
      coutputTreeDs = mgr->CreateContainer(treeDsname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeDs->SetSpecialOutput();
      if(readMC && fillMGgenTrees) {
        coutputTreeGenDs = mgr->CreateContainer(treeGenDsname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeGenDs->SetSpecialOutput();
      }
    }
      
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,coutputEntries);
    mgr->ConnectOutput(task,2,coutputCounter);
    mgr->ConnectOutput(task,3,coutputCuts);
    mgr->ConnectOutput(task,4,coutputNorm);
    mgr->ConnectOutput(task,5,coutputTreeEvChar);
    if(fillTreeD0) {
      mgr->ConnectOutput(task,6,coutputTreeD0);
      if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,7,coutputTreeGenD0);
    }
    if(fillTreeDs) {
      mgr->ConnectOutput(task,8,coutputTreeDs);
      if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,9,coutputTreeGenDs);
    }
    if(fillTreeDplus) {
      mgr->ConnectOutput(task,10,coutputTreeDplus);
      if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,11,coutputTreeGenDplus);
    }
  
    return task;
}
