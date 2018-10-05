AliAnalysisTaskSEHFTreeCreator *AddTaskHFTreeCreator(Bool_t readMC=kTRUE,
                                                     Int_t system=1/*0=pp,1=PbPb*/,
                                                     TString finDirname="HFTreeCreator",
                                                     TString cutsfile="",
                                                     /*Float_t minC=0., Float_t maxC=100.,*/
                                                     Int_t AODProtection = 1,
                                                     Bool_t writeOnlySignalTree=kFALSE,
                                                     Int_t fillTreeD0=1,
                                                     Int_t fillTreeDs=1,
                                                     Int_t fillTreeDplus=1,
                                                     Int_t pidOptD0=AliHFCutOptTreeHandler::kNoPID,
                                                     Int_t pidOptDs=AliHFCutOptTreeHandler::kNsigmaPIDchar,
                                                     Int_t pidOptDplus=AliHFCutOptTreeHandler::kNsigmaPIDchar)
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
    TString cutsname = "coutputCuts";
    TString normname = "coutputNorm";
    TString treename = "coutputTree";

    inname += finDirname.Data();
    histoname += finDirname.Data();
    cutsname += finDirname.Data();
    normname += finDirname.Data();
    treename += finDirname.Data();
   
    
    
    AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_TreeCreator";
    
    AliAnalysisDataContainer *coutputEntries = mgr->CreateContainer(histoname,TH1F::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputCuts    = mgr->CreateContainer(cutsname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputNorm    = mgr->CreateContainer(normname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputTree    = mgr->CreateContainer(treename,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

    
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,coutputEntries);
    mgr->ConnectOutput(task,2,coutputCuts);
    mgr->ConnectOutput(task,3,coutputNorm);
    mgr->ConnectOutput(task,4,coutputTree);
       
    return task;
    
    
    
}

