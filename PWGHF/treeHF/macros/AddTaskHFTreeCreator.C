AliAnalysisTaskSEHFTreeCreator *AddTaskHFTreeCreator(Bool_t readMC=kFALSE,
                                                     Int_t system=0/*0=pp,1=PbPb*/,
                                                     TString finDirname="HFTreeCreator",
                                                     TString cutsfile="",
                                                     Int_t AODProtection = 1,
                                                     Bool_t writeOnlySignalTree=kFALSE,
                                                     Bool_t fillMGgenTrees=kFALSE,
                                                     Int_t fillTreeD0=0,
                                                     Int_t fillTreeDs=0,
                                                     Int_t fillTreeDplus=0,
                                                     Int_t fillTreeDstar=0,
                                                     Int_t fillTreeLctopKpi=0,
                                                     Int_t fillTreeLc2V0bachelor=0,
                                                     Int_t fillTreeBplus=0,
                                                     Int_t fillTreeBs=0,
                                                     Int_t fillTreeLb=0,
                                                     Int_t pidOpt=AliHFTreeHandler::kNsigmaPID,
                                                     Int_t singletrackvarsopt=AliHFTreeHandler::kRedSingleTrackVars,
						     Bool_t fillParticleTree = kFALSE,
                                                     Bool_t fillTrackletTree = kFALSE,
                                                     Int_t fillNJetTrees = 0,
                                                     Bool_t fillJetConstituentTrees = kFALSE,
                                                     Bool_t isITSUpgradeProd = kFALSE,
						     Bool_t fillInclusiveJetTree = kFALSE)
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
    if(!looseCutsD0toKpi && fillTreeD0)         ::Fatal("AddTaskHFTreeCreator", "looseCutsD0toKpi : check your cut file");
    AliRDHFCutsDstoKKpi*     looseCutsDstoKKpi       =(AliRDHFCutsDstoKKpi*)filecuts->Get("DstoKKpiFilteringCuts");
    if(!looseCutsDstoKKpi && fillTreeDs)        ::Fatal("AddTaskHFTreeCreator", "looseCutsDstoKKpi : check your cut file");
    AliRDHFCutsDplustoKpipi* looseCutsDplustoKpipi   =(AliRDHFCutsDplustoKpipi*)filecuts->Get("DplustoKpipiFilteringCuts");
    if(!looseCutsDplustoKpipi && fillTreeDplus)    ::Fatal("AddTaskHFTreeCreator", "looseCutsDplustoKpipi : check your cut file");
    AliRDHFCutsLctopKpi* looseCutsLctopKpi           =(AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiFilteringCuts");
    if(!looseCutsLctopKpi && fillTreeLctopKpi)        ::Fatal("AddTaskHFTreeCreator", "looseCutsLctopKpi : check your cut file");
    AliRDHFCutsD0toKpi* looseCutsBplustoD0pi   =(AliRDHFCutsD0toKpi*)filecuts->Get("BplustoD0piFilteringCuts");
    if(!looseCutsBplustoD0pi && fillTreeBplus)    ::Fatal("AddTaskHFTreeCreator", "looseCutsBplustoD0pi : check your cut file");
    AliRDHFCutsDstoKKpi* looseCutsBstoDspi   =(AliRDHFCutsDstoKKpi*)filecuts->Get("BstoDspiFilteringCuts");
    if(!looseCutsBstoDspi && fillTreeBs)    ::Fatal("AddTaskHFTreeCreator", "looseCutsBstoDspi : check your cut file");
    AliRDHFCutsDStartoKpipi* looseCutsDstartoKpipi   =(AliRDHFCutsDStartoKpipi*)filecuts->Get("DstartoKpipiFilteringCuts");
    if(!looseCutsDstartoKpipi && fillTreeDstar)    ::Fatal("AddTaskHFTreeCreator", "looseCutsDstartoKpipi : check your cut file");
    AliRDHFCutsLctoV0* looseCutsLc2V0bachelor        =(AliRDHFCutsLctoV0*)filecuts->Get("Lc2V0bachelorFilteringCuts");
    if(!looseCutsLc2V0bachelor && fillTreeLc2V0bachelor)    ::Fatal("AddTaskHFTreeCreator", "looseCutsLc2V0bachelor : check your cut file");
    AliRDHFCutsLctopKpi* looseCutsLbtoLcpi   =(AliRDHFCutsLctopKpi*)filecuts->Get("LbtoLcpiFilteringCuts");
    if(!looseCutsLbtoLcpi && fillTreeLb)    ::Fatal("AddTaskHFTreeCreator", "looseCutsLbtoLcpi : check your cut file");
    AliRDHFCutsD0toKpi*      analysisCutsD0toKpi     =(AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiAnalysisCuts");
    if(!analysisCutsD0toKpi && fillTreeD0)      ::Fatal("AddTaskHFTreeCreator", "analysisCutsD0toKpi : check your cut file");
    AliRDHFCutsDstoKKpi*     analysisCutsDstoKKpi    =(AliRDHFCutsDstoKKpi*)filecuts->Get("DstoKKpiAnalysisCuts");
    if(!analysisCutsDstoKKpi && fillTreeDs)     ::Fatal("AddTaskHFTreeCreator", "analysisCutsDstoKKpi : check your cut file");
    AliRDHFCutsDplustoKpipi* analysisCutsDplustoKpipi=(AliRDHFCutsDplustoKpipi*)filecuts->Get("DplustoKpipiAnalysisCuts");
    if(!analysisCutsDplustoKpipi && fillTreeDplus) ::Fatal("AddTaskHFTreeCreator", "analysisCutsDplustoKpipi : check your cut file");
    AliRDHFCutsLctopKpi*     analysisCutsLctopKpi    =(AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiAnalysisCuts");
    if(!analysisCutsLctopKpi && fillTreeLctopKpi) ::Fatal("AddTaskHFTreeCreator", "analysisCutsLctopKpi : check your cut file");
    AliRDHFCutsD0toKpi* analysisCutsBplustoD0pi=(AliRDHFCutsD0toKpi*)filecuts->Get("BplustoD0piAnalysisCuts");
    if(!analysisCutsBplustoD0pi && fillTreeBplus) ::Fatal("AddTaskHFTreeCreator", "analysisCutsBplustoD0pi : check your cut file");
    AliRDHFCutsDstoKKpi* analysisCutsBstoDspi   =(AliRDHFCutsDstoKKpi*)filecuts->Get("BstoDspiAnalysisCuts");
    if(!analysisCutsBstoDspi && fillTreeBs)    ::Fatal("AddTaskHFTreeCreator", "analysisCutsBstoDspi : check your cut file");
    AliRDHFCutsDStartoKpipi* analysisCutsDstartoKpipi=(AliRDHFCutsDStartoKpipi*)filecuts->Get("DstartoKpipiAnalysisCuts");
    if(!analysisCutsDstartoKpipi && fillTreeDstar) ::Fatal("AddTaskHFTreeCreator", "analysisCutsDstartoKpipi : check your cut file");
    AliRDHFCutsLctoV0* analysisCutsLc2V0bachelor=(AliRDHFCutsLctoV0*)filecuts->Get("Lc2V0bachelorAnalysisCuts");
    if(!analysisCutsLc2V0bachelor && fillTreeLc2V0bachelor) ::Fatal("AddTaskHFTreeCreator", "analysisCutsLc2V0bachelor : check your cut file");
    AliRDHFCutsLctopKpi* analysisCutsLbtoLcpi   =(AliRDHFCutsLctopKpi*)filecuts->Get("LbtoLcpiAnalysisCuts");
    if(!analysisCutsLbtoLcpi && fillTreeLb)    ::Fatal("AddTaskHFTreeCreator", "analysisCutsLbtoLcpi : check your cut file");

    TList *cutsList=new TList();
    cutsList->SetOwner(kTRUE);
    cutsList->SetName("cut_objects");
    cutsList->Add(looseCutsD0toKpi);
    cutsList->Add(looseCutsDstoKKpi);
    cutsList->Add(looseCutsDplustoKpipi);
    cutsList->Add(looseCutsLctopKpi);
    cutsList->Add(looseCutsBplustoD0pi);
    cutsList->Add(looseCutsBstoDspi);
    cutsList->Add(looseCutsDstartoKpipi);
    cutsList->Add(looseCutsLc2V0bachelor);
    cutsList->Add(looseCutsLbtoLcpi);
    cutsList->Add(analysisCutsD0toKpi);
    cutsList->Add(analysisCutsDstoKKpi);
    cutsList->Add(analysisCutsDplustoKpipi);
    cutsList->Add(analysisCutsLctopKpi);
    cutsList->Add(analysisCutsBplustoD0pi);
    cutsList->Add(analysisCutsBstoDspi);
    cutsList->Add(analysisCutsDstartoKpipi);
    cutsList->Add(analysisCutsLc2V0bachelor);
    cutsList->Add(analysisCutsLbtoLcpi);

    AliAnalysisTaskSEHFTreeCreator *task = new AliAnalysisTaskSEHFTreeCreator("TreeCreatorTask",cutsList, fillNJetTrees, fillJetConstituentTrees);

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
    task->SetFillLctopKpiTree(fillTreeLctopKpi);
    task->SetFillBplusTree(fillTreeBplus);
    task->SetFillBsTree(fillTreeBs);
    task->SetFillDstarTree(fillTreeDstar);
    task->SetFillLc2V0bachelorTree(fillTreeLc2V0bachelor);
    task->SetFillLbTree(fillTreeLb);
    task->SetFillInclusiveJetTree(fillInclusiveJetTree);
    task->SetFillParticleTree(fillParticleTree);
    task->SetFillTrackletTree(fillTrackletTree);
    task->SetPIDoptD0Tree(pidOpt);
    task->SetPIDoptDsTree(pidOpt);
    task->SetPIDoptDplusTree(pidOpt);
    task->SetPIDoptLctopKpiTree(pidOpt);
    task->SetPIDoptBplusTree(pidOpt);
    task->SetPIDoptBsTree(pidOpt);
    task->SetPIDoptDstarTree(pidOpt);
    task->SetPIDoptLc2V0bachelorTree(pidOpt);
    task->SetPIDoptLbTree(pidOpt);
    task->SetTreeSingleTrackVarsOpt(singletrackvarsopt);
    if(fillTreeBs || fillTreeLb || fillTreeBplus || isITSUpgradeProd){
      task->SetITSUpgradeProduction(kTRUE);
      task->SetITSUpgradePreSelect(kTRUE);
      task->SetStoreOnlyHIJINGBackground(kTRUE);
      task->SetGoodTrackFilterBit(4);
      task->SetGoodTrackEtaRange(0.8);
      task->SetGoodTrackMinPt(0.3);
    }
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
    TString treeLctopKpiname = "coutputTreeLctopKpi";
    TString treeBplusname = "coutputTreeBplus";
    TString treeBsname = "coutputTreeBs";
    TString treeDstarname = "coutputTreeDstar";
    TString treeLc2V0bachelorname = "coutputTreeLc2V0bachelor";
    TString treeLbname = "coutputTreeLb";
    TString treeGenD0name = "coutputTreeGenD0";
    TString treeGenDplusname = "coutputTreeGenDplus";
    TString treeGenDsname = "coutputTreeGenDs";
    TString treeGenLctopKpiname = "coutputTreeGenLctopKpi";
    TString treeGenBplusname = "coutputTreeGenBplus";
    TString treeGenBsname = "coutputTreeGenBs";
    TString treeGenDstarname = "coutputTreeGenDstar";
    TString treeGenLc2V0bachelorname = "coutputTreeGenLc2V0bachelor";
    TString treeGenLbname = "coutputTreeGenLb";
    TString treeParticleName = "coutputTreeParticle";
    TString treeTrackletName = "coutputTreeTracklet";
    TString treeGenParticleName = "coutputTreeGenParticle";
    TString treeJetName = "coutputTreeJet%d";
    TString treeJetConstituentName = "coutputTreeJetConstituent%d";
 
    inname += finDirname.Data();
    histoname += finDirname.Data();
    countername += finDirname.Data();
    cutsname += finDirname.Data();
    normname += finDirname.Data();
    treeevcharname += finDirname.Data();
    treeD0name += finDirname.Data();
    treeDplusname += finDirname.Data();
    treeDsname += finDirname.Data();
    treeLctopKpiname += finDirname.Data();
    treeBplusname += finDirname.Data();
    treeBsname += finDirname.Data();
    treeDstarname += finDirname.Data();
    treeLc2V0bachelorname += finDirname.Data();
    treeLbname += finDirname.Data();
    treeGenD0name += finDirname.Data();
    treeGenDplusname += finDirname.Data();
    treeGenDsname += finDirname.Data();
    treeGenLctopKpiname += finDirname.Data();
    treeGenBplusname += finDirname.Data();
    treeGenBsname += finDirname.Data();
    treeGenDstarname += finDirname.Data();
    treeGenLc2V0bachelorname += finDirname.Data();
    treeGenLbname += finDirname.Data();
    treeParticleName += finDirname.Data();
    treeTrackletName += finDirname.Data();
    treeGenParticleName += finDirname.Data();
    treeJetName += finDirname.Data();
    treeJetConstituentName += finDirname.Data();

    AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_TreeCreator";
    if(fillTreeLb || fillTreeBs || fillTreeBplus || isITSUpgradeProd){
      //Needed to run ITS2 production together with ITS2+ITS3 improver
      outputfile += finDirname;
    }

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
    AliAnalysisDataContainer *coutputTreeLctopKpi = 0x0;
    AliAnalysisDataContainer *coutputTreeGenLctopKpi = 0x0;
    AliAnalysisDataContainer *coutputTreeBplus = 0x0;
    AliAnalysisDataContainer *coutputTreeGenBplus = 0x0;
    AliAnalysisDataContainer *coutputTreeDstar = 0x0;
    AliAnalysisDataContainer *coutputTreeGenDstar = 0x0;
    AliAnalysisDataContainer *coutputTreeLc2V0bachelor = 0x0;
    AliAnalysisDataContainer *coutputTreeGenLc2V0bachelor = 0x0;
    AliAnalysisDataContainer *coutputTreeBs = 0x0;
    AliAnalysisDataContainer *coutputTreeGenBs = 0x0;
    AliAnalysisDataContainer *coutputTreeLb = 0x0;
    AliAnalysisDataContainer *coutputTreeGenLb = 0x0;
    AliAnalysisDataContainer *coutputTreeParticle = 0x0;
    AliAnalysisDataContainer *coutputTreeTracklet = 0x0;
    AliAnalysisDataContainer *coutputTreeGenParticle = 0x0;
    std::vector<AliAnalysisDataContainer*> coutputTreeJet;
    std::vector<AliAnalysisDataContainer*> coutputTreeJetConstituent;

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

    if(fillTreeLctopKpi) {
      coutputTreeLctopKpi = mgr->CreateContainer(treeLctopKpiname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeLctopKpi->SetSpecialOutput();
      if(readMC && fillMGgenTrees) {
        coutputTreeGenLctopKpi = mgr->CreateContainer(treeGenLctopKpiname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeGenLctopKpi->SetSpecialOutput();
      }
    }
    
    if(fillTreeBplus) {
        coutputTreeBplus = mgr->CreateContainer(treeBplusname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeBplus->SetSpecialOutput();
        if(readMC && fillMGgenTrees) {
            coutputTreeGenBplus = mgr->CreateContainer(treeGenBplusname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
            coutputTreeGenBplus->SetSpecialOutput();
        }
    }

    if(fillTreeDstar) {
        coutputTreeDstar = mgr->CreateContainer(treeDstarname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeDstar->SetSpecialOutput();
        if(readMC && fillMGgenTrees) {
            coutputTreeGenDstar = mgr->CreateContainer(treeGenDstarname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
            coutputTreeGenDstar->SetSpecialOutput();
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
  
    if(fillTreeBs) {
      coutputTreeBs = mgr->CreateContainer(treeBsname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeBs->SetSpecialOutput();
      if(readMC && fillMGgenTrees) {
        coutputTreeGenBs = mgr->CreateContainer(treeGenBsname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeGenBs->SetSpecialOutput();
      }
    }
  
    if(fillTreeLb) {
      coutputTreeLb = mgr->CreateContainer(treeLbname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeLb->SetSpecialOutput();
      if(readMC && fillMGgenTrees) {
        coutputTreeGenLb = mgr->CreateContainer(treeGenLbname,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeGenLb->SetSpecialOutput();
      }
    }
  
    if(fillParticleTree) {
      coutputTreeParticle = mgr->CreateContainer(treeParticleName,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeParticle->SetSpecialOutput();
      if(readMC && fillMGgenTrees) {
        coutputTreeGenParticle = mgr->CreateContainer(treeGenParticleName,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
        coutputTreeGenParticle->SetSpecialOutput();
      }
    }
  
    if(fillTrackletTree) {
      coutputTreeTracklet = mgr->CreateContainer(treeTrackletName,TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
      coutputTreeTracklet->SetSpecialOutput();
    }
    for (int i=0; i<fillNJetTrees; i++) {
      coutputTreeJet.push_back(mgr->CreateContainer(Form(treeJetName, i),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()));
      coutputTreeJet.at(i)->SetSpecialOutput();
    }
  
    if (fillJetConstituentTrees) {
      for (int i=0; i<fillNJetTrees; i++) {
        coutputTreeJetConstituent.push_back(mgr->CreateContainer(Form(treeJetConstituentName, i),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()));
        coutputTreeJetConstituent.at(i)->SetSpecialOutput();
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
    if(fillTreeLctopKpi) {
      mgr->ConnectOutput(task,12,coutputTreeLctopKpi);
      if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,13,coutputTreeGenLctopKpi);
    }
    if(fillTreeBplus) {
        mgr->ConnectOutput(task,14,coutputTreeBplus);
        if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,15,coutputTreeGenBplus);
    }
    if(fillTreeDstar) {
        mgr->ConnectOutput(task,16,coutputTreeDstar);
        if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,17,coutputTreeGenDstar);
    }
    if(fillTreeLc2V0bachelor) {
        mgr->ConnectOutput(task,18,coutputTreeLc2V0bachelor);
        if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,19,coutputTreeGenLc2V0bachelor);
    }
    if(fillTreeBs) {
      mgr->ConnectOutput(task,20,coutputTreeBs);
      if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,21,coutputTreeGenBs);
    }
    if(fillTreeLb) {
      mgr->ConnectOutput(task,22,coutputTreeLb);
      if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,23,coutputTreeGenLb);
    }
    if(fillParticleTree) {
      mgr->ConnectOutput(task,24,coutputTreeParticle);
      if(readMC && fillMGgenTrees) mgr->ConnectOutput(task,25,coutputTreeGenParticle);
    }
    if(fillTrackletTree) {
      mgr->ConnectOutput(task,26,coutputTreeTracklet);
    }
    for (int i=0; i<fillNJetTrees; i++) {
      mgr->ConnectOutput(task,27+i,coutputTreeJet.at(i));
    }
    if (fillJetConstituentTrees) {
      for (int i=0; i<fillNJetTrees; i++) {
        mgr->ConnectOutput(task,27+fillNJetTrees+i,coutputTreeJetConstituent.at(i));
      }
    }

    return task;
}
