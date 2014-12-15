Bool_t AliCFMuonSingleTask1(Int_t runmin = 17, Int_t runmax = 17)
{
    Bool_t isMC=kTRUE;
    
    TBenchmark benchmark;
    benchmark.Start("AliMuonSingleTask1");
    AliLog::SetGlobalDebugLevel(0);
    Load() ; // load the required libraries
    TChain * analysisChain ;
    analysisChain = new TChain("esdTree");

    Char_t RunFile[256];
    for(Int_t i=runmin; i<=runmax; i++){
	sprintf(RunFile,"/home/lopez/alice/data/TEST/run%d-10/AliESDs.root",i);
	analysisChain->Add(RunFile);
    }

    enum             { kEta,  kY, kPhi, kPt, kP3, kHit, kChi2Fit, kTrM, kChi2TrM,  kContrN,  kVt,  kVz, kTrig, kDCA, kZcoor, kRabs, kCharge, kTheta,kNVars };
    Int_t  nBins[] = {   5 , 5  ,  45 , 60 ,150 ,  20 ,      20 ,  4  ,      20 ,     202 , 100 ,  100 ,  10 , 500 ,   1000,    7 ,     3 , 100      };
    Double_t min[] = {  -4.,-4. ,-180.,  0.,  0.,   0.,       0., -0.5,       0.,     -2.5,   0., -100.,   0.,   0., -3000.,  171.,  -1.5 , 2.95      };
    Double_t max[] = {-2.5.,-2.5, 180., 30.,150.,  20.,      20.,  3.5,      10.,    199.5, 200.,  100.,  10., 500.,  1000.,  178.,   1.5 , 3.15     };
    
    Double_t *binLimits = 0;
    Int_t nSteps=1; if (isMC) nSteps=2;
    AliCFContainer* contCF = new AliCFContainer("container", "", nSteps, kNVars, nBins);
    for (Int_t var=kNVars; var--;) {
	binLimits = new Double_t[nBins[var]+1];
	for (Int_t i=0; i<=nBins[var]; i++) binLimits[i]=min[var]+i*(max[var]-min[var])/nBins[var];
	contCF->SetBinLimits(var, binLimits);
	delete [] binLimits; binLimits=0;
    }
    
    TList *qaList = new TList();
    TObjArray *genList = new TObjArray(0);
    AliCFTrackKineCuts *genKineCuts = new AliCFTrackKineCuts("genKineCuts", "GenKineCuts");
    genKineCuts->SetPtRange(min[kPt], max[kPt]);
    genKineCuts->SetRapidityRange(min[kY], max[kY]);
    genKineCuts->SetQAOn(qaList);
    genList->AddLast(genKineCuts);

    TObjArray *recList = new TObjArray(0);
    AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts", "RecKineCuts");
    recKineCuts->SetPtRange(min[kPt], max[kPt]);
    recKineCuts->SetRapidityRange(min[kY], max[kY]);
    recKineCuts->SetQAOn(qaList);
    recList->AddLast(recKineCuts);
//__

    AliCFManager* managerCF = new AliCFManager() ;
    managerCF->SetParticleContainer(contCF);
    managerCF->SetParticleCutsList(AliCFManager::kPartGenCuts, genList);
    managerCF->SetParticleCutsList(AliCFManager::kPartAccCuts, recList);

    AliCFMuonSingleTask1 *taskMuonCF = new AliCFMuonSingleTask1("AliMuonSingleTask1");
    taskMuonCF->SetCFManager(managerCF);
    taskMuonCF->SetQAList(qaList);
    taskMuonCF->SetUseMC(isMC);

    AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
    mgr->SetAnalysisType(AliAnalysisManager::kLocalAnalysis);
    
    AliMCEventHandler*  mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);

//  AliInputEventHandler* dataHandler ;
//  dataHandler = new AliESDInputHandler();
//  mgr->SetInputEventHandler(dataHandler);

    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdHandler);


    AliAnalysisDataContainer *cinput0  = mgr->CreateContainer("cchain0",TChain::Class(),AliAnalysisManager::kInputContainer);
    cinput0->SetData(analysisChain);

    mgr->AddTask(taskMuonCF);
    mgr->ConnectInput(taskMuonCF, 0, mgr->GetCommonInputContainer());  

    Char_t fileName[256];
    sprintf(fileName,"muonCF_%d_%d.root",runmin,runmax);
    printf("Analysis output in %s \n",fileName);

    mgr->ConnectOutput(taskMuonCF,1,mgr->CreateContainer("chist",TH1I::Class(),AliAnalysisManager::kOutputContainer,fileName));
    mgr->ConnectOutput(taskMuonCF,2,mgr->CreateContainer("ccont",AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,fileName));
 
    printf("READY TO RUN\n");
    //RUN !!!
    if (mgr->InitAnalysis()) {
	mgr->PrintStatus();
	mgr->StartAnalysis("local",analysisChain);
  }
    benchmark.Stop("AliMuonSingleTask1");
    benchmark.Show("AliMuonSingleTask1");

    return kTRUE ;
}

void Load() {

  //load the required aliroot libraries
    gSystem->Load("libANALYSIS") ;
    gSystem->Load("libANALYSISalice") ;
    gSystem->Load("$ALICE_ROOT/lib/tgt_linux/libCORRFW") ;

  //compile online the task class
    gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_ROOT/MUON -I$ALICE_ROOT/STEER -I$ROOTSYS/include");
    gROOT->LoadMacro("./AliCFMuonSingleTask1.cxx+");
}
