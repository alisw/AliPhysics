AliAnalysisTaskHFJetIPQA* AddTaskHFJetIPQA(
                                           const char *ntracks            = "Tracks",
                                           const char *nclusters           = "",
                                           const char *njets              = "Jets",
                                           const char *nrho               = "",
                                           Double_t jetradius =0.4,
                                           Bool_t isMC = kFALSE,
                                           const char * type = "TPC",
                                           const char *taskname           = "AliAnalysisTaskEmcalJetBJetTaggingIP",
                                           const char *njetsMC              = "Jets",
                                           const char *nrhoMC               = "RhoMC",
                                           TString PathToWeights = 	"alien:///alice/cern.ch/user/l/lfeldkam/Weights.root",
                                           TString PathToRunwiseCorrectionParameters = "alien:///alice/cern.ch/user/l/lfeldkam/MeanSigmaImpParFactors.root",
                                           TString PathToJetProbabilityInput = "alien:///alice/cern.ch/user/l/lfeldkam/dummy_ResFct_XYSignificance_pp7TeV.root",
                                           Bool_t GenerateMeanSigmaCorrectionTable=kTRUE,
                                           Int_t nITSReq=6,
                                           const char* suffix = ""
                                           )
{
    
    Bool_t IsESD = kFALSE;
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
        return NULL;
    }
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
        return NULL;
    }
    //Determine analysis Type
    //==============================================================================
    if (mgr) {
        AliVEventHandler *evhand = mgr->GetInputEventHandler();
        if (evhand) {
            if (evhand->InheritsFrom("AliESDInputHandler"))
            {
                IsESD = kTRUE;
                Printf("%s :: Running on ESD files",taskname);
            }
            else {
                IsESD = kFALSE;
                Printf("%s :: Running on AOD files",taskname);
            }
        }
    }
    // Load and setup Monte Carlo composition correction factors from file
    //==============================================================================
    TFile* filecorrectionfactors;
    if(isMC){
        if( PathToWeights.EqualTo("") ) {
        } else {
            filecorrectionfactors=TFile::Open(PathToWeights.Data());
            if(!filecorrectionfactors ||(filecorrectionfactors&& !filecorrectionfactors->IsOpen())){
                AliFatal("%s :: Input weight file not found",taskname);
                return 0x0;
            }
        }
    }
    
    
    Printf("%s :: File %s successfully loaded, setting up correction factors.",taskname,PathToWeights.Data());
    
    TString name(taskname),combinedName;
    combinedName.Form("%s%s", name.Data(),suffix);
    AliAnalysisTaskHFJetIPQA* jetTask = new AliAnalysisTaskHFJetIPQA(combinedName);
    if(isMC && filecorrectionfactors){
        TH1F * h[20] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
        const char * nampart[20] = {"pi0","eta","etap","rho","phi","omega","k0s","lambda","pi","kaon","proton","D0","Dp","Dsp","Ds","lambdac","bplus","b0","lambdab","bsp"};
        for(int i = 0;i<20;++i){
            h[i]=(TH1F*)filecorrectionfactors->Get(nampart[i]);
        }
        Printf("%s :: Reading correction factors completed.",taskname);
        for (int i = 0 ; i<20;++i) if(h[i]==0) return 0x0;
        jetTask->SetUseMonteCarloWeighingLinus(h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8],h[9],h[10],h[11],h[12],h[13],h[14],h[15],h[16],h[17],h[18],h[19]);
        Printf("%s :: Weights written to analysis task.",taskname);
        if(filecorrectionfactors) filecorrectionfactors->Close();
    }
    
    
    // Setup input containers
    //==============================================================================
    Printf("%s :: Setting up input containers.",taskname);
    AliParticleContainer *trackCont  = jetTask->AddParticleContainer(ntracks);
    AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);
    TString strType(type);
    AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
    
    if(jetCont) {
        jetCont->SetRhoName(nrho);
        jetCont->ConnectParticleContainer(trackCont);
        jetCont->ConnectClusterContainer(clusterCont);
    }
    
    if(isMC)
    {
        AliJetContainer *jetContMC = jetTask->AddJetContainer(njetsMC,strType,jetradius);
        
        if(jetContMC) {
            jetContMC->SetRhoName(nrhoMC);
            jetContMC->SetIsParticleLevel(kTRUE);
            jetContMC->SetMaxTrackPt(1000);
        }
    }
    //==============================================================================
    TH1::AddDirectory(0);
    if( PathToJetProbabilityInput.EqualTo("") ) {
    } else {
        jetTask->SetResFunctionPID(PathToJetProbabilityInput.Data());
    }
    // Set Monte Carlo / Data status
    //==============================================================================
    jetTask->SetIsPythia(isMC);
    // Setup jet cuts
    //==============================================================================
    DefineCutsTaskpp(jetTask,-1.,100);
    if(IsESD) {
        // Setup initial ESD track cuts
        //==============================================================================
        jetTask->SetRunESD();
        AliESDtrackCuts * trackCuts = new AliESDtrackCuts();
        // trackCuts->SetMaxFractionSharedTPCClusters(0.4);
        // TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
        //trackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
        trackCuts->SetMinNClustersTPC(100);
        trackCuts->SetMaxChi2PerClusterTPC(20);
        //  trackCuts->SetRequireITSStandAlone(kTRUE);
        // trackCuts->SetRequireITSPureStandAlone(kTRUE);
        trackCuts->SetAcceptKinkDaughters(kFALSE);
        trackCuts->SetRequireTPCRefit(kTRUE);
        trackCuts->SetRequireITSRefit(kTRUE);
        trackCuts->SetMaxChi2PerClusterITS(4);
        trackCuts->SetEtaRange(-0.9, 0.9);
        //trackCuts->SetMaxRel1PtUncertainty(100.);
        trackCuts->SetPtRange(1., 1000000.0);
        // trackCuts->SetDCAToVertex2D(kTRUE);
        // trackCuts->SetMaxDCAToVertexZ(1E10);
        // trackCuts->SetMaxDCAToVertexXY(1E10);
        jetTask->SetESDCuts(new AliESDtrackCuts(*trackCuts));
    }
    
    //  Final settings, pass to manager and set the containers
    //==============================================================================
    mgr->AddTask(jetTask);
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
    TString contname(combinedName);
    contname += "_histos";
 

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                              TList::Class(),AliAnalysisManager::kOutputContainer,
                                                              Form("%s", AliAnalysisManager::GetCommonFileName()));
 

    
    mgr->ConnectInput  (jetTask, 0,  cinput1 );
    mgr->ConnectOutput (jetTask, 1, coutput1 );

    
    
    
    
    return jetTask;
}


// Setter jet cuts
//==============================================================================

Bool_t DefineCutsTaskpp(AliAnalysisTaskHFJetIPQA *task, Float_t minC, Float_t maxC)
{
    // define cuts for task
    AliRDHFJetsCuts *cuts=task->GetJetCutsHF();
    // jets
    cuts->SetJetRadius(0.4); // this cut does nothing
    cuts->SetMaxEtaJet(0.5);//0.9-R
    cuts->SetMinPtJet(5.);
    cuts->SetMaxPtJet(250.);
    // Set centrality
    cuts->SetMinCentrality(minC);
    cuts->SetMaxCentrality(maxC);
    cuts->SetUsePhysicsSelection(kFALSE);
    cuts->SetOptPileup(1);
    cuts->ConfigurePileupCuts(5,0.8);
    cuts->SetTriggerClass("");
    cuts->SetTriggerMask(AliVEvent::kMB);
    cuts->PrintAll();
    // pPb minbias only
    return kTRUE;
}


//
