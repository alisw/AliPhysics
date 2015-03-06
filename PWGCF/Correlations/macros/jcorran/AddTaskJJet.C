// $Id$

AliJJetTask* AddTaskJJet(
  Int_t       trigger            = AliVEvent::kEMCEJE,
  const char *taskname           = "AliJJetTask",
  int       debug 		 = 1	
)
{  

    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskJJet", "No analysis manager to connect to.");
        return NULL;
    }  

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AddTaskJJet", "This task requires an input event handler");
        return NULL;
    }

    //-------------------------------------------------------
    // Init the task and do settings
    //-------------------------------------------------------

    TString name(taskname);
    //name += trigger;

    //TString tracksName = "PicoTracks";
    TString tracksName = "PicoTracks";
    TString clustersName = "EmcCaloClusters";
    TString clustersCorrName = "CaloClustersCorr";
    TString rhoName = "";


    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");

    const int nJetFinder=6;	
    AliEmcalJetTask* jetFinderTask[nJetFinder];
    double aConeSizes[nJetFinder]={0.4,0.5,0.6,0.4,0.5,0.6};
    int aJetType[nJetFinder]={0,0,0,1,1,1}; // 0 :FullJet  1:Charged
    for (int i=0; i<nJetFinder; i++){
        if(i<3) jetFinderTask[i] = AddTaskEmcalJet(tracksName,clustersCorrName,1,aConeSizes[i],aJetType[i],0.15,0.300,0.005,1,"Jet",5.); // anti-kt
        else jetFinderTask[i] = AddTaskEmcalJet(tracksName,"",1,aConeSizes[i],aJetType[i],0.15,0.300,0.005,1,"Jet",5.); // anti-kt
    }

    char *ntracks = tracksName;
    char *nclusters = clustersCorrName;
    char *nrho = rhoName;
    Printf("name: %s",name.Data());

    AliJJetTask* jtTask = new AliJJetTask(name, nJetFinder);
    jtTask->SetDebug(debug);
    jtTask->SelectCollisionCandidates(trigger);


    AliParticleContainer *trackCont  = jtTask->AddParticleContainer(ntracks);
    trackCont->SetClassName("AliVTrack");
    AliClusterContainer *clusterCont = jtTask->AddClusterContainer(nclusters);

    //char *type="EMCAL";  
    char *type[]={"EMCAL","EMCAL","EMCAL","TPC","TPC","TPC"};  
    AliJetContainer *jetCont[nJetFinder];

    for(int i=0; i<nJetFinder; i++){
        jetCont[i] = jtTask->AddJetContainer(jetFinderTask[i]->GetName(),type[i],aConeSizes[i]);
        if(jetCont[i]) {
            jetCont[i]->SetRhoName(nrho);
            jetCont[i]->ConnectParticleContainer(trackCont);
            if (i<3) jetCont[i]->ConnectClusterContainer(clusterCont);
            jetCont[i]->SetZLeadingCut(0.98,0.98);
            jetCont[i]->SetPercAreaCut(0.6);
            jetCont[i]->SetJetPtCut(5);    
            jetCont[i]->SetLeadingHadronType(0);
        }
    }

    //-------------------------------------------------------
    // Final settings, pass to manager and set the containers
    //-------------------------------------------------------

    mgr->AddTask(jtTask);



    // Create containers for input/output
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
   TString contname(name);
    contname += "_histos";
    cout << "Create container " << contname << endl;
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(contname.Data(),
			TList::Class(),AliAnalysisManager::kOutputContainer,
			Form("%s", AliAnalysisManager::GetCommonFileName()));

	mgr->ConnectInput  (jtTask, 0,  cinput );
	mgr->ConnectOutput (jtTask, 1, coutput ); // MUST HAVE IT, DON"T KNOW WHY ??? maybe from EMCALJET code

	return jtTask;
}
