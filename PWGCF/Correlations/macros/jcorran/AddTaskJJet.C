AliJJetTask* AddTaskJJet(
    // FIXME: Really need default values? Is this needed by AliROOT? Or, let's remove it to reduce mistakes
    Int_t       trigger           = AliVEvent::kEMCEJE,
    const char *taskname          = "AliJJetTask",
    int         debug             = 1,
    int         doRecoMCPartleJet = 0,  // if do Jet Reconstruction with MC Particle. only for MC
    int         doRecoTrackJet    = 1,  // if do Jet Reconstruction with Reconstructed track. both of MC,Data
    int         isMC              = 0,   // If this event is MC ( both of particle, track )
    int         nR                = 3 //Number(1-3) of R parameters in order from list 0.4, 0.5, 0.6
    )
{  

  //==============================================================================
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJJet", "No analysis manager to connect to.");
    return NULL;
  }  

  //==============================================================================
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
  //name += trigger; // FIXME: Why not? Could be

  //TString tracksName = "PicoTracks";
  //TString tracksName = "PicoTracks";
  //TString clustersName = "EmcCaloClusters";
  //TString clustersCorrName = "CaloClustersCorr";
  TString tracksName        = "tracks";
  TString clustersCorrName  = "caloClusters"; 
  TString tracksNameMC      = "mcparticles";  //Check these
  TString rhoName           = "";             // FIXME: Comment me  

  if(nR < 1) nR = 1;
  if(nR > 3) nR = 3;

  //-------------------------------------------------------
  // LOAD EMCALJetTasks and Variables
  //-------------------------------------------------------
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  if(doRecoTrackJet){
    const int nTrackJetFinder       = 2*nR; 
  }else{
    const int nTrackJetFinder       = 0;
  }
  if(doRecoMCPartleJet){
    const int nMCParticleJetFinder  = 2*nR;
  }else{
    const int nMCParticleJetFinder  = 0;

  }
  const int nJetFinder              = nMCParticleJetFinder + nTrackJetFinder;


  int countJetFinder                = 0;  // Counter for number of current Jet Finders

  //-------------------------------------------------------
  // AliJJetTask , AliEmcalJetTask, AliJetContainer
  //-------------------------------------------------------
  AliJJetTask *jtTask = new AliJJetTask(name, nJetFinder);
  jtTask->SetMC( isMC ) ; // Set isMC explicitly. 
  jtTask->SetDebug( debug );
  jtTask->SelectCollisionCandidates( trigger ); // WARNING: default is AliVEvent::kEMCEJE. Double check it!

  AliEmcalJetTask *jetFinderTask[nJetFinder];
  for( int i=0;i<nJetFinder;i++ ) jetFinderTask[i] = NULL;

  AliJetContainer *jetCont[nJetFinder];
  for( int i=0;i<nJetFinder;i++ ) jetCont[i] = NULL;


  //-------------------------------------------------------
  // Reconstructed Track Jets : Both of Data, MC
  //-------------------------------------------------------
  if( doRecoTrackJet ){
    //================= Variables for array
    int iStart = countJetFinder;
    int iEnd   = countJetFinder += nTrackJetFinder;
    //================= AddTaskEmcalJet
    // We have 6 type of recos
    double   aConeSizes[6]={     0.4,     0.5,     0.6,   0.4,   0.5,   0.6 };
    int      aJetType  [6]={       0,       0,       0,     1,     1,     1 }; // 0 :FullJet  1:Charged
    TString  aType     [6]={ "EMCAL", "EMCAL", "EMCAL", "TPC", "TPC", "TPC" };

    double bConeSizes[3] = {0.4,0.5,0.6};
    int bJetType[2] = {0,1};
    int bType[2] = {"EMCAL","TPC"};

    //================= Containers : Track, Cluster
    AliTrackContainer* partCont = new AliTrackContainer(tracksName);
    partCont->SetFilterHybridTracks(kFALSE);
    //partCont->SetFilterHybridTracks(kTRUE);
    AliParticleContainer *trackCont = partCont; 
    trackCont->SetClassName("AliVTrack");
    //AliParticleContainer *trackCont  = jtTask->AddParticleContainer( tracksName );
    if (trackCont) jtTask->AdoptParticleContainer(trackCont);
    AliClusterContainer *clusterCont = jtTask->AddClusterContainer( clustersCorrName );
    clusterCont->SetClassName("AliAODCaloCluster");

    //================= LOOP for configuraion
    //for(int i=iStart ; i<iEnd; i++){
    for(int itype = 0; itype < 2; itype++){
      for(int iR = 0; iR < nR; iR++){
        int iF = iStart + itype*nR + iR;
        cout << "iF: " << iF << endl;
        //== Variables
        //int iHere = iStart - countJetFinder; // array index based 0 in a block
        int iHere = i-iStart; // array index based 0 in a block
        //double consize = aConeSizes[iHere];
        //int jettype = aJetType  [iHere];
        //int type    = aType     [iHere];
        double consize = bConeSizes[iR];
        int jettype = bJetType[itype];
        int type = bType[itype];

        //== JetFinderTask, JetContainer
        TString _clustersCorrName = ( type == "EMCAL" ? clustersCorrName : "" ); // Check if it is EMCAL
        jetFinderTask[iF] = AddTaskEmcalJet( tracksName, _clustersCorrName, 1, consize, jettype, 0.15, 0.300, 0.005, 1, "Jet", 5. ); // anti-kt
        jtTask->SetTrackOrMCParticle( iF, AliJJetTask::kJRecoTrack ); // TODO: Move this to AddJetContainer
        cout << jetFinderTask[iF]->GetName() << endl;
        jetCont[iF] = jtTask->AddJetContainer( jetFinderTask[iF]->GetName(), type, consize ); // TODO: SetTrackOrMCParticle here

        if( jetCont[iF] ) {
          jetCont[iF]->SetRhoName( rhoName );
          jetCont[iF]->ConnectParticleContainer( trackCont );
          if ( type == "EMCAL" ) jetCont[iF]->ConnectClusterContainer( clusterCont );
          jetCont[iF]->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
          jetCont[iF]->SetPercAreaCut( 0.6 );
          jetCont[iF]->SetJetPtCut( 5 );    
          jetCont[iF]->SetLeadingHadronType( 0 );
        }
      }
    } // for i
  } // if doRecoTrackJet

  //-------------------------------------------------------
  // MC Particles Jets
  //-------------------------------------------------------
  if( doRecoMCPartleJet ){
    //================= Variables for array
    int iStart = countJetFinder;
    int iEnd   = countJetFinder += nMCParticleJetFinder; // WARN: this is a real end + 1

    //================= AddTaskEmcalJet
    // We have 3 type of recos
    double   aConeSizesMC[6]={     0.4,     0.5,     0.6,   0.4,   0.5,   0.6 };
    int      aJetTypeMC  [6]={       0,       0,       0,     1,     1,     1 }; // 0 :FullJet  1:Charged
    TString  aTypeMC     [6]={ "TPC", "TPC", "TPC", "TPC", "TPC", "TPC" };  

    double bConeSizes[3] = {0.4,0.5,0.6};
    int bJetType[2] = {0,1};
    int bType[2] = {"EMCAL","TPC"};

    char *nrho = rhoName;

    //================= Containers : MC Particle
    AliMCParticleContainer *mcTrackCont = jtTask->AddMCParticleContainer(tracksNameMC);
    mcTrackCont->SelectPhysicalPrimaries(kTRUE);
    mcTrackCont->SetClassName("AliAODMCParticle");

    //================= LOOP for configuraion
    for(int itype = 0; itype < 2; itype++){
      for(int iR = 0; iR < nR; iR++){
        int iF = iStart + itype*nR + iR;
        cout << "iF: " << iF << endl;
        //== Variables
        //int iHere = iStart - countJetFinder; // array index based 0 in a block
        //int iHere = i - iStart; // array index based 0 in a block
        //double consizeMC = aConeSizesMC[iHere];
        //int jettype = aJetTypeMC  [iHere];
        //int type    = aTypeMC     [iHere];
        double consizeMC = bConeSizes[iR];
        int jettype = bJetType[itype];
        int type = bType[itype];

        //== JetFinderTask, JetContainer
        TString _clustersCorrName = ( type == "EMCAL" ? clustersCorrName : "" ); // Check if it is EMCAL // FIXME: Good for MC particles also? //No clusters in mc particles?
        jetFinderTask[iF] = AddTaskEmcalJet( tracksNameMC, "", 1, consizeMC, jettype, 0.15, 0.300, 0.005, 1, "Jet", 5. ); // anti-kt
        cout << jetFinderTask[iF]->GetName() << endl;
        //jetFinderTask[i]->GetParticleContainer(0)->SelectPhysicalPrimaries(kTRUE); 
        // FIXME: DUPLICATION? not needed here, I think // FIXME: Why before ConnectParticleContainer? //I don't think the order matters
        // FIXME: Is this means, we don't need to make particle container manually?
        jtTask->SetTrackOrMCParticle( iF, AliJJetTask::kJMCParticle ); // TODO: Move this to AddJetContainer FIXME: Maybe not possible in easy
        jetCont[iF] = jtTask->AddJetContainer( jetFinderTask[iF]->GetName(), type, consizeMC ); // TODO: SetTrackOrMCParticle here

        if(jetCont[iF]) {
          jetCont[iF]->SetRhoName( rhoName );
          jetCont[iF]->ConnectParticleContainer( mcTrackCont );
          // if ( aType[iHere] == "EMCAL" ) jetCont[i]->ConnectClusterContainer( clusterCont );
          jetCont[iF]->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
          jetCont[iF]->SetPercAreaCut( 0.6 );
          jetCont[iF]->SetJetPtCut( 5 );    
          jetCont[iF]->SetLeadingHadronType( 0 );
        }
      }
    } // for i
  } // if doRecoMCPartleJet

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  Printf("name: %s",name.Data());

  cout << "Add task" << endl;
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
