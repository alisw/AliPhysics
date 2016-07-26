// $Id$

AliJJetTask* AddTaskJJet(
    // FIXME: Really need default values? Is this needed by AliROOT? Or, let's remove it to reduce mistakes
    Int_t       trigger           = AliVEvent::kEMCEJE,
    const char *taskname          = "AliJJetTask",
    int         debug             = 1,
    int         doRecoMCPartleJet = 0,  // if do Jet Reconstruction with MC Particle. only for MC
    int         doRecoTrackJet    = 1,  // if do Jet Reconstruction with Reconstructed track. both of MC,Data
    int         isMC              = 0   // If this event is MC ( both of particle, track )
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


  //-------------------------------------------------------
  // LOAD EMCALJetTasks and Variables
  //-------------------------------------------------------
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  if(doRecoTrackJet){
    const int nTrackJetFinder       = 6; 
  }else{
    const int nTrackJetFinder       = 0;
  }
  if(doRecoMCPartleJet){
    const int nMCParticleJetFinder  = 3;
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
    double   aConeSizes[nTrackJetFinder]={     0.4,     0.5,     0.6,   0.4,   0.5,   0.6 };
    int      aJetType  [nTrackJetFinder]={       0,       0,       0,     1,     1,     1 }; // 0 :FullJet  1:Charged
    TString  aType     [nTrackJetFinder]={ "EMCAL", "EMCAL", "EMCAL", "TPC", "TPC", "TPC" };  

    //================= Containers : Track, Cluster
    AliParticleContainer *trackCont  = jtTask->AddParticleContainer( tracksName );
    trackCont->SetClassName("AliVTrack");
    AliClusterContainer *clusterCont = jtTask->AddClusterContainer( clustersCorrName );
    clusterCont->SetClassName("AliAODCaloCluster");

    //================= LOOP for configuraion
    for(int i=iStart ; i<iEnd; i++){
      //== Variables
      //int iHere = iStart - countJetFinder; // array index based 0 in a block
      int iHere = i-iStart; // array index based 0 in a block
      double consize = aConeSizes[iHere];
      int jettype = aJetType  [iHere];
      int type    = aType     [iHere];

      //== JetFinderTask, JetContainer
      TString _clustersCorrName = ( type == "EMCAL" ? clustersCorrName : "" ); // Check if it is EMCAL
      jetFinderTask[i] = AddTaskEmcalJet( tracksName, _clustersCorrName, 1, consize, jettype, 0.15, 0.300, 0.005, 1, "Jet", 5. ); // anti-kt
      jtTask->SetTrackOrMCParticle( i, AliJJetTask::kJRecoTrack ); // TODO: Move this to AddJetContainer
      cout << jetFinderTask[i]->GetName() << endl;
      jetCont[i] = jtTask->AddJetContainer( jetFinderTask[i]->GetName(), type, consize ); // TODO: SetTrackOrMCParticle here

      if( jetCont[i] ) {
        jetCont[i]->SetRhoName( rhoName );
        jetCont[i]->ConnectParticleContainer( trackCont );
        if ( type == "EMCAL" ) jetCont[i]->ConnectClusterContainer( clusterCont );
        jetCont[i]->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
        jetCont[i]->SetPercAreaCut( 0.6 );
        jetCont[i]->SetJetPtCut( 5 );    
        jetCont[i]->SetLeadingHadronType( 0 );
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
    double   aConeSizesMC[nMCParticleJetFinder] = {   0.4,   0.5,   0.6 };
    int      aJetTypeMC  [nMCParticleJetFinder] = {     1,     1,     1 }; // 0 :FullJet  1:Charged
    TString *aTypeMC     [nMCParticleJetFinder] = { "TPC", "TPC", "TPC" };  

    char *nrho = rhoName;

    //================= Containers : MC Particle
    AliMCParticleContainer *mcTrackCont = jtTask->AddMCParticleContainer(tracksNameMC);
    mcTrackCont->SelectPhysicalPrimaries(kTRUE);
    mcTrackCont->SetClassName("AliAODMCParticle");

    //================= LOOP for configuraion
    for(int i=iStart ; i<iEnd ; i++){
      //== Variables
      //int iHere = iStart - countJetFinder; // array index based 0 in a block
      int iHere = i - iStart; // array index based 0 in a block
      double consizeMC = aConeSizesMC[iHere];
      int jettype = aJetTypeMC  [iHere];
      int type    = aTypeMC     [iHere];

      //== JetFinderTask, JetContainer
      TString _clustersCorrName = ( type == "EMCAL" ? clustersCorrName : "" ); // Check if it is EMCAL // FIXME: Good for MC particles also? //No clusters in mc particles?
      jetFinderTask[i] = AddTaskEmcalJet( tracksNameMC, _clustersCorrName, 1, consizeMC, jettype, 0.15, 0.300, 0.005, 1, "Jet", 5. ); // anti-kt
      cout << jetFinderTask[i]->GetName() << endl;
      //jetFinderTask[i]->GetParticleContainer(0)->SelectPhysicalPrimaries(kTRUE); 
      // FIXME: DUPLICATION? not needed here, I think // FIXME: Why before ConnectParticleContainer? //I don't think the order matters
      // FIXME: Is this means, we don't need to make particle container manually?
      jtTask->SetTrackOrMCParticle( i, AliJJetTask::kJMCParticle ); // TODO: Move this to AddJetContainer FIXME: Maybe not possible in easy
      jetCont[i] = jtTask->AddJetContainer( jetFinderTask[i]->GetName(), type, consizeMC ); // TODO: SetTrackOrMCParticle here

      if(jetCont[i]) {
        jetCont[i]->SetRhoName( rhoName );
        jetCont[i]->ConnectParticleContainer( mcTrackCont );
        // if ( aType[iHere] == "EMCAL" ) jetCont[i]->ConnectClusterContainer( clusterCont );
        jetCont[i]->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
        jetCont[i]->SetPercAreaCut( 0.6 );
        jetCont[i]->SetJetPtCut( 5 );    
        jetCont[i]->SetLeadingHadronType( 0 );
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
