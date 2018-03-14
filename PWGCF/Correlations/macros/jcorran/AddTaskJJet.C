AliJJetTask* AddTaskJJet(
    Int_t       trigger           = AliVEvent::kEMCEJE,
    const char *taskname          = "AliJJetTask",
    int         debug             = 1,
    int         doRecoMCPartleJet = 0,  // if do Jet Reconstruction with MC Particle. only for MC
    int         doRecoTrackJet    = 1,  // if do Jet Reconstruction with Reconstructed track. both of MC,Data
    int         isMC              = 0,  // If this event is MC ( both of particle, track )
    int         nR                = 3,  //Number(1-3) of R parameters in order from list 0.4, 0.5, 0.6
    int         doBackgroundEst   = 0   // If user wants to do background estimation with kt-algorithm
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

  if(doBackgroundEst){
      if(doRecoTrackJet){
          const int nktTrackFinders = 1;
      }
      else{
          const int nktTrackFinders = 0;
      }
      if(nMCParticleJetFinder){
          const int nktMCParticleFinders = 1;
      }
      else{
          const int nktMCParticleFinders = 0;
      }
      const int nktFinders = nktTrackFinders + nktMCParticleFinders;
  }
  else{
      const int nktFinders = 0;
  }


  int countJetFinder                = 0;  // Counter for number of current Jet Finders

  //-------------------------------------------------------
  // AliJJetTask , AliEmcalJetTask, AliJetContainer
  //-------------------------------------------------------
  AliJJetTask *jtTask = new AliJJetTask(name, nJetFinder + nktFinders);
  jtTask->SetMC( isMC ) ; // Set isMC explicitly. 
  jtTask->SetnR( nR ) ; // Set number of jet resolution parameters. 
  jtTask->Setnkt( nktFinders ) ; // Set information about if kt algorithm reconstruction was done or not. 
  jtTask->SetDebug( debug );
  jtTask->SelectCollisionCandidates( trigger ); // WARNING: default is AliVEvent::kEMCEJE. Double check it!

  AliEmcalJetTask *jetFinderTask[nJetFinder];
  for( int i=0;i<nJetFinder;i++ ) jetFinderTask[i] = NULL;

  AliJetContainer *jetCont[nJetFinder];
  for( int i=0;i<nJetFinder;i++ ) jetCont[i] = NULL;

  AliEmcalJetTask *ktFinderTask[nktFinders];
  for( int i=0;i<nktFinders;i++) ktFinderTask[i] = NULL;
  AliJetContainer *ktCont[nktFinders]; //0 for Data and 1 for MC. Implementation for different R left for later.
  for( int i=0;i<nktFinders;i++) ktCont[i] = NULL; 
  double  ktConeSizes[1]={  0.2};
  int     ktJetTypes [1]={    1};// 0:FullJet 1:Charged 
  TString ktTypes    [1]={"TPC"};// Check if 100% sure about this

  //-------------------------------------------------------
  // Reconstructed Track Jets : Both of Data, MC
  //-------------------------------------------------------
  if( doRecoTrackJet ){
    //================= Variables for array
    int iStart = countJetFinder;
    int iEnd   = countJetFinder += nTrackJetFinder;
    //================= AddTaskEmcalJet
    double bConeSizes[3] = {0.4,0.3,0.2};
    int bJetType[2] = {0,1};
    TString bType[2] = {"EMCAL","TPC"};

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
        //== Variables
        double consize = bConeSizes[iR];
        int jettype = bJetType[itype];
        TString type = bType[itype];

        //== JetFinderTask, JetContainer
        TString _clustersCorrName = ( type == "EMCAL" ? clustersCorrName : "" ); // Check if it is EMCAL
        jetFinderTask[iF] = AddTaskEmcalJet( tracksName, _clustersCorrName, AliJetContainer::antikt_algorithm, consize, jettype, 0.15, 0.300, 0.005, 1, "Jet", 5. ); // anti-kt
        jetFinderTask[iF]->SelectCollisionCandidates(trigger);
        jtTask->SetTrackOrMCParticle( iF, AliJJetTask::kJRecoTrack ); 
        cout << jetFinderTask[iF]->GetName() << endl;
        jetCont[iF] = jtTask->AddJetContainer( jetFinderTask[iF]->GetName(), type, consize ); 

        if( jetCont[iF] ) {
          jetCont[iF]->SetRhoName( rhoName );
          jetCont[iF]->ConnectParticleContainer( trackCont );
          if ( type == "EMCAL" ) jetCont[iF]->ConnectClusterContainer( clusterCont );
            jetCont[iF]->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
            jetCont[iF]->SetPercAreaCut( 0.6 );
            jetCont[iF]->SetJetPtCut( 5 );    
            jetCont[iF]->SetLeadingHadronType( 0 );
        }
      } // for iR
    } // for itype
  } // if doRecoTrackJet

  //-------------------------------------------------------
  // MC Particles Jets
  //-------------------------------------------------------
  if( doRecoMCPartleJet ){
    //================= Variables for array
    int iStart = countJetFinder;
    int iEnd   = countJetFinder += nMCParticleJetFinder; // WARN: this is a real end + 1

    //================= AddTaskEmcalJet
    double bConeSizes[3] = {0.4,0.3,0.2};
    int bJetType[2] = {0,1};
    TString bType[2] = {"EMCAL","TPC"};
    char *nrho = rhoName;

    //================= Containers : MC Particle
    AliMCParticleContainer *mcTrackCont = jtTask->AddMCParticleContainer(tracksNameMC);
    mcTrackCont->SelectPhysicalPrimaries(kTRUE);
    mcTrackCont->SetClassName("AliAODMCParticle");

    //================= LOOP for configuraion
    for(int itype = 0; itype < 2; itype++){
      for(int iR = 0; iR < nR; iR++){
        int iF = iStart + itype*nR + iR;
        double consizeMC = bConeSizes[iR];
        int jettype = bJetType[itype];
        TString type = bType[itype];

        //== JetFinderTask, JetContainer
        TString _clustersCorrName = ( type == "EMCAL" ? clustersCorrName : "" ); // Check if it is EMCAL 
        jetFinderTask[iF] = AddTaskEmcalJet( tracksNameMC, "", AliJetContainer::antikt_algorithm, consizeMC, jettype, 0.15, 0.300, 0.005, 1, "Jet", 5. ); // anti-kt
        jetFinderTask[iF]->SelectCollisionCandidates(trigger);
        cout << jetFinderTask[iF]->GetName() << endl;
        //jetFinderTask[i]->GetParticleContainer(0)->SelectPhysicalPrimaries(kTRUE); 
        jtTask->SetTrackOrMCParticle( iF, AliJJetTask::kJMCParticle ); 
        jetCont[iF] = jtTask->AddJetContainer( jetFinderTask[iF]->GetName(), type, consizeMC );

        if(jetCont[iF]) {
          jetCont[iF]->SetRhoName( rhoName );
          jetCont[iF]->ConnectParticleContainer( mcTrackCont );
          jetCont[iF]->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
          jetCont[iF]->SetPercAreaCut( 0.6 );
          jetCont[iF]->SetJetPtCut( 5 );    
          jetCont[iF]->SetLeadingHadronType( 0 );
        }
      } // for iR
    } // for itype
  } // if doRecoMCPartleJet

  //-------------------------------------------------------
  // kt-algorithm
  //-------------------------------------------------------
  if (doBackgroundEst) {
    double ktConeSize = ktConeSizes[0];
    int ktJetType  = ktJetTypes[0];
    TString ktType = ktTypes[0];
    if (doRecoTrackJet) {
      TString _clustersCorrName = ( ktType == "EMCAL" ? clustersCorrName : "" );
      ktFinderTask[0] = AddTaskEmcalJet( tracksName, _clustersCorrName, AliJetContainer::kt_algorithm, ktConeSize, ktJetType, 0.15, 0.300, 0.005, 1, "Jet", 5. ); // kt
      ktFinderTask[0]->SelectCollisionCandidates(trigger);
      jtTask->SetTrackOrMCParticle( iEnd, AliJJetTask::kJRecoTrack );
      cout << ktFinderTask[0]->GetName() << endl;
      ktCont[0] = jtTask->AddJetContainer( ktFinderTask[0]->GetName(), ktType, ktConeSize );

      if( ktCont[0] ) {
        ktCont[0]->SetRhoName( rhoName );
        ktCont[0]->ConnectParticleContainer( trackCont );
        if ( ktType == "EMCAL" ) ktCont[0]->ConnectClusterContainer( clusterCont );
          ktCont[0]->SetZLeadingCut( 0.98, 0.98 );
          ktCont[0]->SetPercAreaCut( 0.6 );
          ktCont[0]->SetJetPtCut( 5 );    
          ktCont[0]->SetLeadingHadronType( 0 );
      }
    } // if doRecoTrackJet
    if (doRecoMCPartleJet) {
      TString _clustersCorrName = ( ktType == "EMCAL" ? clustersCorrName : "" );
      ktFinderTask[1] = AddTaskEmcalJet( tracksNameMC, "", AliJetContainer::kt_algorithm, ktConeSize, ktJetType, 0.15, 0.300, 0.005, 1, "Jet", 5. ); // kt
      ktFinderTask[1]->SelectCollisionCandidates(trigger);
      jtTask->SetTrackOrMCParticle( iEnd+1, AliJJetTask::kJMCParticle );
      ktCont[1] = jtTask->AddJetContainer( ktFinderTask[1]->GetName(), ktType, ktConeSize );

      if( ktCont[1] ) {
        ktCont[1]->SetRhoName( rhoName );
        ktCont[1]->ConnectParticleContainer( mcTrackCont );
        if ( ktType == "EMCAL" ) ktCont[0]->ConnectClusterContainer( clusterCont );
          ktCont[1]->SetZLeadingCut( 0.98, 0.98 );
          ktCont[1]->SetPercAreaCut( 0.6 );
          ktCont[1]->SetJetPtCut( 5 );    
          ktCont[1]->SetLeadingHadronType( 0 );
      }
    } // if doRecoTrackJet
  } // if doBackgroundEst

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  Printf("name: %s",name.Data());

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
