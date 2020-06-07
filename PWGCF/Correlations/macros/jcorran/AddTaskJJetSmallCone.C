#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C>
#endif
AliJJetTask* AddTaskJJetSmallCone(
    Int_t       trigger           = AliVEvent::kEMCEJE,
    const char *taskname          = "AliJJetTaskSmallCone",
    int         debug             = 1,
    int         doRecoMCPartleJet = 0,  // if do Jet Reconstruction with MC Particle. only for MC
    int         doRecoTrackJet    = 1,  // if do Jet Reconstruction with Reconstructed track. both of MC,Data
    int         isMC              = 0,  // If this event is MC ( both of particle, track )
    int         nR                = 3,  //Number(1-3) of R parameters in order from list 0.4, 0.5, 0.6
    int         doBackgroundEst   = 0,   // If user wants to do background estimation with kt-algorithm
    int         doFullJets        = 1   //If user wants to use full jets, 0=No full jets, 1=Include full jets
    )
{  

  //==============================================================================
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJJetSmallCone", "No analysis manager to connect to.");
    return NULL;
  }  

  //==============================================================================
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJJetSmallCone", "This task requires an input event handler");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(taskname);
  //name += trigger;

  //TString tracksName = "PicoTracks";
  //TString clustersName = "EmcCaloClusters";
  //TString clustersCorrName = "CaloClustersCorr";
  TString tracksName        = "tracks";
  TString clustersCorrName  = "caloClusters";
  TString tracksNameMC      = "mcparticles";
  TString rhoName           = "";

  if(nR < 1) nR = 1;
  if(nR > 3) nR = 3;

  //-------------------------------------------------------
  // LOAD EMCALJetTasks and Variables
  //-------------------------------------------------------
  #if !defined(__CLING__)
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  #endif
  int nTrackJetFinder, nMCParticleJetFinder;
  if(doRecoTrackJet){ //If detector level jets are to be used
    nTrackJetFinder       = (1+doFullJets)*nR;
  }else{
    nTrackJetFinder       = 0;
  }
  if(doRecoMCPartleJet){ //If particle level jets are to be used
    nMCParticleJetFinder  = (1+doFullJets)*nR;
  }else{
    nMCParticleJetFinder  = 0;

  }
  const int nJetFinder              = nMCParticleJetFinder + nTrackJetFinder; //Total number of jet finders is the sum of detector and particle level finders


  //kT jet finders for background estimation
  int nktTrackFinders, nktMCParticleFinders;
  if(doBackgroundEst){
      if(doRecoTrackJet){
          nktTrackFinders = 1;
      }
      else{
          nktTrackFinders = 0;
      }
      if(nMCParticleJetFinder){
          nktMCParticleFinders = 1;
      }
      else{
          nktMCParticleFinders = 0;
      }
  }
  else{
      nktTrackFinders = nktMCParticleFinders =  0;
  }
  const int nktFinders = nktTrackFinders + nktMCParticleFinders;

  int countJetFinder                = 0;  // Counter for number of current Jet Finders

  //-------------------------------------------------------
  // AliJJetTask , AliEmcalJetTask, AliJetContainer
  //-------------------------------------------------------
  AliJJetTask *jtTask = new AliJJetTask(name, nJetFinder + nktFinders);
  jtTask->SetMC( isMC ) ; // Set isMC explicitly.
  jtTask->SetnR( nR ) ; // Set number of jet resolution parameters.
  jtTask->SetIncludeFullJets(doFullJets); //Set whether full jets are included
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
  AliJetContainer::EJetType_t ktJetTypes [1]={  AliJetContainer::kChargedJet};
  TString ktTypes    [1]={"TPC"};// Check if 100% sure about this

  //-------------------------------------------------------
  // Reconstructed Track Jets : Both of Data, MC
  //-------------------------------------------------------
  double bConeSizes[3] = {0.2,0.3,0.4};
  AliJetContainer::EJetType_t bJetType[2] = {AliJetContainer::kFullJet,AliJetContainer::kChargedJet};
  int iTypeStart = doFullJets ? 0 : 1; //Start from index 1 if full jets are to be ignored
  const AliJetContainer::ERecoScheme_t reco  = AliJetContainer::pt_scheme;
  TString bType[2] = {"EMCAL","TPC"};
  int iEnd;
  AliParticleContainer *trackCont;
  AliClusterContainer *clusterCont;
  if( doRecoTrackJet ){
    //================= Variables for array
    int iStart = countJetFinder;
    iEnd   = countJetFinder += nTrackJetFinder;
    //================= AddTaskEmcalJet

    //================= Containers : Track, Cluster
    AliTrackContainer* partCont = new AliTrackContainer(tracksName);
    partCont->SetFilterHybridTracks(kFALSE);
    //partCont->SetFilterHybridTracks(kTRUE);
    trackCont = partCont;
    trackCont->SetClassName("AliVTrack");
    //AliParticleContainer *trackCont  = jtTask->AddParticleContainer( tracksName );
    if (trackCont) jtTask->AdoptParticleContainer(trackCont);
    clusterCont = jtTask->AddClusterContainer( clustersCorrName );
    clusterCont->SetClassName("AliAODCaloCluster");

    //================= LOOP for configuraion
    //for(int i=iStart ; i<iEnd; i++){
    for(int itype = iTypeStart; itype < 2; itype++){
      for(int iR = 0; iR < nR; iR++){
        int iF = iStart + (itype-iTypeStart)*nR + iR; //Jet Finder index
        //== Variables
        double consize = bConeSizes[iR];
        AliJetContainer::EJetType_t jettype = bJetType[itype];
        TString type = bType[itype];

        //== JetFinderTask, JetContainer
        TString _clustersCorrName = ( type == "EMCAL" ? clustersCorrName : "" ); // Check if it is EMCAL
        jetFinderTask[iF] = AddTaskEmcalJet( tracksName, _clustersCorrName, AliJetContainer::antikt_algorithm, consize, jettype, 0.15, 0.300, 0.005, reco, "Jet", 5. ); // anti-kt
        jetFinderTask[iF]->SelectCollisionCandidates(trigger);
        jtTask->SetTrackOrMCParticle( iF, AliJJetTask::kJRecoTrack ); 
        jtTask->SetConeSize( iF, consize );
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
      } //  iR loop
    } // itype loop
  } // if doRecoTrackJet

  //-------------------------------------------------------
  // MC Particles Jets
  //-------------------------------------------------------
  AliMCParticleContainer *mcTrackCont;
  if( doRecoMCPartleJet ){
    //================= Variables for array
    int iStart = countJetFinder;
    iEnd   = countJetFinder += nMCParticleJetFinder; // WARN: this is a real end + 1

    //================= AddTaskEmcalJet
    TString nrho = rhoName;

    //================= Containers : MC Particle
    mcTrackCont = jtTask->AddMCParticleContainer(tracksNameMC);
    mcTrackCont->SelectPhysicalPrimaries(kTRUE);
    mcTrackCont->SetClassName("AliAODMCParticle");

    //================= LOOP for configuraion
    for(int itype = iTypeStart; itype < 2; itype++){
      for(int iR = 0; iR < nR; iR++){
        int iF = iStart + (itype-iTypeStart)*nR + iR;
        double consizeMC = bConeSizes[iR];
        AliJetContainer::EJetType_t jettype = bJetType[itype];
        TString type = bType[itype];

        //== JetFinderTask, JetContainer
        TString _clustersCorrName = ( type == "EMCAL" ? clustersCorrName : "" ); // Check if it is EMCAL 
        jetFinderTask[iF] = AddTaskEmcalJet( tracksNameMC, "", AliJetContainer::antikt_algorithm, consizeMC, jettype, 0.15, 0.300, 0.005, reco, "Jet", 5. ); // anti-kt
        jetFinderTask[iF]->SelectCollisionCandidates(trigger);
        cout << jetFinderTask[iF]->GetName() << endl;
        //jetFinderTask[i]->GetParticleContainer(0)->SelectPhysicalPrimaries(kTRUE); 
        jtTask->SetTrackOrMCParticle( iF, AliJJetTask::kJMCParticle ); 
        jtTask->SetConeSize( iF, consizeMC );
        jetCont[iF] = jtTask->AddJetContainer( jetFinderTask[iF]->GetName(), type, consizeMC );

        if(jetCont[iF]) {
          jetCont[iF]->SetRhoName( rhoName );
          jetCont[iF]->ConnectParticleContainer( mcTrackCont );
          jetCont[iF]->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
          jetCont[iF]->SetPercAreaCut( 0.6 );
          jetCont[iF]->SetJetPtCut( 5 );
          jetCont[iF]->SetLeadingHadronType( 0 );
        }
      } // iR loop
    } // itype loop
  } // if doRecoMCPartleJet

  //-------------------------------------------------------
  // kt-algorithm
  //-------------------------------------------------------
  if (doBackgroundEst) {
    double ktConeSize = ktConeSizes[0];
    AliJetContainer::EJetType_t ktJetType  = ktJetTypes[0];
    TString ktType = ktTypes[0];
    if (doRecoTrackJet) {
      TString _clustersCorrName = ( ktType == "EMCAL" ? clustersCorrName : "" );
      ktFinderTask[0] = AddTaskEmcalJet( tracksName, _clustersCorrName, AliJetContainer::kt_algorithm, ktConeSize, ktJetType, 0.15, 0.300, 0.005, reco, "Jet", 0.0 ); // kt
      ktFinderTask[0]->SelectCollisionCandidates(trigger);
      jtTask->SetTrackOrMCParticle( iEnd, AliJJetTask::kJRecoTrack );
      jtTask->SetConeSize( iEnd, ktConeSize );
      cout << ktFinderTask[0]->GetName() << endl;
      ktCont[0] = jtTask->AddJetContainer( ktFinderTask[0]->GetName(), ktType, ktConeSize );

      if( ktCont[0] ) {
        ktCont[0]->SetRhoName( rhoName );
        ktCont[0]->ConnectParticleContainer( trackCont );
        if ( ktType == "EMCAL" ) ktCont[0]->ConnectClusterContainer( clusterCont );
          //ktCont[0]->SetZLeadingCut( 0.98, 0.98 );
          //ktCont[0]->SetPercAreaCut( 0.6 );
          //ktCont[0]->SetJetPtCut( 5 );
          //ktCont[0]->SetLeadingHadronType( 0 );
      }
    } // if doRecoTrackJet
    if (doRecoMCPartleJet) {
      TString _clustersCorrName = ( ktType == "EMCAL" ? clustersCorrName : "" );
      ktFinderTask[1] = AddTaskEmcalJet( tracksNameMC, "", AliJetContainer::kt_algorithm, ktConeSize, ktJetType, 0.15, 0.300, 0.005, reco, "Jet", 0.0 ); // kt
      ktFinderTask[1]->SelectCollisionCandidates(trigger);
      jtTask->SetTrackOrMCParticle( iEnd+1, AliJJetTask::kJMCParticle );
      jtTask->SetConeSize( iEnd+1, ktConeSize );
      ktCont[1] = jtTask->AddJetContainer( ktFinderTask[1]->GetName(), ktType, ktConeSize );

      if( ktCont[1] ) {
        ktCont[1]->SetRhoName( rhoName );
        ktCont[1]->ConnectParticleContainer( mcTrackCont );
        if ( ktType == "EMCAL" ) ktCont[0]->ConnectClusterContainer( clusterCont );
          //ktCont[1]->SetZLeadingCut( 0.98, 0.98 );
          //ktCont[1]->SetPercAreaCut( 0.6 );
          //ktCont[1]->SetJetPtCut( 5 );
          //ktCont[1]->SetLeadingHadronType( 0 );
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
  mgr->ConnectOutput (jtTask, 1, coutput );

  return jtTask;
  }
