#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C>
#endif
AliJJetTask* AddTaskJJetSingle(
    Int_t       trigger           = AliVEvent::kEMCEJE,
    const char *taskname          = "AliJJetTask",
    int         debug             = 1,
    int         isMC              = 0,    // 0:Reco 1:MC(kineOnly)
    double      consize           = 0.4,  //
    int         doFullJets        = 1   //If user wants to use full jets, 0=No full jets, 1=Include full jets
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
  // isMC 0:1 
  TString tracksName[2]        = {"usedefault","usedefault"};//"mcparticles"};
  TString clustersCorrName[2]  = {"usedefault","mcparticles"};
  TString DataClassName[2]  = {"AliAODCaloCluster","AliAODMCParticle"};
  TString rhoName              = "rho";
  //-------------------------------------------------------
  // LOAD EMCALJetTasks and Variables
  //-------------------------------------------------------
  #if !defined(__CLING__)
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  #endif 
  //-------------------------------------------------------
  // AliJJetTask , AliEmcalJetTask, AliJetContainer
  //-------------------------------------------------------
  AliJJetTask *jtTask = new AliJJetTask(name, 2);
  jtTask->SetMC( isMC ) ; // Set isMC explicitly.
  jtTask->SetnR( 1 ) ; // Set number of jet resolution parameters.
  jtTask->SetIncludeFullJets(doFullJets); //Set whether full jets are included
  jtTask->Setnkt( 1 ) ; // Set information about if kt algorithm reconstruction was done or not.
  jtTask->SetDebug( debug );
  jtTask->SelectCollisionCandidates( trigger ); // WARNING: default is AliVEvent::kEMCEJE. Double check it!

  //EMCAL jets
  AliEmcalJetTask *jetFinderTask;
  AliJetContainer *jetCont;
  AliEmcalJetTask *ktFinderTask;
  AliJetContainer *ktCont;
  double  ktConeSize = 0.2;
  //-------------------------------------------------------
  // Reconstructed Track Jets : Both of Data, MC
  //-------------------------------------------------------
  AliJetContainer::EJetType_t bJetType[2] = {AliJetContainer::kFullJet,AliJetContainer::kChargedJet};
  int iType = doFullJets ? 0 : 1; //Start from index 1 if full jets are to be ignored
  const AliJetContainer::ERecoScheme_t reco  = AliJetContainer::pt_scheme;
  TString bType[2] = {"EMCALfid","TPCfid"};
  AliJetContainer::EJetType_t jettype = bJetType[iType];
  TString type = bType[iType];
  
  //== JetFinderTask, JetContainer
  jetFinderTask = AddTaskEmcalJet( tracksName[isMC], clustersCorrName[isMC], AliJetContainer::antikt_algorithm, consize, jettype, 0.15, 0.300, 0.005, reco, "Jet", 0.15 ); // anti-kt
  jetFinderTask->SelectCollisionCandidates(trigger);
  jtTask->SetTrackOrMCParticle( isMC, isMC ? AliJJetTask::kJMCParticle : AliJJetTask::kJRecoTrack );  //check with AliJJet??
  jtTask->SetConeSize(0 , consize );
  cout << jetFinderTask->GetName() << endl;
  jetCont = jtTask->AddJetContainer( jetFinderTask->GetName(), type, consize ); 
  
//-------------------------------------------------------
  // kt-algorithm
  //-------------------------------------------------------
  ktFinderTask = AddTaskEmcalJet( tracksName[isMC], clustersCorrName[isMC], AliJetContainer::kt_algorithm, ktConeSize, jettype, 0.15, 0.300, 0.005, reco, "Jet", 0.15 ); // kt
  ktFinderTask->SelectCollisionCandidates(trigger);
  jtTask->SetTrackOrMCParticle(isMC, isMC ? AliJJetTask::kJMCParticle : AliJJetTask::kJRecoTrack ); // CHECK!!!
  jtTask->SetConeSize( 1, ktConeSize );
  cout << ktFinderTask->GetName() << endl;
  ktCont = jtTask->AddJetContainer( ktFinderTask->GetName(), jettype, ktConeSize );
  

  //------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  Printf("name: %s",name.Data());

  mgr->AddTask(jtTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
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
