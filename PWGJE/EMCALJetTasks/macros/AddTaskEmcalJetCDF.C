/// \file AddTaskEmcalJetCDF.C
/// \brief Adds a AliAnalysisTaskEmcalJetCDF analysis task and coresponding containers
///
/// Analysis of leading jets distribution of pt and multiplicity, R distribution
/// N80 and Pt80 and Toward, Transverse, Away UE histos
///
/// \author Adrian SEVCENCO <Adrian.Sevcenco@cern.ch>, Institute of Space Science, Romania
/// \date Mar 19, 2015

/// Add a AliAnalysisTaskEmcalJetCDF task - detailed signature
/// \param const char *ntracks
/// \param const char *nclusters
/// \param const char *njets
/// \param const char *nrho
/// \param Double_t jetradius
/// \param Double_t jetptcut
/// \param Double_t jetareacut
/// \param const char *type ; either TPC, EMCAL or USER
/// \param Int_t leadhadtype ; 0 = charged, 1 = neutral, 2 = both
/// \param const char *taskname
/// \return AliAnalysisTaskEmcalJetCDF* task
AliAnalysisTaskEmcalJetCDF *AddTaskEmcalJetCDF (
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nrho               = "",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1.,
  Double_t    jetptcutmax        = 250.,
  Double_t    jetareacut         = 0.001,
  const char *type               = "TPC",      // EMCAL, TPC
  Int_t       leadhadtype        = 0,          // AliJetContainer :: Int_t fLeadingHadronType;  0 = charged, 1 = neutral, 2 = both
  const char *taskname           = "JetCDF"
//   Int_t       nCentBins          = 1,
)
  {
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if ( !mgr ) { ::Error ( "AddTaskEmcalJetCDF", "No analysis manager to connect to." );  return NULL; }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if ( !mgr->GetInputEventHandler() ) { ::Error ( "AddTaskEmcalJetCDF", "This task requires an input event handler" ); return NULL; }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name ( taskname ); TString tracks ( ntracks );
  TString clusters ( nclusters ); TString jets ( njets );  TString rho ( nrho );

  TString acctype = type; acctype.ToUpper();
  if ( acctype.Contains ("TPC") )   { acctype = "TPC"; }
  if ( acctype.Contains ("EMCAL") ) { acctype = "EMCAL"; }
  if ( acctype.Contains ("USER") )  { acctype = "USER"; }

  if ( jetptcut < 1. ) { jetptcut = 1.; }

  TString jetstr = "jpt";
  jetstr += ( ULong_t ) ( jetptcut * 1000 );

  if ( !jets.IsNull() )    { name += "_" + jets + "_" + jetstr; }
  if ( !rho.IsNull() )     { name += "_" + rho; }
  if ( !acctype.IsNull() ) { name += "_" + acctype; }

  AliAnalysisTaskEmcalJetCDF *jetTask = new AliAnalysisTaskEmcalJetCDF ( name );
  jetTask->SetCentRange ( 0., 100. );
  jetTask->SetNCentBins ( 1 );

  AliParticleContainer *trackCont  = jetTask->AddParticleContainer ( ntracks );
  trackCont->SetClassName ( "AliVTrack" );

  AliClusterContainer *clusterCont = jetTask->AddClusterContainer ( nclusters );
//     clusterCont->SetClassName("AliVCluster");

  AliJetContainer *jetCont = jetTask->AddJetContainer ( njets, acctype, jetradius );

  if ( jetCont )
      {
      if ( !rho.IsNull() ) { jetCont->SetRhoName ( nrho ); }
      jetCont->ConnectParticleContainer ( trackCont );
      jetCont->ConnectClusterContainer ( clusterCont );
      jetCont->SetPercAreaCut ( jetareacut );
      jetCont->SetJetPtCut ( jetptcut );
      jetCont->SetJetPtCutMax ( jetptcutmax );
      jetCont->SetLeadingHadronType ( leadhadtype ); // Int_t fLeadingHadronType;  0 = charged, 1 = neutral, 2 = both
      jetCont->SetZLeadingCut ( 0., 1. );
      }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask ( jetTask );

  TString contname = name + "_histos";
  TString outfile = AliAnalysisManager::GetCommonFileName();

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer ( contname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data() );

  mgr->ConnectInput ( jetTask, 0,  cinput1 );
  mgr->ConnectOutput ( jetTask, 1, coutput1 );

  return jetTask;
  }


/// Add a AliAnalysisTaskEmcalJetCDF task - info from AliEmcalJetTask*
/// \param AliEmcalJetTask *jetFinderTask
/// \param Double_t jetptcut
/// \param Double_t jetareacut
/// \param const char *type ; either TPC, EMCAL or USER
/// \param Int_t leadhadtype ; 0 = charged, 1 = neutral, 2 = both
/// \param const char *nrho
/// \param const char *taskname
/// \return AliAnalysisTaskEmcalJetCDF* task
AliAnalysisTaskEmcalJetCDF *AddTaskEmcalJetCDF ( AliEmcalJetTask *jetFinderTask,
    Double_t     jetptcut     = 1.,
    Double_t     jetptcutmax  = 250.,
    Double_t     jetareacut   = 0.001,
    const char  *type         = "TPC",     // EMCAL, TPC
    Int_t        leadhadtype  = 0,         // AliJetContainer :: Int_t fLeadingHadronType;  0 = charged, 1 = neutral, 2 = both
    const char  *nrho         = "",
    const char  *taskname     = "JetCDF" )
  {
  if ( !jetFinderTask->InheritsFrom ( AliEmcalJetTask::Class() ) )
    { AliError("AddTaskEmcalJetSample :: task is not/ does not inherits from AliEmcalJetTask"); }

  const char *ntracks            = jetFinderTask->GetTracksName();
  const char *nclusters          = jetFinderTask->GetClusName();
  const char *njets              = jetFinderTask->GetJetsName();
  Double_t    jetradius          = jetFinderTask->GetRadius();

  return AddTaskEmcalJetCDF ( ntracks , nclusters, njets, nrho, jetradius, jetptcut, jetptcutmax, jetareacut, type, leadhadtype, taskname );
  }


/// Add a AliAnalysisTaskEmcalJetCDF task - info from char* taskname
/// \param const char* taskname ; to be retrieved from list of tasks
/// \param Double_t jetptcut
/// \param Double_t jetareacut
/// \param const char *type ; either TPC, EMCAL or USER
/// \param Int_t leadhadtype ; 0 = charged, 1 = neutral, 2 = both
/// \param const char *nrho
/// \param const char *taskname
/// \return AliAnalysisTaskEmcalJetCDF* task
AliAnalysisTaskEmcalJetCDF *AddTaskEmcalJetCDF ( const char* taskname,
    Double_t     jetptcut     = 1.,
    Double_t     jetptcutmax  = 250.,
    Double_t     jetareacut   = 0.001,
    const char  *type         = "TPC",     // EMCAL, TPC
    Int_t        leadhadtype  = 0,         // AliJetContainer :: Int_t fLeadingHadronType;  0 = charged, 1 = neutral, 2 = both
    const char  *nrho         = "",
    const char  *taskname     = "JetCDF" )
  {
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { ::Error("AddTaskEmcalJetCDF", "No analysis manager to connect to."); }

  AliEmcalJetTask* jf = dynamic_cast<AliEmcalJetTask*>(mgr->GetTask(taskname));
  if (!jf) { AliError("AddTaskEmcalJetCDF :: task is not EmcalJetTask");}

  const char *ntracks            = jf->GetTracksName();
  const char *nclusters          = jf->GetClusName();
  const char *njets              = jf->GetJetsName();
  Double_t    jetradius          = jf->GetRadius();

  return AddTaskEmcalJetCDF ( ntracks , nclusters, njets, nrho, jetradius, jetptcut, jetptcutmax, jetareacut, type, leadhadtype, taskname );
  }


// kate: indent-mode none; indent-width 2; replace-tabs on;
