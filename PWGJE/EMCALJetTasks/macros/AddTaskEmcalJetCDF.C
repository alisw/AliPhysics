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
AliAnalysisTaskEmcalJetCDF* AddTaskEmcalJetCDF (
  AliEmcalJetTask* jf_task,
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char *njets              = "Jets",
  const char *nrho               = "",
  Double_t    jetptcut           = 1.,
  Double_t    jetptcutmax        = 250.,
  Double_t    jetareacut         = 0.001,
  AliJetContainer::JetAcceptanceType acceptance = AliJetContainer::kTPCfid,
  Int_t       leadhadtype        = 0,          // AliJetContainer :: Int_t fLeadingHadronType;  0 = charged, 1 = neutral, 2 = both
  const char *taskname           = "CDF"
)
  {
//   const char* ncells             = "usedefault",

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if ( !mgr ) { ::Error ( "AddTaskEmcalJetCDF", "No analysis manager to connect to." );  return NULL; }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) { ::Error ( "AddTaskEmcalJetCDF", "This task requires an input event handler" ); return NULL; }

  if ( (!jf_task) && (!jf_task->InheritsFrom ( "AliEmcalJetTask")) ) { ::Error ( "AddTaskEmcalJetCDF", "Invalid AliEmcalJetTask* pointer"); return NULL; }

  enum EDataType_t { kUnknown, kESD, kAOD };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) { dataType = kESD; }
  else
  if (handler->InheritsFrom("AliAODInputHandler")) { dataType = kAOD; }

  AliJetContainer::EJetType_t       jettype_t = jf_task->GetJetType();
  AliJetContainer::ERecoScheme_t     recomb_t = jf_task->GetRecombScheme();
  AliJetContainer::EJetAlgo_t       jetalgo_t = AliJetContainer::antikt_algorithm;
  Double_t                             radius = jf_task->GetRadius();

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name ( taskname );
  TString tracks ( ntracks );
  TString clusters ( nclusters );
//TString cellName(ncells);
  TString jets ( njets );
  TString rho ( nrho );

  if ( tracks.EqualTo("usedefault") )
    {
    if ( dataType == kESD ) { tracks = "Tracks"; }
    else
    if ( dataType == kAOD ) { tracks = "tracks"; }
    else
      { tracks = ""; }
    }

  if ( clusters.EqualTo("usedefault") )
    {
    if ( dataType == kESD ) { clusters = "CaloClusters"; }
    else
    if ( dataType == kAOD ) { clusters = "caloClusters"; }
    else
      { clusters = ""; }
    }

//   if (cellName.EqualTo ("usedefault"))
//     {
//     if (dataType == kESD) { cellName = "EMCALCells"; }
//     else cellName = "emcalCells"; }
//     }

  TString acctype = "acc"; acctype += TString::Itoa(acceptance,10);
  if ( jetptcut < 1. ) { jetptcut = 1.; }

  TString jetstr = "jpt";
  jetstr += ( ULong_t ) ( jetptcut * 1000 );

  if ( !jets.IsNull() )    { name += "_" + jets + "_" + jetstr; }
  if ( !rho.IsNull() )     { name += "_" + rho; }
  if ( !acctype.IsNull() ) { name += "_" + acctype; }

  AliAnalysisTaskEmcalJetCDF* jetTask = new AliAnalysisTaskEmcalJetCDF ( name );
  jetTask->SetVzRange(-10,10);
  jetTask->SetNeedEmcalGeom(kFALSE);

  if ( tracks.EqualTo("mcparticles") )
    {
    AliMCParticleContainer* mcpartCont = jetTask->AddMCParticleContainer ( tracks.Data() );
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
    }
  else
  if ( tracks.EqualTo("tracks") || tracks.EqualTo("Tracks") ) {
    AliTrackContainer* trackCont = jetTask->AddTrackContainer( tracks.Data() );
    trackCont->SetFilterHybridTracks(kTRUE);
  }
  else
  if ( !tracks.IsNull())
    {
    jetTask->AddParticleContainer(tracks.Data());
    }

  AliParticleContainer* partCont  = jetTask->GetParticleContainer(0);
  if (partCont) { partCont->SetParticlePtCut(0.15); }

  AliClusterContainer* clusterCont = NULL;
  if ( !clusters.IsNull() )
    {
    clusterCont = jetTask->AddClusterContainer ( clusters.Data() );
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(0.30);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    }

  AliJetContainer *jetCont = jetTask->AddJetContainer(
                              jettype_t,
                              jetalgo_t,
                              recomb_t,
                              radius,
                              acceptance,
                              partCont,
                              clusterCont,
                              name);

  if ( !jetCont ) { std::cout << "AddTaskEmcalJetCDF.C :: could not add jetCont" << std::endl; return NULL; }

  if ( !rho.IsNull() ) { jetCont->SetRhoName ( nrho ); }
  jetCont->SetPercAreaCut ( jetareacut );
  jetCont->SetJetPtCut ( jetptcut );
  jetCont->SetJetPtCutMax ( jetptcutmax );
  jetCont->SetLeadingHadronType ( leadhadtype ); // Int_t fLeadingHadronType;  0 = charged, 1 = neutral, 2 = both
  jetCont->SetMaxTrackPt(1000);
  jetCont->SetNLeadingJets(1);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask ( jetTask );

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();

  TString contname = name + "_histos";
  TString outfile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer ( contname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data() );

  mgr->ConnectInput  ( jetTask, 0,  cinput1 );
  mgr->ConnectOutput ( jetTask, 1, coutput1 );

  return jetTask;
  }

AliAnalysisTaskEmcalJetCDF* AddTaskEmcalJetCDF (
  const char* taskname,
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char *njets              = "Jets",
  const char *nrho               = "",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1.,
  Double_t    jetptcutmax        = 250.,
  Double_t    jetareacut         = 0.001,
  AliJetContainer::JetAcceptanceType type = AliJetContainer::kTPCfid,
  Int_t       leadhadtype        = 0,          // AliJetContainer :: Int_t fLeadingHadronType;  0 = charged, 1 = neutral, 2 = both
  const char *taskname           = "CDF"
)
  {
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { ::Error("AddTaskEmcalJetCDF", "No analysis manager to connect to."); }

  AliEmcalJetTask* jf = dynamic_cast<AliEmcalJetTask*>(mgr->GetTask(taskname));
  if (!jf) { AliError("AddTaskEmcalJetCDF :: task is not EmcalJetTask");}

  return AddTaskEmcalJetCDF ( jf, ntracks, nclusters, njets, nrho, jetptcut, jetptcutmax, jetareacut, type, leadhadtype, taskname);
  }

// kate: indent-mode none; indent-width 2; replace-tabs on;
