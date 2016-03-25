/// \file AddTaskEmcalJetCDF.C
/// \brief Adds a AliAnalysisTaskEmcalJetCDF analysis task and coresponding containers
///
/// Analysis of jets (all and leading) distribution of pt and multiplicity, R distribution, N{70,75,80,85,90} and Pt{70,75,80,85,90}
///
/// \author Adrian SEVCENCO <Adrian.Sevcenco@cern.ch>, Institute of Space Science, Romania
/// \date Mar 19, 2015

/// Add a AliAnalysisTaskEmcalJetCDF task - detailed signature
/// \param const char* ntracks : name of tracks collection
/// \param const char* nclusters : name of clusters collection
/// \param const char* taskname
/// \param Double_t track pt cut
/// \param Double_t cluster E cut
/// \param Bool_t debug flag
/// \return AliAnalysisTaskEmcalJetCDF* task
AliAnalysisTaskEmcalJetCDF* AddTaskEmcalJetCDF (
  const char* ntracks                      = "usedefault",
  const char* nclusters                    = "usedefault",
  const char* taskname                     = "CDF",
  Double_t    trackPtCut                   = 0.15,
  Double_t    clusECut                     = 0.30,
  Bool_t      debug                        = kFALSE
)
  {
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if ( !mgr ) { ::Error ( "AddTaskEmcalJetCDF", "No analysis manager to connect to." );  return NULL; }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) { ::Error ( "AddTaskEmcalJetCDF", "This task requires an input event handler" ); return NULL; }

  enum EDataType_t { kUnknown, kESD, kAOD }; EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) { dataType = kESD; }
  else
  if (handler->InheritsFrom("AliAODInputHandler")) { dataType = kAOD; }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name ( taskname );
  TString tracks ( ntracks );
  TString clusters ( nclusters );

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

  if (debug) { std::cout << "AliAnalysisTaskEmcalJetCDF :: Name of CDF task is : " << name.Data() << std::endl;}

  AliAnalysisTaskEmcalJetCDF* cdfTask = new AliAnalysisTaskEmcalJetCDF ( name.Data() );
  cdfTask->SetVzRange(-10,10);
  cdfTask->SetNeedEmcalGeom(kFALSE);

  if ( tracks.EqualTo("mcparticles") )
    {
    AliMCParticleContainer* mcpartCont = cdfTask->AddMCParticleContainer ( tracks.Data() );
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
    }
  else
  if ( tracks.EqualTo("tracks") || tracks.EqualTo("Tracks") )
    {
    AliTrackContainer* trackCont = cdfTask->AddTrackContainer( tracks.Data() );
    }
  else
  if ( !tracks.IsNull())
    { cdfTask->AddParticleContainer(tracks.Data()); }

  AliParticleContainer* partCont  = cdfTask->GetParticleContainer(0);
  if (partCont) { partCont->SetParticlePtCut(trackPtCut); }

  AliClusterContainer* clusterCont = cdfTask->AddClusterContainer( clusters.Data() );
  if (clusterCont)
    {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(clusECut);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask ( cdfTask );

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();

  TString contname = name + "_histos";
  TString outfile ("CDFhistos.root");
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer ( contname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data() );

  mgr->ConnectInput  ( cdfTask, 0,  cinput1 );
  mgr->ConnectOutput ( cdfTask, 1, coutput1 );

  return cdfTask;
  }


// kate: indent-mode none; indent-width 2; replace-tabs on;
