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
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char *njets              = "Jets",
  const char *nrho               = "",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1.,
  Double_t    jetptcutmax        = 250.,
  Double_t    jetareacut         = 0.001,
  const char *type               = "TPC",      // EMCAL, TPC
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

  enum EDataType_t { kUnknown, kESD, kAOD };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) { dataType = kESD; }
  else
  if (handler->InheritsFrom("AliAODInputHandler")) { dataType = kAOD; }


  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name ( taskname );
  TString tracks ( ntracks );
//TString cellName(ncells);
  TString clusters ( nclusters );
  TString jets ( njets );
  TString rho ( nrho );

  if (tracks.EqualTo("usedefault"))
    {
    if (dataType == kESD) { tracks = "Tracks"; }
    else { tracks = "tracks"; }
    }

  if (clusters.EqualTo("usedefault"))
    {
    if (dataType == kESD) { clusters = "CaloClusters"; }
    else { clusters = "caloClusters"; }
    }

//   if (cellName.EqualTo ("usedefault"))
//     {
//     if (dataType == kESD) { cellName = "EMCALCells"; }
//     else cellName = "emcalCells"; }
//     }

  TString acctype = type;
  if ( jetptcut < 1. ) { jetptcut = 1.; }

  TString jetstr = "jpt";
  jetstr += ( ULong_t ) ( jetptcut * 1000 );

  if ( !jets.IsNull() )    { name += "_" + jets + "_" + jetstr; }
  if ( !rho.IsNull() )     { name += "_" + rho; }
  if ( !acctype.IsNull() ) { name += "_" + acctype; }

  AliAnalysisTaskEmcalJetCDF* jetTask = new AliAnalysisTaskEmcalJetCDF ( name );
  jetTask->SetVzRange(-10,10);
  jetTask->SetNeedEmcalGeom(kFALSE);

  AliParticleContainer *trackCont  = jetTask->AddParticleContainer ( tracks.Data() );
  AliClusterContainer *clusterCont = jetTask->AddClusterContainer ( clusters.Data() );

  AliJetContainer *jetCont = jetTask->AddJetContainer ( njets, acctype.Data(), jetradius );
  if ( jetCont )
      {
      if ( !rho.IsNull() ) { jetCont->SetRhoName ( nrho ); }
      jetCont->ConnectParticleContainer ( trackCont );
      jetCont->ConnectClusterContainer ( clusterCont );
      jetCont->SetPercAreaCut ( jetareacut );
      jetCont->SetJetPtCut ( jetptcut );
      jetCont->SetJetPtCutMax ( jetptcutmax );
      jetCont->SetLeadingHadronType ( leadhadtype ); // Int_t fLeadingHadronType;  0 = charged, 1 = neutral, 2 = both
      jetCont->SetMaxTrackPt(1000);
      }

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

// kate: indent-mode none; indent-width 2; replace-tabs on;
