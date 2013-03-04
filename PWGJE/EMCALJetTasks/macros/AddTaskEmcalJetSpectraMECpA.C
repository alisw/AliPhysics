// $Id$

AliAnalysisTaskEmcalJetSpectraMECpA* AddTaskEmcalJetSpectraMECpA(
   const char *outfilename    = "AnalysisOutput.root",
   UInt_t type                = AliAnalysisTaskEmcal::kTPC,
   const char *nRhosCh        = "rhoCh",
   const char *nRhosChEm      = "rhoChEm",
   TF1 *sfunc                 = 0,
   const Double_t radius      = 0.2,
   const Double_t minPhi      = 1.8,
   const Double_t maxPhi      = 2.74,
   const Double_t minEta      = -0.3,
   const Double_t maxEta      = 0.3,
   const char* usedTracks     = "PicoTracks",
   const char* outClusName    = "CaloClustersCorr",
   const Double_t minTrackPt  = 0.15,
   const Double_t minClusterPt = 0.30,
   const char *CentEst         = "V0A"
   )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTasEmcalJetSpectraMECpA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetSpectraMECpA", "This task requires an input event handler");
    return NULL;
  }
  
  //Run the jet finder and rho tasks first

  // Some constants for the jet finders
  const Int_t cKT                 = 0;
  const Int_t cANTIKT             = 1;
  const Int_t cFULLJETS           = 0;
  const Int_t cCHARGEDJETS        = 1;
  const Int_t cNEUTRALJETS        = 2;
    
  float AreaCut = radius*radius*TMath::Pi();

  sfunc=new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
  sfunc->SetParameter(0,0.0);
  sfunc->SetParameter(1,0.0);
  sfunc->SetParameter(2,1.5);


  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  
  AliEmcalJetTask* jetFinderTaskChBack = AddTaskEmcalJet(usedTracks,"",cKT,radius,cCHARGEDJETS,minTrackPt, minClusterPt);
  
  AliEmcalJetTask* jetFinderTaskChEmBack = AddTaskEmcalJet(usedTracks,outClusName,cKT,radius,cFULLJETS,minTrackPt, minClusterPt);

  AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(usedTracks,outClusName,cANTIKT,radius, cFULLJETS,minTrackPt,minClusterPt);

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRho.C");

  AliAnalysisTaskRho *rhochtask = AddTaskRho(jetFinderTaskChBack->GetName(),usedTracks,outClusName,nRhosCh,0.2,0,0.01,0,sfunc,2,kTRUE,nRhosCh);
  rhochtask->SetCentralityEstimator(CentEst);
  rhochtask->SetCentralityEstimator(CentEst);

  AliAnalysisTaskRho *rhochemtask = AddTaskRho(jetFinderTaskChEmBack->GetName(),usedTracks,outClusName,nRhosChEm,0.2,0,0.01,0,0,1,kTRUE,nRhosChEm);
  rhochemtask->SetCentralityEstimator(CentEst);

  const char *nJets = jetFinderTask->GetName();

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskDeltaPt.C");

  AliAnalysisTaskDeltaPt* deltapt = AddTaskDeltaPt(usedTracks,outClusName,nJets,"","","","","","",radius,1,0.557,minTrackPt,minClusterPt,AliAnalysisTaskEmcal::kTPC,"DeltaPtTask");
  deltapt->SetCentralityEstimator(CentEst);


  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------



  TString name(Form("SpectraMECpA_%s", nJets));
  AliAnalysisTaskEmcalJetSpectraMECpA *spectratask = new AliAnalysisTaskEmcalJetSpectraMECpA(name);
  spectratask->SetJetsName(jetFinderTask->GetName());
  spectratask->SetCentralityEstimator(CentEst);
  spectratask->SetAnaType(type);
  spectratask->SetRhoName(nRhosCh);
  spectratask->SetJetPhiLimits(minPhi,maxPhi);
  spectratask->SetJetEtaLimits(minEta,maxEta);
  spectratask->SetJetAreaCut(AreaCut);
  spectratask->SetTracksName(usedTracks);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(spectratask);

  // Create containers for input/output
  mgr->ConnectInput (spectratask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *cospectra = mgr->CreateContainer(name,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);
  mgr->ConnectOutput(spectratask,1,cospectra);

  return spectratask;
}
