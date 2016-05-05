// $Id$

AliAnalysisTaskEmcalJetSpectraMECpA* AddTaskEmcalJetSpectraMECpA(
   const char *outfilename    = "AnalysisOutput.root",
   const char *nRhosCh        = "rhoCh",
   const char *nRhosChEm      = "rhoChEm",
   const char *nRhosEm        = "rhoEm",
   const Double_t scale       = 1.0,
   const Double_t radius      = 0.2,
   const Double_t minPhi      = 1.8,
   const Double_t maxPhi      = 2.74,
   const Double_t minEta      = -0.3,
   const Double_t maxEta      = 0.3,
   const char* usedTracks     = "PicoTracks",
   const char* outClusName    = "caloClustersCorr",
   const Double_t minTrackPt  = 0.15,
   const Double_t minClusterPt = 0.30,
   const char *CentEst         = "V0A",
   const Int_t type            =1
   )
{  

  //cout << "used tracks: " << usedTracks <<endl;

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

  //UInt_t typeTPC                = 0;
  //UInt_t typeEMC                = AliAnalysisTaskEmcal::kEMCAL;
  //UInt_t typeTPC                = 0;
  //UInt_t typeEMC                = 1;

  char *typeTPC                = "TPC";
  char *typeEMC                = "EMC";

    
  float AreaCut = 0.6*radius*radius*TMath::Pi();

  TF1 *sfunc=new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
  sfunc->SetParameter(0,0.0);
  sfunc->SetParameter(1,0.0);
  sfunc->SetParameter(2,scale);


  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");

  //const char *nJets;
  TString nJets("");
  TString nJetsCh("");

  TString scaledname(Form("%s_Scaled", nRhosCh));
  TString newrhoname(Form("%s_All", nRhosCh));
  //TString scaledname(Form("%s_Scaled", newrhoname));


  if(!(usedTracks=="")){
    //cout << "USEDTRACKS EXISTS" << usedTracks <<endl;
    AliEmcalJetTask* jetFinderTaskChBack = AddTaskEmcalJet(AliEmcalJetTask::kAKT | AliEmcalJetTask::kChargedJet, usedTracks,"",minTrackPt, minClusterPt, 0.005, radius, 1, "Jet", 0.1);
  
    AliEmcalJetTask* jetFinderTaskChBackall = AddTaskEmcalJet(AliEmcalJetTask::kAKT | AliEmcalJetTask::kChargedJet, usedTracks,"",minTrackPt, minClusterPt, 0.005, radius, 1, "Jets_allpt", 0.0);
    //jetFinderTaskChBackall->SetMinJetPt(0);

    AliEmcalJetTask* jetFinderTaskChSig = AddTaskEmcalJet(AliEmcalJetTask::kAKT | AliEmcalJetTask::kChargedJet, usedTracks,"",minTrackPt, minClusterPt, 0.005, radius, 1, "Jet", 0.1);

    AliEmcalJetTask* jetFinderTaskChEmBack = AddTaskEmcalJet(AliEmcalJetTask::kKT | AliEmcalJetTask::kChargedJet, usedTracks,outClusName,minTrackPt, minClusterPt, 0.005, radius, 1, "Jet", 0.0);
    //jetFinderTaskChEmBack->SetMinJetPt(0);

  AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(AliEmcalJetTask::kAKT | AliEmcalJetTask::kFullJet, usedTracks,outClusName,minTrackPt,minClusterPt, 0.005, radius, 1, "Jet", 0.1);

  AliAnalysisTaskRhoSparse *rhochtask = AddTaskRhoSparse(jetFinderTaskChBack->GetName(),jetFinderTaskChSig->GetName(),usedTracks,outClusName,nRhosCh,radius,typeTPC,0.01,0,0,sfunc,0,kTRUE,nRhosCh);
  rhochtask->SetCentralityEstimator(CentEst);

  AliAnalysisTaskRhoSparse *rhochalltask = AddTaskRhoSparse(jetFinderTaskChBackall->GetName(),jetFinderTaskChSig->GetName(),usedTracks,outClusName,newrhoname,radius,0,0.0,0,0,sfunc,0,kTRUE,newrhoname);
  rhochtask->SetCentralityEstimator(CentEst);
						   

  AliAnalysisTaskRhoSparse *rhochemtask = AddTaskRhoSparse(jetFinderTaskChEmBack->GetName(),jetFinderTask->GetName(),usedTracks,outClusName,nRhosChEm,radius,typeEMC,0.01,0,0,0,0,kTRUE,nRhosChEm);
  rhochemtask->SetCentralityEstimator(CentEst);

  //nJets=jetFinderTask->GetName();
  nJets+=jetFinderTask->GetName();
  nJetsCh+=jetFinderTaskChSig->GetName();

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskDeltaPt.C");

  TString deltaname(Form("DeltaPt_%s_Scaled", nRhosCh));
  AliAnalysisTaskDeltaPt* deltapt = AddTaskDeltaPt(usedTracks,outClusName,nJets,"","","","","",scaledname,radius,AreaCut,minTrackPt,minClusterPt,typeEMC,deltaname);
  deltapt->SetCentralityEstimator(CentEst);

  TString chdeltaname(Form("DeltaPt_%s", nRhosCh));
  AliAnalysisTaskDeltaPt* deltaptch = AddTaskDeltaPt(usedTracks,"",nJetsCh,"","","","","",nRhosCh,radius,AreaCut,minTrackPt,minClusterPt,typeTPC,chdeltaname);
  deltaptch->SetCentralityEstimator(CentEst);

  TString emcdeltaname(Form("DeltaPt_%s", nRhosChEm));
  AliAnalysisTaskDeltaPt* deltaptEMC = AddTaskDeltaPt(usedTracks,outClusName,nJets,"","","","","",nRhosChEm,radius,AreaCut,minTrackPt,minClusterPt,typeEMC,emcdeltaname);
  deltaptEMC->SetCentralityEstimator(CentEst);

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskScale.C");


  Int_t radlabel=(Int_t)floor(radius*100+0.5);
  Int_t mincluslabel=(Int_t)floor(minClusterPt*1000+0.5);
  Int_t mintracklabel=(Int_t)floor(minTrackPt*1000+0.5);
  TString scalename(Form("Scale_R0%d", radlabel));

  AliAnalysisTaskScale* scaletask = AddTaskScale(usedTracks,outClusName,minTrackPt,minClusterPt,scalename);
  scaletask->SetCentralityEstimator(CentEst);
  scaletask->SetScaleFunction(sfunc);

  
  }

  //cout << "Running non charged jet finders..." <<endl;
  AliEmcalJetTask* jetFinderTaskEm = AddTaskEmcalJet(AliEmcalJetTask::kAKT | AliEmcalJetTask::kChargedJet,"",outClusName,minTrackPt,minClusterPt, 0.005, radius, 1, "Jet", 0.1);

  AliEmcalJetTask* jetFinderTaskEmBack = AddTaskEmcalJet(AliEmcalJetTask::kKT | AliEmcalJetTask::kNeutralJet,"",outClusName,minTrackPt, minClusterPt, 0.005, radius, 1, "Jet", 0.1);

  //cout << "Running non charged rho task..." <<endl;

  AliAnalysisTaskRhoSparse *rhoemtask = AddTaskRhoSparse(jetFinderTaskEmBack->GetName(),jetFinderTaskEm->GetName(),usedTracks,outClusName,nRhosEm,radius,0,0.01,0,0,0,1,kTRUE,nRhosEm);
  rhoemtask->SetCentralityEstimator(CentEst);

  if(usedTracks=="") nJets +=jetFinderTaskEm->GetName();





  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------


  //cout << "Ready to run my task..." << nJets <<endl;
  
  TString name(Form("SpectraMECpA_%s", nJets.Data()));
  AliAnalysisTaskEmcalJetSpectraMECpA *spectratask = new AliAnalysisTaskEmcalJetSpectraMECpA(name);
  spectratask->AddJetContainer(nJets.Data());
  spectratask->SetCentralityEstimator(CentEst);
  spectratask->SetVzRange(-10,10);

  if(type==0){
    //spectratask->SetAnaType(typeTPC);
    spectratask->SetAnaType(0);
    spectratask->SetRhoName(nRhosCh);
  }else{
    //spectratask->SetAnaType(typeEMC);
    spectratask->SetAnaType(1);
    if(!(usedTracks=="")) spectratask->SetRhoName(scaledname);
    else spectratask->SetRhoName(nRhosEm);
  }
  spectratask->SetJetPhiLimits(minPhi,maxPhi);
  spectratask->SetJetEtaLimits(minEta,maxEta);
  spectratask->SetJetAreaCut(AreaCut);
  spectratask->AddParticleContainer(usedTracks);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------


  mgr->AddTask(spectratask);

  // Create containers for input/output
  mgr->ConnectInput (spectratask, 0, mgr->GetCommonInputContainer());
  AliAnalysisDataContainer *cospectra = mgr->CreateContainer(name,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);
  mgr->ConnectOutput(spectratask,1,cospectra);


							 
  return spectratask;
}
