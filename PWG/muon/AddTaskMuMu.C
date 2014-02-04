///
/// Configuration example of a task to get invariant mass spectrum of dimuons
///
/// \author: L. Aphecetche (Subatech) (laurent.aphecetche - at - subatech.in2p3.fr)
///

AliAnalysisTask* AddTaskMuMu(const char* outputname, 
                             TList* triggerClassesToConsider,
                             const char* beamYear,
                             Bool_t simulations)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuMu", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMuMu", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  cout << "inputDataType=" << inputDataType.Data() << endl;
  
  // Configure analysis
  //===========================================================================  

  if (simulations && triggerClassesToConsider )
  {
//    triggerClassesToConsider->Add(new TObjString("CMULLO-B-NOPF-MUON"));
//    triggerClassesToConsider->Add(new TObjString("CMSNGL-B-NOPF-MUON"));
    triggerClassesToConsider->Add(new TObjString("ANY"));

//    triggerClassesToConsider->Add(new TObjString("MB1"));
//    triggerClassesToConsider->Add(new TObjString("C0T0A"));

    //    triggerClassesToConsider->Add(new TObjString("MULow"));
//    triggerClassesToConsider->Add(new TObjString("V0L"));
//    triggerClassesToConsider->Add(new TObjString("V0R"));
//
// for dpmjet simulations (at least) we have the following "triggers" :
//    C0T0A,C0T0C,MB1,MBBG1,V0L,V0R,MULow,EMPTY,MBBG3,MULL,MULU,MUHigh
  }

  AliAnalysisTaskMuMu* task = new AliAnalysisTaskMuMu;
  
  AliAnalysisMuMuEventCutter* eventCutter = new AliAnalysisMuMuEventCutter(triggerClassesToConsider);
  
  AliAnalysisMuMuCutRegistry* cr = task->CutRegistry();
  
  AliAnalysisMuMuCutElement* eventTrue = cr->AddEventCut(*eventCutter,"IsTrue","const AliVEvent&","");

  AliAnalysisMuMuCutElement* ps = eventTrue;
  
  if (!simulations)
  {
    ps = cr->AddEventCut(*eventCutter,"IsPhysicsSelected","const AliInputEventHandler&","");
  }

  AliAnalysisMuMuCutElement* triggerSelection = cr->AddTriggerClassCut(*eventCutter,"SelectTriggerClass","const TString&,TString&,UInt_t,UInt_t,UInt_t","");

  cr->AddCutCombination(triggerSelection);
  
//  AliAnalysisMuMuCutElement* cutVDM = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedVDM","const AliVEvent&","");
//
//  AliAnalysisMuMuCutElement* cutTZEROPileUp = cr->AddEventCut(*eventCutter,"IsTZEROPileUp","const AliVEvent&","");
//
//  AliAnalysisMuMuCutElement* notTZEROPileUp = cr->Not(*cutTZEROPileUp);

  cr->AddCutCombination(ps);

  task->SetBeamYear(beamYear);

  AliAnalysisMuMuGlobal* globalAnalysis = 0x0;//new AliAnalysisMuMuGlobal;
  
  AliAnalysisMuMuSingle* singleAnalysis = new AliAnalysisMuMuSingle;
  
  AliAnalysisMuMuMinv* minvAnalysis = new AliAnalysisMuMuMinv;

//  TFile f("$HOME/Downloads/ResponseMatrix_QASPDZPSALL_ANY.root");
//  TH2* h = static_cast<TH2*>(f.Get("ResponseMatrix"));
  
  AliAnalysisMuMuNch* nchAnalysis = new AliAnalysisMuMuNch; // (0x0,-2,2,-40,40);
  
//  AliAnalysisMuMuNch(TH2* spdCorrection=0x0,
//                     Int_t nbinsEta=10, Double_t etaMin=-0.5, Double_t etaMax=0.5,
//                     Int_t nbinsZ=320, Double_t zmin=-40, Double_t zmax=40);


  if ( globalAnalysis )
  {
    AliAnalysisMuMuCutElement* triggerAll = cr->AddTriggerClassCut(*globalAnalysis,"SelectAnyTriggerClass","const TString&,TString&","");
    
    cr->AddCutCombination(triggerAll);
    
    task->AdoptSubAnalysis(globalAnalysis);
  }
  
  if ( nchAnalysis )
  {
//    AliAnalysisMuMuCutElement* trketa2 = cr->AddEventCut(*nchAnalysis,"HasAtLeastNTrackletsInEtaRange","const AliVEvent&,Int_t, Double_t& , Double_t& ","1,-2.0,2.0");
//    AliAnalysisMuMuCutElement* trketa1 = cr->AddEventCut(*nchAnalysis,"HasAtLeastNTrackletsInEtaRange","const AliVEvent&,Int_t, Double_t& , Double_t& ","1,-1.0,1.0");
//    AliAnalysisMuMuCutElement* trk5eta1 = cr->AddEventCut(*nchAnalysis,"HasAtLeastNTrackletsInEtaRange","const AliVEvent&,Int_t, Double_t& , Double_t& ","5,-1.0,1.0");
//    
//    cr->AddCutCombination(ps,trketa1);
//    cr->AddCutCombination(ps,trk5eta1);
//    cr->AddCutCombination(ps,trketa2);
    
    AliAnalysisMuMuCutElement* cutZ18= cr->AddEventCut(*eventCutter,"IsAbsZBelowValue","const AliVEvent&,Double_t","18");
    
    AliAnalysisMuMuCutElement* cutZ10 = cr->AddEventCut(*eventCutter,"IsAbsZBelowValue","const AliVEvent&,Double_t","10");
    
    cr->AddCutCombination(ps,cutZ18);
    cr->AddCutCombination(ps,cutZ10);
    
    task->AdoptSubAnalysis(nchAnalysis);
  }
  
  if ( singleAnalysis )
  {
    AliAnalysisMuMuCutElement* trackTrue = cr->AddTrackCut(*cr,"AlwaysTrue","const AliVParticle&","");
    AliAnalysisMuMuCutElement* rabs = cr->AddTrackCut(*singleAnalysis,"IsRabsOK","const AliVParticle&","");
    AliAnalysisMuMuCutElement* eta = cr->AddTrackCut(*singleAnalysis,"IsEtaInRange","const AliVParticle&,Double_t&,Double_t&","-4.0,-2.5");
    
    cr->AddCutCombination(trackTrue);

    cr->AddCutCombination(rabs);
    
    cr->AddCutCombination(eta);

    cr->AddCutCombination(rabs,eta);

    task->AdoptSubAnalysis(singleAnalysis);

    if ( minvAnalysis )
    {
      AliAnalysisMuMuCutElement* pairTrue = cr->AddTrackPairCut(*cr,"AlwaysTrue","const AliVParticle&, const AliVParticle&","");
      
      AliAnalysisMuMuCutElement* pairy = cr->AddTrackPairCut(*minvAnalysis,"IsRapidityInRange","const AliVParticle&,const AliVParticle&,Double_t&,Double_t&","-4,-2.5");
  
      cr->AddCutCombination(rabs,eta,pairy);

      cr->AddCutCombination(rabs,pairy);

      cr->AddCutCombination(pairy);

      cr->AddCutCombination(pairTrue);

      task->AdoptSubAnalysis(minvAnalysis);
    }
  }
  
  /// below are the kind of configurations that can be performed :
  /// - adding cuts (at event, track or pair level)
  /// - adding bins (in pt, y, centrality, etc...) for minv (and meanpt)

  AliAnalysisMuMuBinning* binning = task->Binning();
  
  if (minvAnalysis)
  {
    binning->AddBin("psi","pt");
    
    // pt binning
    
    binning->AddBin("psi","pt", 0.0, 1.0,"IGOR");
    binning->AddBin("psi","pt", 1.0, 2.0,"IGOR");
    binning->AddBin("psi","pt", 2.0, 3.0,"IGOR");
    binning->AddBin("psi","pt", 3.0, 4.0,"IGOR");
    binning->AddBin("psi","pt", 4.0, 5.0,"IGOR");
    binning->AddBin("psi","pt", 5.0, 6.0,"IGOR");
    binning->AddBin("psi","pt", 6.0, 7.0,"IGOR");
    binning->AddBin("psi","pt", 7.0, 9.0,"IGOR");
    binning->AddBin("psi","pt", 9.0,11.0,"IGOR");
    binning->AddBin("psi","pt",11.0,15.0,"IGOR");
    binning->AddBin("psi","pt",15.0,25.0,"IGOR");
    
    // y binning
    
    binning->AddBin("psi","y",-4,-2.5,"ILAB");
    
    for ( Int_t i = 0; i < 6; ++i )
    {
      Double_t y = -4+i*0.25;
      
      binning->AddBin("psi","y",y,y+0.25,"6PACK");
    }
    
    // nch binning
    if (nchAnalysis)
    {
      binning->AddBin("psi","nch",0.0,30.0,"JAVI"); // 0 - 29
      binning->AddBin("psi","nch",29.0,42.0,"JAVI"); // 30 - 41
      binning->AddBin("psi","nch",41.0,56.0,"JAVI"); // 42 - 55
      binning->AddBin("psi","nch",55.0,72.0,"JAVI"); // 56 - 72
      binning->AddBin("psi","nch",72.0,301.0,"JAVI"); // 73 - 300
    }
  }
  
  // v0 centrality binning
  
  binning->AddBin("centrality","v0a");
//  binning->AddBin("centrality","v0a",0,5);
//  binning->AddBin("centrality","v0a",5,10);
//  binning->AddBin("centrality","v0a",10,20);
//  binning->AddBin("centrality","v0a",20,40);
//  binning->AddBin("centrality","v0a",40,60);
//  binning->AddBin("centrality","v0a",60,80);
//  binning->AddBin("centrality","v0a",80,100);

  // disable some histograms if we don't want them
  task->DisableHistograms("^V02D");
  task->DisableHistograms("^dca");
  task->DisableHistograms("^Chi12");
  task->DisableHistograms("^Rabs12");

  // add the configured task to the analysis manager
  mgr->AddTask(task);  
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();  

  AliAnalysisDataContainer *coutputHC = 
  mgr->CreateContainer("OC",AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,outputname);

  AliAnalysisDataContainer *coutputCC = 
  mgr->CreateContainer("CC",AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer,outputname);
  
  AliAnalysisDataContainer* cparam = 
  mgr->CreateContainer("BIN", AliAnalysisMuMuBinning::Class(),AliAnalysisManager::kParamContainer,outputname);
  
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputHC);
  mgr->ConnectOutput(task, 2, coutputCC);
  mgr->ConnectOutput(task, 3, cparam);
  
  return task;
}

