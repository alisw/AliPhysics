///
/// Configuration example of a task to correct charged particle multiplicity
///

AliAnalysisTaskMuMu* AddTaskMuMuDT(const char* outputname,
				 const char* beamYear = "PP",
				 Bool_t simulations = kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager(); // Get the manager
  if (!mgr) {
    ::Error("AddTaskMuMu", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMuMu", "This task requires an input event handler");
    return NULL;
  }

  // =========== list of trigger =========

  // --- Your trigger list, could be what you want ---
  TList* triggers = new TList;
  triggers->SetOwner(kTRUE);
  triggers->Add(new TObjString("CMUL7-B-NOPF-MUFAST"));
  triggers->Add(new TObjString("CINT7-B-NOPF-MUFAST"));


  // ========= Configure inputmaps (Default on is in AliMuonEventCuts). Do not touch it for pp@13TeV Analysis ========

  TList* triggerInputsMap = new TList();
  triggerInputsMap->SetOwner(kTRUE);
  triggerInputsMap->Add(new TObjString("0MSL:17"));
  triggerInputsMap->Add(new TObjString("0MSH:18"));
  triggerInputsMap->Add(new TObjString("0MLL:19"));
  triggerInputsMap->Add(new TObjString("0MUL:20"));
  triggerInputsMap->Add(new TObjString("0V0M:3"));
  triggerInputsMap->Add(new TObjString("0TVX:2"));

  // =========

  // ========= Create task =========

  // File to plug after the first run of the task.
  // You need to give it the MeanTrackletsVsZVertex TProfile to use the data-driven multiplicity correction method

  /*  TProfile* hT =0x0;
  
  const char * filename1 ="TrackletsCorr.root";
  AliAnalysisMuMu ana1(filename1,"","","mumu.pp2015.config");
  //=====-10-10========================================
  hT = static_cast<TProfile*>(ana1.OC()->GetObject("/PSMULHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00/CMUL7-B-NOPF-MUFAST/PP/MeanTrackletsVsZVertex")->Clone()); 
  */ 
  

  
  // cut to fill histo from the task
  Double_t etaMin = -1.;
  Double_t etaMax = 1.;

  Double_t zMin = -10.;
  Double_t zMax = 10.;


  //General task that contains all the sub-task analysis
  AliAnalysisTaskMuMu* task            = new AliAnalysisTaskMuMu;

  // Raw Tracklets sub-analysis. The first configuration is meant to produce multiplicity QA plots and stuff (First run)
  // For the second run, you can use the second conf. providing the MeanTrackletsVsZVertex from the first run output file.

  AliAnalysisMuMuNch* nch              = new AliAnalysisMuMuNch(0x0,0x0,-1.,etaMin,etaMax,zMin,zMax,kFALSE);//

  //==========Data Driven method===================
  //AliAnalysisMuMuNch* nch              = new AliAnalysisMuMuNch(0x0,hT,22.0908,etaMin,etaMax,zMin,zMax,kFALSE);
 
  
  // Inv. Mass sub-analysis. For the first run, you can leave it to NULL, and only plug it for the second run.
  AliAnalysisMuMuMinv* minvAnalysis    =  0x0; // Call the task
  //AliAnalysisMuMuMinv* minvAnalysis    =  new AliAnalysisMuMuMinv; // Call the task

  // Single muon sub-analysis. Provides the cut methods for single muons, never unplug
  AliAnalysisMuMuSingle* singleAnalysis= new AliAnalysisMuMuSingle;
  task->SetBeamYear(beamYear);


  // ========= Configure cut on event =========

  // --- standard Analysis ---

  AliAnalysisMuMuCutRegistry* cr                  = task->CutRegistry(); // Set CutRegistry

  // Class where all the event cut methods are stored
  AliAnalysisMuMuEventCutter* eventCutter         = new AliAnalysisMuMuEventCutter(triggers,triggerInputsMap);

  // Call several event cuts
  
  AliAnalysisMuMuCutElement* eventTrue           = cr->AddEventCut(*eventCutter,"IsTrue","const AliVEvent&",""); // ALL Events
  AliAnalysisMuMuCutElement* cutSPDZ10           = cr->AddEventCut(*eventCutter,"IsAbsZSPDBelowValue","const AliVEvent&,Double_t","10");
  AliAnalysisMuMuCutElement* triggerSelection    = cr->AddTriggerClassCut(*eventCutter,"SelectTriggerClass","const TString&,TString&,UInt_t,UInt_t,UInt_t","");
  AliAnalysisMuMuCutElement* cutHasSPD           = cr->AddEventCut(*eventCutter,"HasSPDVertex","const AliVEvent&","");
  AliAnalysisMuMuCutElement* cutZQA              = cr->AddEventCut(*eventCutter,"IsSPDzQA","const AliVEvent&,Double_t,Double_t","0.25,0.5");
  AliAnalysisMuMuCutElement* cutTZEROPileUp      = cr->AddEventCut(*eventCutter,"IsTZEROPileUp","const AliVEvent&","");
  AliAnalysisMuMuCutElement* SPDPileUp           = cr->AddEventCut(*eventCutter,"IsSPDPileUp","AliVEvent&","");
  AliAnalysisMuMuCutElement* NoSPDPileUp         = cr->Not(*SPDPileUp);
 
 

  // Here you can change the physic selction as you want, see AliAnalysisMuMuEventCutter.h. It is similar to SetCollisionCandidates()
  AliAnalysisMuMuCutElement * ps1                 = eventTrue;
  //AliAnalysisMuMuCutElement * ps2                 = eventTrue;
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(kFALSE,kTRUE,0,kFALSE); 
 

  if (!simulations)
    {
      ps1 = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedMUL","const AliInputEventHandler&","");
    } 

  // Add a combination to the cut registry. This is your actual event cuts. You can use several combination at the same time
  cr->AddCutCombination(ps1,cutHasSPD,cutZQA,cutSPDZ10,triggerSelection);
  //cr->AddCutCombination(ps2,cutHasSPD,cutZQA,cutSPDZ10,triggerSelection,NoSPDPileUp);

  // cr->AddCutCombination(ps1,cutHasSPD,cutZQA,triggerSelection);


  // ========= Configure sub-analysis =========

  if ( singleAnalysis ){

    // --- Create combination of cuts to apply on single muons ---
    AliAnalysisMuMuCutElement* trackTrue = cr->AddTrackCut(*cr,"AlwaysTrue","const AliVParticle&",""); // Apply "AlwaysTrue" cut on AliVParticle derived from AliAnalysisMuMuSingle
    AliAnalysisMuMuCutElement* rabs      = cr->AddTrackCut(*singleAnalysis,"IsRabsOK","const AliVParticle&","");
    AliAnalysisMuMuCutElement* matchlow  = cr->AddTrackCut(*singleAnalysis,"IsMatchingTriggerLowPt","const AliVParticle&","");
    AliAnalysisMuMuCutElement* eta       = cr->AddTrackCut(*singleAnalysis,"IsEtaInRange","const AliVParticle&","");
    cr->AddCutCombination(trackTrue,rabs,matchlow,eta);

    // --- Disable some histo. Comment the one you want ---
    // singleAnalysis->DisableHistograms("BCX");
    // singleAnalysis->DisableHistograms("EtaRapidityMu*");
    // singleAnalysis->DisableHistograms("PtEtaMu*");
    // singleAnalysis->DisableHistograms("PtRapidityMu*");
    // singleAnalysis->DisableHistograms("PEtaMu*");
    // singleAnalysis->DisableHistograms("PtPhiMu*");
    // singleAnalysis->DisableHistograms("Chi2Mu*");
    // singleAnalysis->DisableHistograms("dcaP23Mu*");
    // singleAnalysis->DisableHistograms("dcaPwPtCut23Mu*");
    // singleAnalysis->DisableHistograms("dcaP310Mu*");
    // singleAnalysis->DisableHis∏git stattograms("dcaPwPtCut310Mu*");

    // --- separate positive and negative muons ---
    // singleAnalysis->ShouldSeparatePlusAndMinus(kTRUE);

    // Adding the sub analysis
    task->AdoptSubAnalysis(singleAnalysis);

    if ( minvAnalysis ){

      // --- Array of cut elements ---
      TObjArray cutElements;
      TObjArray cutElementsmix;

      // --- Cuts on muon pairs level ---
      AliAnalysisMuMuCutElement* pairTrue = cr->AddTrackPairCut(*cr,"AlwaysTrue","const AliVParticle&,const AliVParticle&","");// Apply "AlwaysTrue" cut on AliVParticle derived from AliAnalysisMuMuMinv
      AliAnalysisMuMuCutElement* pairy    = cr->AddTrackPairCut(*minvAnalysis,"IsRapidityInRange","const AliVParticle&,const AliVParticle&","");
      AliAnalysisMuMuCutElement* ptRange  = cr->AddTrackPairCut(*minvAnalysis,"IsPtInRange","const AliVParticle&,const AliVParticle&,Double_t&,Double_t&","0.,12.");
      cutElements.Add(pairTrue);
      cutElements.Add(pairy);
      cutElements.Add(ptRange);
      cutElements.Add(rabs);
      cutElements.Add(matchlow);
      cutElements.Add(eta);
      cr->AddCutCombination(cutElements);

      // --- Disable some histo. Comment the one you want ---
      minvAnalysis->DisableHistograms("Eta");
      //minvAnalysis->DisableHistograms("Pt");
      //minvAnalysis->DisableHistograms("Y");

      // Adding the sub analysis
      task->AdoptSubAnalysis(minvAnalysis);
    }
  }


  // Addinf raw Tracklet analysis
  //task->AdoptSubAnalysis(nch1);
  task->AdoptSubAnalysis(nch);
  /*task->AdoptSubAnalysis(nch3);
    task->AdoptSubAnalysis(nch4);
    task->AdoptSubAnalysis(nch5);
    task->AdoptSubAnalysis(nch6);
    task->AdoptSubAnalysis(nch7);
    task->AdoptSubAnalysis(nch8);
    task->AdoptSubAnalysis(nch9);
    task->AdoptSubAnalysis(nch10);*/


  // ========= Configure binning =========

  // Class that hold the binnings.
  AliAnalysisMuMuBinning* binning = task->Binning();

  // --- Histogram binning ---

  // Binning versus "psi" == binning for track Histogram (see AliAnalysisMuMuMinv or AliAnalysisMuMuSingle to see the available binning )
  // binning->AddBin("psi","pt", 0.0, 1.0,"");
  // binning->AddBin("psi","pt", 1.0, 2.0,"");
  // binning->AddBin("psi","pt", 2.0, 3.0,"");
  // binning->AddBin("psi","pt", 3.0, 4.0,"");
  // binning->AddBin("psi","pt", 4.0, 5.0,"");
  // binning->AddBin("psi","pt", 5.0, 6.0,"");
  // binning->AddBin("psi","pt", 6.0, 7.0,"");
  // binning->AddBin("psi","pt", 7.0, 9.0,"");
  // binning->AddBin("psi","pt", 9.0,11.0,"");
  // binning->AddBin("psi","pt",11.0,15.0,"");
  // binning->AddBin("psi","pt",15.0,25.0,"");

  // Last arguments ("DH2" here) is just to add a term to Histo name to help differentiation in case of several binnings.
  // binning->AddBin("psi","ntrcorr",-1.5,-0.5,"D2H"); // 1 - 8
  // binning->AddBin("psi","ntrcorr",-0.5,0.5,"D2H"); // 1 - 8
  // binning->AddBin("psi","ntrcorr",0.5,8.5,"D2H"); // 1 - 8
  // binning->AddBin("psi","ntrcorr",8.5,13.5,"D2H"); // 9 - 14
  // binning->AddBin("psi","ntrcorr",13.5,16.5,"D2H"); // 20 - 24
  // binning->AddBin("psi","ntrcorr",16.5,20.5,"D2H"); // 25 - 33
  // binning->AddBin("psi","ntrcorr",20.5,24.5,"D2H"); // 34 - 41
  // binning->AddBin("psi","ntrcorr",24.5,28.5,"D2H"); // 34 - 41
  // binning->AddBin("psi","ntrcorr",28.5,32.5,"D2H"); // 42 - 50
  // binning->AddBin("psi","ntrcorr",32.5,38.5,"D2H"); // 51 - 59
  // binning->AddBin("psi","ntrcorr",38.5,44.5,"D2H"); // 51 - 59
  // binning->AddBin("psi","ntrcorr",44.5,54.5,"D2H"); // 51 - 59
  // binning->AddBin("psi","ntrcorr",54.5,74.5,"D2H"); // 100 - 199
  // binning->AddBin("psi","ntrcorr",74.5,140.5,"D2H"); // 100 - 199

  // binning->AddBin("psi","ntrcorr",140.5,200.0,"D2H"); // 100 - 199 // Only to check that FNorm matches with previous analysis

  // --- Centrality binning ---
  // You won't have to modify it i guess
  binning->AddBin("centrality","pp");


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
