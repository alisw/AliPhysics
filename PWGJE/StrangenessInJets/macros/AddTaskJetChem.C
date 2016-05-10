AliAnalysisTaskJetChem *AddTaskJetChem(const char* recJetsBranch = "clustersAOD_ANTIKT02_B2_Filter00768_Cut00150_Skip00", TString cutVar = "noCutVar", Int_t eventClass = 1, Int_t K0type = AliAnalysisTaskJetChem::kOffl, Int_t Latype = AliAnalysisTaskJetChem::kOffl, Int_t ALatype = AliAnalysisTaskJetChem::kOffl, Bool_t IsArmenterosSelected = kTRUE, Bool_t IsJetPtBiasSelected = kTRUE, Double_t jetradius = 0.2, Double_t V0EtaCut = 0.7, Double_t jetEtaCut = 0.5, Bool_t IsMC = kFALSE, Double_t DeltaVtxZCut = 0.1, Int_t filtermask = 768, Int_t fdebug = -1, Bool_t useExtraOnlyTracks = 0, Bool_t useExtraTracks = 0, Bool_t useExtraJetPt = kFALSE)
{
  // Creates a JetChem task,
  // configures it and adds it to the analysis manager.
  //for cut variations, set TString cutVar either to "noCutVar", "CPAMin", "CPAMax", "DCADMin", "DCADMax", "LifetMin" or "LifetMax"
   
  
  // ** Parameters **
  // (char) recJetsBranch: branch in AOD for (reconstructed) jets
  // (char) genJetsBranch: branch in AOD for (generated) jets
  // (char) jetType: "AOD"   jets from recJetsBranch
  //                 "AODMC" jets from genJetsBranch
  //                 "KINE"  jets from PYCELL
  //                  +"b" (e.g. "AODb") jets with acceptance cuts
  // (char) trackType: "AOD"     reconstructed tracks from AOD filled by ESD filter (choose filter mask!)
  //                   "AODMC"   MC tracks from AOD filled by kine filter
  //                   "KINE"    kine particles from MC event 
  //                   +"2" (e.g. "AOD2")  charged tracks only
  //                   +"b" (e.g. "AOD2b") with acceptance cuts
  // (UInt_t) filterMask: select filter bit of ESD filter task
  
  


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskJetChem", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetChem", "This task requires an input event handler");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD", or "MC"?
  Printf("Data Type: %s", type.Data());
  
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskJetChem *task = new AliAnalysisTaskJetChem("TaskJetChem");

  Int_t debug = fdebug; // debug level
  if(debug>=0) task->SetDebugLevel(debug);
  
  char* genJetsBranch = "";
  //recJetsBranch = "";

  TString branchRecJets(recJetsBranch);
  TString branchGenJets(genJetsBranch);

  task->SetBranchEmbeddedJets("");//insert string for embedding in wagon configuration, method declared and defined in FragmentationFunction task
  task->SetBranchGenJets("");//insert string for embedding in wagon configuration, method declared and defined in FragmentationFunction task
  
  if(useExtraTracks)task->SetMatchMode(1);//default = 1: is det.level rec. - det.level PYTHIA matching, '2' is matching from 'det.level rec. extra' jets to 'particle level PYTHIA' jets
  task->SetUseStandardV0s(kTRUE);//fill extra branch also with standard v0s for UE V0 subtraction

  fJetAreaMin = 0.6*TMath::Pi()*jetradius*jetradius;//calculate jetareamin cut value for FF task
  task->SetJetMinArea(fJetAreaMin);//cut on jet area, applied together with all other jet cuts in jet finding by AliAnalysisTaskFragmentationFunction.cxx
  task->SetCutJetEta(jetEtaCut);
  task->SetDeltaZVertexCut(DeltaVtxZCut);
  task->SetBranchRecBackClusters("clustersAOD_KT04_B0_Filter00768_Cut00150_Skip00"); 

  //task->SetEventSelectionMask(AliVEvent::kMB); //for 2010 Pb-Pb data !!
  task->SetEventSelectionMask(AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral); //event selection for 2011 Pb-Pb data
  task->SetEventClass(eventClass);
  task->SetK0Type(K0type);
  task->SetLaType(Latype); 
  task->SetALaType(ALatype); 
  task->SetSelectArmenteros(IsArmenterosSelected);
  task->SetAnalysisMC(IsMC); // 0: real data, 1: MC data
  task->SetCutV0Rap(0.);//not applied as default
  task->SetCutDeltaREmbedded(0.); //for standard tracks
  task->SetCutFractionPtEmbedded(1.); //for standard tracks
  task->SetUseNJEvents(kFALSE);//Embedding into all events (kFALSE) or only in no-jet events (kTRUE)

  if(K0type == AliAnalysisTaskJetChem::kOnFlyPrim || AliAnalysisTaskJetChem::kOfflPrim) task->SetFilterMaskK0(768);
  if(Latype == AliAnalysisTaskJetChem::kOnFlyPrim || AliAnalysisTaskJetChem::kOfflPrim) task->SetFilterMaskLa(768);
  if(ALatype == AliAnalysisTaskJetChem::kOnFlyPrim || AliAnalysisTaskJetChem::kOfflPrim) task->SetFilterMaskALa(768);

  task->SetFFRadius(jetradius); //jet cone size
  task->SetFilterMask(filtermask);//2011 Track FilterMask

  //Cuts---------------------------------

  task->SetTrackCuts(0.15, -0.9, 0.9, 0., 2*TMath::Pi());// (pt Cut, daughtertrack rap's, phi min max cuts)
  task->SetJetCuts(5., (-1)*jetEtaCut, jetEtaCut, 0., 2*TMath::Pi());//(jet pt Cut, jet acceptance, phi min max cuts)
  task->SetCuttrackPosEta(0.8);
  task->SetCuttrackNegEta(0.8);
  task->SetCutV0Eta(V0EtaCut); //pseudorapidity cut, don't use 0.8, because too many tracks would fall out of the acceptance; recommended cut for jet analysis of strange particles: 0.75
  task->SetCosOfPointingAngleK0(0.998);//as default the same value for all V0 types
  task->SetCosOfPointingAngleLa(0.998);
  task->SetAcceptKinkDaughters(kFALSE);//accept kink daughters -> dont use this cut anymore
  task->SetRequireTPCRefit(kTRUE);
  task->SetCutV0DecayMin(0.);//multiples of ctau, cut on 2D decay distance over transverse mom. (for K0s, Lambda, Antilambda)
  task->SetCutV0DecayMax(5.);//multiples of ctau (for K0s, Lambda, Antilambda) Lee Barnby uses 3.0, use 5.0!!!!!
  task->SetCutDcaV0Daughters(1.);//cut value in multiples of sigma default: 1.
  task->SetCutDcaPosToPrimVertex(0.1); //cut value in cm 
  task->SetCutDcaNegToPrimVertex(0.1); //cut value in cm
  task->SetCutV0RadiusMin(5.);//in cm previous value was 0.9 cm
  task->SetCutV0RadiusMax(100.);//in cm
  task->SetCutBetheBloch(3.);//in units of sigma


  if(useExtraTracks)task->UseExtraTracks();
  if(useExtraOnlyTracks)task->UseExtraonlyTracks();
  if(useExtraJetPt)task->SetUseExtraJetPt(kTRUE);//Use smeared jet pt for MC truth reference
  task->SetUseEmbeddedJetPt(kFALSE);

  //task->SetCutRatioTPC(0.8);//Cut on Ratio of crossed Rows over findable clusters in TPC -> not used anymore by Strangeness PAG group
   //task->SetCuttrackPosNcls(70);
  //task->SetCuttrackNegNcls(70);
  //task->SetCuttrackPosRap(100000.0);
  //task->SetCuttrackNegRap(100000.0);
  //task->SetCutV0Rap(0.5);
  //task->SetCutV0totMom(10000.);//tot Mom of V0s 

  //Armenteros Cut:
 
  if(IsArmenterosSelected == 1){
    task->SetCutArmenteros(0.2);
  } else {
    task->SetCutArmenteros(0.);
  }
  
  if(IsJetPtBiasSelected == 1){
    task->SetFFMinLTrackPt(5.);//at least one track must have a track-pt higher or equal than this JetMinPt value
  } else {
    task->SetFFMinLTrackPt(-1.);
  }
  
  //------------------------------------
  // Define histo bins
  task->SetFFHistoBins();
  task->SetQAJetHistoBins();
  task->SetQATrackHistoBins();
  task->SetFFInvMassHistoBins();
  task->SetFFInvMassLaHistoBins();
  
  mgr->AddTask(task);
  
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   TString strK0type;
   if(K0type ==  AliAnalysisTaskJetChem::kOnFly)     strK0type = "OnFly";
   if(K0type ==  AliAnalysisTaskJetChem::kOnFlyPID)  strK0type = "OnFlyPID";
   if(K0type ==  AliAnalysisTaskJetChem::kOnFlydEdx) strK0type = "OnFlydEdx";
   if(K0type ==  AliAnalysisTaskJetChem::kOnFlyPrim) strK0type = "OnFlyPrim";
   if(K0type ==  AliAnalysisTaskJetChem::kOffl)      strK0type = "Offl";
   if(K0type ==  AliAnalysisTaskJetChem::kOfflPID)   strK0type = "OfflPID";
   if(K0type ==  AliAnalysisTaskJetChem::kOffldEdx)  strK0type = "OffldEdx";
   if(K0type ==  AliAnalysisTaskJetChem::kOfflPrim)  strK0type = "OfflPrim";

   TString strLatype;
   if(Latype ==  AliAnalysisTaskJetChem::kOnFly)     strLatype = "OnFly";
   if(Latype ==  AliAnalysisTaskJetChem::kOnFlyPID)  strLatype = "OnFlyPID";
   if(Latype ==  AliAnalysisTaskJetChem::kOnFlydEdx) strLatype = "OnFlydEdx";
   if(Latype ==  AliAnalysisTaskJetChem::kOnFlyPrim) strLatype = "OnFlyPrim";
   if(Latype ==  AliAnalysisTaskJetChem::kOffl)      strLatype = "Offl";
   if(Latype ==  AliAnalysisTaskJetChem::kOfflPID)   strLatype = "OfflPID";
   if(Latype ==  AliAnalysisTaskJetChem::kOffldEdx)  strLatype = "OffldEdx";
   if(Latype ==  AliAnalysisTaskJetChem::kOfflPrim)  strLatype = "OfflPrim";

   TString strALatype;
   if(ALatype ==  AliAnalysisTaskJetChem::kOnFly)     strALatype = "OnFly";
   if(ALatype ==  AliAnalysisTaskJetChem::kOnFlyPID)  strALatype = "OnFlyPID";
   if(ALatype ==  AliAnalysisTaskJetChem::kOnFlydEdx) strALatype = "OnFlydEdx";
   if(ALatype ==  AliAnalysisTaskJetChem::kOnFlyPrim) strALatype = "OnFlyPrim";
   if(ALatype ==  AliAnalysisTaskJetChem::kOffl)      strALatype = "Offl";
   if(ALatype ==  AliAnalysisTaskJetChem::kOfflPID)   strALatype = "OfflPID";
   if(ALatype ==  AliAnalysisTaskJetChem::kOffldEdx)  strALatype = "OffldEdx";
   if(ALatype ==  AliAnalysisTaskJetChem::kOfflPrim)  strALatype = "OfflPrim";
   
   if((IsArmenterosSelected == 0) && (IsJetPtBiasSelected == 0)){
       TString listName1(Form("PWG4_JetChem_%s_%s_cl%d_%s",branchRecJets.Data(),strK0type.Data(),eventClass,cutVar.Data()));
       AliAnalysisDataContainer *coutput_JetChem = mgr->CreateContainer(listName1, 
									TList::Class(),
									AliAnalysisManager::kOutputContainer,
									Form("%s:PWG4_zimmerma_JetChem",AliAnalysisManager::GetCommonFileName()));
       if(useExtraOnlyTracks){listName1 += "_exonly";} 
       if(useExtraTracks){listName1 += "_extra";}
       //if(useExtraTracks){listName1 += "_extra";}
       if(useExtraJetPt){listName1 += "_extraJetPt";}
   }
   
   if((IsArmenterosSelected == 1) && (IsJetPtBiasSelected == 0)){
     TString listName2(Form("PWG4_JetChem_%s_%s_cl%d_Armenteros_%s",branchRecJets.Data(),strK0type.Data(),eventClass,cutVar.Data()));
       AliAnalysisDataContainer *coutput_JetChem = mgr->CreateContainer(listName2, 
									TList::Class(),
									AliAnalysisManager::kOutputContainer,
									Form("%s:PWG4_zimmerma_JetChem_Armenteros",AliAnalysisManager::GetCommonFileName()));
       if(useExtraOnlyTracks){listName2 += "_exonly";} 
       if(useExtraTracks){listName2 += "_extra";}   
       if(useExtraJetPt){listName2 += "_extraJetPt";}
   }
      
   
   if((IsArmenterosSelected == 0) && (IsJetPtBiasSelected == 1)) {
     TString listName3(Form("PWG4_JetChem_%s_%s_cl%d_JetPtbias_%s",branchRecJets.Data(),strK0type.Data(),eventClass,cutVar.Data()));
     AliAnalysisDataContainer *coutput_JetChem = mgr->CreateContainer(listName3, 
								      TList::Class(),
								      AliAnalysisManager::kOutputContainer,
								      Form("%s:PWG4_zimmerma_JetChem_JetPtBias",AliAnalysisManager::GetCommonFileName()));
     if(useExtraOnlyTracks){listName3 += "_exonly";} 
     if(useExtraTracks){listName3 += "_extra";}
     if(useExtraJetPt){listName3 += "_extraJetPt";}
   }
   
   if((IsArmenterosSelected == 1) && (IsJetPtBiasSelected == 1)) {
     TString listName4(Form("PWG4_JetChem_%s_%s_cl%d_Armenteros_JetPtBias_%s",branchRecJets.Data(),strK0type.Data(),eventClass,cutVar.Data()));
     AliAnalysisDataContainer *coutput_JetChem = mgr->CreateContainer(listName4, 
								      TList::Class(),
								      AliAnalysisManager::kOutputContainer,
								      Form("%s:PWG4_zimmerma_JetChem_Armenteros_JetPtBias",AliAnalysisManager::GetCommonFileName()));
     if(useExtraOnlyTracks){listName4 += "_exonly";} 
     if(useExtraTracks){listName4 += "_extra";}
     if(useExtraJetPt){listName4 += "_extraJetPt";}
   } 
   


   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, coutput_JetChem);
   
   return task;
}
