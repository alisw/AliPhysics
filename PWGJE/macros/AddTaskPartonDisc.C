AliAnalysisTaskPartonDisc* AddTaskPartonDisc(const char* bRec = "jetsAOD_UA104_B0_Filter00128_Cut01000",const char* bRec2 = "jetsAOD_SISCONE04_B0_Filter00128_Cut00150",const char* bGen = "jetsAODMC2_UA104_B0_Filter00000_Cut01000",UInt_t filterMask = 128, Int_t iPhysicsSelectionFlag = AliVEvent::kMB, Int_t option=1, Int_t ntx=90, Double_t jetrad = 0.4, Double_t trackpTcut = 1.0, Double_t incrad = 0.0, Double_t minpTUM = 0.150, Double_t maxpTUM = 0.900, Bool_t skipsingletr = kFALSE, Bool_t notextendExcl = kFALSE,  Double_t sqrts = 7000., Double_t minpTMc = 0.150, Bool_t enaEtaRest = kFALSE)
{
  ////options////
  // 1: Real pp data
  // 2: MC
  // 3: PbPb

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddTaskPartonDisc", "No analysis manager to connect to.");
      return NULL;
    }  
  
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddTaskPartonDisc", "This task requires an input event handler");
      return NULL;
    }
  
  AliAnalysisTaskPartonDisc* taskPD = new  AliAnalysisTaskPartonDisc(Form("PartonDisc%s-%s",bRec,bGen));
  if(option==1)
    {
      // Default settings in constructor
    }
  if(option==2)
    {
      taskPD->SetAODwithMC(kTRUE);
      taskPD->SetPhojetMC(kTRUE); 
      taskPD->SetMCBranch(bGen);
      taskPD->SetFlavorRadius(0.3); //default is 0.3 
      taskPD->SetSqrtS(sqrts);
      taskPD->SetMinPtMC(minpTMc); // pT min for mc particles (V0 mips)
    }
  if(option==3)
    taskPD->SetHIEvent(kTRUE);

  if(iPhysicsSelectionFlag)
    taskPD->SelectCollisionCandidates(iPhysicsSelectionFlag);
  taskPD->SetAODMCInput(kTRUE); // true when running over AODs
  taskPD->SetRecBranch(bRec);
  taskPD->SetSecondRecBranch(bRec2);
  taskPD->SetXNtX(ntx);  
  taskPD->SetJetRadius(jetrad); 
  taskPD->SetFilterBitTracks(filterMask);
  taskPD->SetMinPtTrackCut(trackpTcut);
  taskPD->SetMinPtUE(minpTUM); // for underlying event multiplicity counting
  taskPD->SetMaxPtUE(maxpTUM); // for underlying event multiplicity counting
  taskPD->SetIncreaseOfExclusionR(incrad); // increase to R+incrad
  taskPD->ForceNotUseTrackRefs(kFALSE); // kTRUE para ignorar las tracks refs(diferente bit del usado en jet finding) 
  taskPD->NotExtendDiJetExclusion(notextendExcl); // kFALSE->extend radius in dijet area, kTRUE->don't extend the radius in dijet area 
  taskPD->ForceSkipSingleTrackJets(skipsingletr); // kTRUE to force to skip single track jets
  task->SetEnableJetEtaRestriction(enaEtaRest); // If increase of exclusion radius =!0 -> kTRUE, if not kFALSE
  mgr->AddTask(taskPD);
    
  AliAnalysisDataContainer *coutput_PartDisc = mgr->CreateContainer(Form("taskPD_%s_%s",bRec,bGen),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWJE_taskPD_%s_%s",AliAnalysisManager::GetCommonFileName(),bRec,bGen));

  mgr->ConnectInput  (taskPD, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (taskPD, 0, mgr->GetCommonOutputContainer()); // comment to run local
  mgr->ConnectOutput (taskPD, 1, coutput_PartDisc );
  
  return taskPD;
}
