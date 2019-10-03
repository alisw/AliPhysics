AliAnalysisTaskPrepareInputForEmbedding *AddTaskPrepareInputForEmbedding(
   TString njet = "Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme",
   /*TString njetClosest = "Jet_AKTChargedR040_MCParticlesSelected_pT0150_E_scheme",*/
   Double_t R = 0.4,
   Double_t jetptcut = 1,
   Double_t jetareacut = 0.6,
   TString ntrack = "PicoTracks", 
   TString accType = "TPC", 
   Bool_t  useLeadingJet = kFALSE,
   Bool_t  usedHCtagging = kFALSE){
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr){
      Error("AddTaskPrepareInputForEmbedding","No analysis manager found.");
      return 0;
   }
   if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskPrepareInputForEmbedding", "This task requires an input event handler");
      return NULL;
    }
    TString wagonName = Form("PrepareInputForEmbedding%s%s", useLeadingJet ? "L" : "", usedHCtagging ? "HC" : "");
    AliAnalysisTaskPrepareInputForEmbedding *task = new AliAnalysisTaskPrepareInputForEmbedding(wagonName);
    
    task->SetLeadingJetOnly(useLeadingJet);
    task->SetDetHardCoreTagging(usedHCtagging);
    
    AliParticleContainer *partCont  = task->AddParticleContainer(ntrack);
       
    AliJetContainer *jetCont  = task->AddJetContainer(njet, accType, R);
    if(jetCont){
       //jetCont->SetRhoName(rhoname);
       //jetCont->SetRhoMassName(rhoMname);
       jetCont->ConnectParticleContainer(partCont);
       jetCont->SetJetPtCut(jetptcut);
       jetCont->SetPercAreaCut(jetareacut);
       
    }
    //AliJetContainer *jetContClosest  = task->AddJetContainer(njetClosest, accType, R);
    
    mgr->AddTask(task);
    
    //Connnect input
    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
    //Connect output
    TString contName(wagonName);
    TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    contName.Append("Tree");
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    
    mgr->ConnectOutput(task,1,coutput1);
    mgr->ConnectOutput(task,2,coutput2);

    return task;
}
