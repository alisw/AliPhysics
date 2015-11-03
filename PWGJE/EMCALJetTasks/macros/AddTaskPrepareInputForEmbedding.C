AliAnalysisTaskPrepareInputForEmbedding *AddTaskPrepareInputForEmbedding(
   TString njetAreaSub = "Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme", 
   TString njetConstSub = "Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeConstSub", 
   Double_t R = 0.4, 
   TString ntrack = "PicoTracks", 
   TString accType = "TPC", 
   TString rhoname = "", 
   TString rhoMname = "",
   Bool_t  useLeadingJet = kFALSE){
   
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
    TString wagonName = Form("PrepareInputForEmbedding%s", useLeadingJet ? "L" : "");
    AliAnalysisTaskPrepareInputForEmbedding *task = new AliAnalysisTaskPrepareInputForEmbedding(wagonName);
    
    task->SetLeadingJetOnly(useLeadingJet);
    
    AliParticleContainer *partCont  = task->AddParticleContainer(ntrack);
    AliParticleContainer *partContConst  = task->AddParticleContainer(Form("%s_%s",ntrack.Data(), njetConstSub.Data()));
       
    AliJetContainer *jetContA  = task->AddJetContainer(njetAreaSub, accType, R);
    jetContA->SetRhoName(rhoname);
    jetContA->SetRhoMassName(rhoMname);
    jetContA->ConnectParticleContainer(partCont);
    
    AliJetContainer *jetContC  = task->AddJetContainer(njetConstSub, accType, R);
    jetContC->SetRhoName(rhoname);
    jetContC->SetRhoMassName(rhoMname);
    jetContC->ConnectParticleContainer(partContConst);
    
    //Connnect input
    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
    //Connect output
    TString contName(wagonName);
    TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    mgr->ConnectOutput(task,1,coutput1);

    return task;
}
