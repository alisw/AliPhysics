AliAnalysisTaskJetClusterKine *AddTaskJetClusterKine(char* bGen = "KINECHARGED",Char_t *jf = "ANTIKT", Float_t radius = 0.4, Int_t kWriteAOD = 1, char* deltaFile = "", Float_t ptTrackCut = 0.15, Float_t etaTrackWindow = 0.9, Float_t vertexWindow = 10.);

Float_t kPtTrackCutCl     = 0.15;
Float_t kTrackEtaWindowCl = 0.8;
Float_t kVertexWindowCl   = 10.0;


AliAnalysisTaskJetClusterKine *AddTaskJetClusterKine(char* bGen, Char_t *jf, Float_t radius, Int_t kWriteAOD, char *deltaFile, Float_t ptTrackCut, Float_t etaTrackWindow, Float_t vertexWindow){

 // Creates a jet fider task, configures it and adds it to the analysis manager.
   kPtTrackCutCl     = ptTrackCut;
   kTrackEtaWindowCl = etaTrackWindow;
   kVertexWindowCl   = vertexWindow;

   TString outputFile(deltaFile);
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr){
       ::Error("AddTaskJetClusterKine", "No analysis manager to connect to.");
       return NULL;
    }  

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if(!mgr->GetMCtruthEventHandler()){
      ::Error("AddTaskJetClusterKine", "This task requires an input MC event handler");
       return NULL;
    }

    TString typeGen(bGen);
    typeGen.ToUpper();

    // Create the task and configure it.
    //===========================================================================

    TString cAdd = "";
    cAdd += Form("%02d_",TMath::Nint(radius*10.));
    cAdd += Form("Cut%05d",TMath::Nint(1000.*kPtTrackCutCl));
    
    Printf("%s %s%s", typeGen.Data(), jf, cAdd.Data());

    AliAnalysisTaskJetClusterKine* clus = new  AliAnalysisTaskJetClusterKine(Form("JetCluster%s_%s%s",bGen,jf,cAdd.Data()));
      
    // or a config file
    clus->SetVtxCuts(kVertexWindowCl);

    if(typeGen.Contains("KINECHARGED")){
       clus->SetTrackTypeGen(AliAnalysisTaskJetClusterKine::kTrackKineCharged);
       clus->SetTrackPtCut(kPtTrackCutCl);
       clus->SetTrackEtaWindow(kTrackEtaWindowCl);

    }else if(typeGen.Contains("KINEFULL")){
       clus->SetTrackTypeGen(AliAnalysisTaskJetClusterKine::kTrackKineAll);
       clus->SetTrackPtCut(kPtTrackCutCl);
       clus->SetTrackEtaWindow(kTrackEtaWindowCl);
    }

   clus->SetRparam(radius);
   clus->SetGhostArea(0.005);
   clus->SetGhostEtamax(kTrackEtaWindowCl);

   clus->SetDebugLevel(0);

   switch (jf) {
   case "ANTIKT":
     clus->SetAlgorithm(2); // antikt from fastjet/JetDefinition.hh
     break;
   case "CA":
     clus->SetAlgorithm(1); // CA from fastjet/JetDefinition.hh
     break;
   case "KT":
     clus->SetAlgorithm(0); // kt from fastjet/JetDefinition.hh
     break;
   default:
     ::Error("AddTaskJetClusterKine", "Wrong jet finder selected\n");
     return 0;
   }

   TString nameOutArray =  Form("clusters%s_%s%s",bGen,jf,cAdd.Data()); //FF//
   if(kWriteAOD){
     if(outputFile.Length())clus->SetJetOutputFile(outputFile);
     Printf("Output branch: %s",nameOutArray.Data());//FF//  
     clus->SetJetOutputBranch(nameOutArray.Data());//FF//
     clus->SetJetOutputMinPt(0); // store only jets / clusters above a certain threshold
   }
   clus->SetJetOutputContainer(kWriteAOD); //0=no output 1=AOD 2=Exchange

   mgr->AddTask(clus);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_clus = mgr->CreateContainer(Form("cluster_%s_%s%s",bGen,jf,cAdd.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_cluster_%s_%s%s",AliAnalysisManager::GetCommonFileName(),bGen,jf,cAdd.Data()));

   mgr->ConnectInput  (clus, 0, mgr->GetCommonInputContainer());

   if(kWriteAOD==1){//FF//
      mgr->ConnectOutput (clus, 0, mgr->GetCommonOutputContainer());

   }

   mgr->ConnectOutput (clus, 1, coutput1_clus );



   if(kWriteAOD==2){//FF//
      AliAnalysisDataContainer *coutput2_clus = mgr->CreateContainer( nameOutArray.Data(), //??
                                           TClonesArray::Class(), 
                                           AliAnalysisManager::kExchangeContainer);
      mgr->ConnectOutput (clus, 2, coutput2_clus); //FF//
   }


   return clus;
}
