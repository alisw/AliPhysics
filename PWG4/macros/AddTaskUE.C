

void ConfigTaskUE(AliAnalysisTaskUE * ueana );          // common config, extend with different cases
void SetTrackCuts(AliAnalysisTaskUE * ueana, Char_t *ct);              //can be extended                 
void SetAnaTopology(AliAnalysisTaskUE * ueana, Char_t *at);    //can be extended                  
void SetRegionType(AliAnalysisTaskUE * ueana, Char_t *rt);     //can be extended                  

AliAnalysisTaskUE *AddTaskUE(Char_t *jetBranch = "jets",Char_t *cuts = "ALICE", Char_t *anaTopology="LJ", Char_t *regionType="TRANSV")
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskUE", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskUE", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   
   TString jb(jetBranch);
   TString ct(cuts);
   TString at(anaTopology);
   TString rt(regionType);
   TString name(Form("%s_%s_%s_%s",jb.Data(),ct.Data(),at.Data(),rt.Data()));
   
   AliAnalysisTaskUE* ueana = new  AliAnalysisTaskUE(Form("UEAnalysis_%s",name.Data()));
   ueana->SelectAODBranch(jb.Data());
   ConfigTaskUE(ueana);
   SetTrackCuts(ueana,cuts);
   SetAnaTopology(ueana,anaTopology);
   SetRegionType(ueana,regionType);

   //  ***** to be fixed *******
   /*
   if(jb.Length()>0){
     ueana->ReadDeltaAOD(kTRUE);
     ueana->SelectDeltaAODBranch(jb.Data());
   }
   */

   mgr->AddTask(ueana);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_UE = 0;
   if(jb.Length()==0)coutput1_UE = mgr->CreateContainer("histosUE", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_UE",AliAnalysisManager::GetCommonFileName()));
   else coutput1_UE = mgr->CreateContainer(Form("histosUE_%s",name.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_UE_%s",AliAnalysisManager::GetCommonFileName(),name.Data()));
   
   mgr->ConnectInput  (ueana, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (ueana,     0, coutput1_UE );
   
   return ueana;
}

void ConfigTaskUE(AliAnalysisTaskUE * ueana){
  // common config,
  ueana->SelectTrigger(AliAnalysisHelperJetTasks::kMB1);
  ueana->SetDebugLevel(0); 
  ueana->SetPtRangeInHist(25, 0., 50.);
  ueana->SetPtSumOrdering(2);
}

//------------------------------------------------------------------------
SetTrackCuts(AliAnalysisTaskUE * ueana, Char_t *ct){
  
  ueana->SetFilterBit(16); // ITS refit

  switch (ct) {
       case "ALICE":
       ueana->SetTrackPtCut(0.1);
       ueana->SetTrackEtaCut(0.9);
       break;
  
       case "CDF":
       ueana->SetTrackPtCut(0.4);
       ueana->SetTrackEtaCut(1.2); // meaningful only for pure MC
        break;

       default:
       Printf("\n >>>>>>> AddTaskUE Error Wrong Set of Track Cuts selected\n");
       break;
  }
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
SetAnaTopology(AliAnalysisTaskUE * ueana, Char_t *at){

  switch (at) {
       case "LJ":  // leading jet
        ueana->SetAnaTopology(1);      
        ueana->SetJet1EtaCut(0.5);     
       break;
  
       case "BB":  // back-to-back
        ueana->SetAnaTopology(2);      
        ueana->SetJet1EtaCut(0.5);     
        ueana->SetJet2DeltaPhiCut(2.616);
        break;

       default:
       Printf("\n >>>>>>> AddTaskUE Error Wrong Analysis Topology selected\n");
       break;
  }

}

//------------------------------------------------------------------------
SetRegionType(AliAnalysisTaskUE * ueana, Char_t *rt){

  switch (rt) {
       case "TRANSV":  
        ueana->SetRegionType(1);        
       break;
  
       case "CONE":
        ueana->SetRegionType(2);        
        ueana->SetConeRadius(0.4);   
        break;

       default:
       Printf("\n >>>>>>> AddTaskUE Error Wrong Region Type selected\n");
       break;
  }
}  

