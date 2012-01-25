

void ConfigTaskUE(AliAnalysisTaskUE * ueana );          // common config, extend with different cases
void SetTrackCuts(AliAnalysisTaskUE * ueana, Char_t *ct);      //can be extended
void SetAnaTopology(AliAnalysisTaskUE * ueana, Char_t *at);    //can be extended                  
void SetRegionType(AliAnalysisTaskUE * ueana, Char_t *rt);     //can be extended        
void SetSorting(AliAnalysisTaskUE * ueana, Char_t *sort);

AliAnalysisTaskUE *AddTaskUE(Char_t *jetBranch = "jets",Char_t *cuts = "ALICE", Char_t *anaTopology="LJ", Char_t *regionType="TRANSV", Char_t *sortBy="MSP")
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
   TString sort(sortBy);
   TString name(Form("%s_%s_%s_%s_%s",jb.Data(),ct.Data(),at.Data(),rt.Data(),sort.Data()));
   
   AliAnalysisTaskUE* ueana = new  AliAnalysisTaskUE(Form("UEAnalysis_%s",name.Data()));
   ueana->SelectAODBranch(jb.Data());
   ConfigTaskUE(ueana);
   SetTrackCuts(ueana,cuts);
   SetAnaTopology(ueana,anaTopology);
   SetRegionType(ueana,regionType);
   SetSorting(ueana,sortBy);

   if( jb.Contains("ICDF") ) { // Use internal "Charged Particle Jet CDF" instead of jets from AOD
      ueana->SetUseChPartJet( kTRUE );
      ueana->SetPtMinChPartJet(0.5);
   }

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
  ueana->SetDebugLevel(0); 
  ueana->SetPtRangeInHist(100, 0., 100.);
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
       ueana->SetTrackPtCut(0.5);
       ueana->SetTrackEtaCut(0.9); // meaningful only for pure MC
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

       case "LJCC":  // leading jet with charged cluster
        ueana->SetAnaTopology(1);
        ueana->SetJet1EtaCut(0.9);  // CDF allow jet radius extend ouside |eta|<1
       break;
  
       case "BB":  // back-to-back
        ueana->SetAnaTopology(2);
        ueana->SetJet1EtaCut(0.5);
        ueana->SetJet2DeltaPhiCut(2.616);
       break;

       case "BE":  // back-to-back exclusive
        ueana->SetAnaTopology(3);
        ueana->SetJet1EtaCut(0.5);
        ueana->SetJet2DeltaPhiCut(2.616);
        ueana->SetJet3PtCut(10.);
       break;

       case "MP_eta05":
        ueana->SetAnaTopology(4);
        ueana->SetJet1EtaCut(0.9);
       break;

       case "MP_eta09":  // max Pt charged particle
         ueana->SetAnaTopology(4);
         ueana->SetJet1EtaCut(0.9);
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

//------------------------------------------------------------------------
SetSorting(AliAnalysisTaskUE * ueana, Char_t *sort){
   
   switch (sort){
      // Minimum zone is selected per quantity type
      case "MQA":
       ueana->SetPtSumOrdering(1);             
      break;
      
      // Minimum Zone is selected at the region having lowest pt summation of all tracks 
      case "MSP":
       ueana->SetPtSumOrdering(2);
      break;
      
      // Minimum Zone is selected at the region having lowest average pt per track
      case "MAP":
       ueana->SetPtSumOrdering(3);
      break;
      
      default:
       Printf("\n >>>>>>> AddTaskUE Error Wrong Sorting selected\n");     
       Printf("\n >>>>>>> AddTaskUE Sort by Set to MSP \n");
       ueana->SetPtSumOrdering(2);
      break;
   }
}
