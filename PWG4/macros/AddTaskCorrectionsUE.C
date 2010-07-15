void ConfigTaskCorrectionsUE(AliAnalysisTaskCorrectionsUE * uecorr );          // common config, extend with different cases
void SetTrackCuts(AliAnalysisTaskCorrectionsUE * uecorr, Char_t *ct);      //can be extended
void SetAnaTopology(AliAnalysisTaskCorrectionsUE * uecorr, Char_t *at);    //can be extended                  
void SetRegionType(AliAnalysisTaskCorrectionsUE * uecorr, Char_t *rt);     //can be extended        
void SetSorting(AliAnalysisTaskCorrectionsUE * uecorr, Char_t *sort);

AliAnalysisTaskCorrectionsUE *AddTaskCorrectionsUE(Char_t *jetBranch = "jets",Char_t *cuts = "ALICE", Char_t *anaTopology="LJ", Char_t *regionType="TRANSV", Char_t *sortBy="MSP", Bool_t proof=0)
{

   gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/lib/tgt_linux -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/JETAN  -I$ALICE_ROOT/JETAN/fastjet -I$ALICE_ROOT/PWG4 -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/CORRFW");
  
  // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCorrectionsUE", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCorrectionsUE", "This task requires an input event handler");
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
   
   AliAnalysisTaskCorrectionsUE* uecorr = new  AliAnalysisTaskCorrectionsUE(Form("UECorrections_%s",name.Data()));
   uecorr->SelectAODBranch(jb.Data());
   ConfigTaskCorrectionsUE(uecorr);
   SetTrackCuts(uecorr,cuts);
   SetAnaTopology(uecorr,anaTopology);
   SetRegionType(uecorr,regionType);
   SetSorting(uecorr,sortBy);

   if( jb.Contains("ICDF") ) { // Use internal "Charged Particle Jet CDF" instead of jets from AOD
      uecorr->SetUseChPartJet( kTRUE );
      uecorr->SetPtMinChPartJet(0.5);
   }

   //  ***** to be fixed *******
   /*
   if(jb.Length()>0){
     ueana->ReadDeltaAOD(kTRUE);
     ueana->SelectDeltaAODBranch(jb.Data());
   }
   */

   mgr->AddTask(uecorr);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_UE = 0;
   AliAnalysisDataContainer *coutput2_UE = 0;

   if (!proof){
   	if(jb.Length()==0)coutput1_UE = mgr->CreateContainer("histosCorrectionsUE", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_CorrectionsUE",AliAnalysisManager::GetCommonFileName()));
   	else coutput1_UE = mgr->CreateContainer(Form("histosCorrectionsUE_%s",name.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_CorrectionsUE_%s",AliAnalysisManager::GetCommonFileName(),name.Data()));
   
   	if(jb.Length()==0)coutput2_UE = mgr->CreateContainer("containersCorrectionsUE", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_CorrectionsUE",AliAnalysisManager::GetCommonFileName()));
   	else coutput2_UE = mgr->CreateContainer(Form("containersCorrectionsUE_%s",name.Data()), AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_CorrectionsUE_%s",AliAnalysisManager::GetCommonFileName(),name.Data()));
  } else {//simplify output structure in case of PROOF or it crashes while transferring files from master)
   	if(jb.Length()==0)coutput1_UE = mgr->CreateContainer("histosCorrectionsUE", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s",AliAnalysisManager::GetCommonFileName()));
   	else coutput1_UE = mgr->CreateContainer(Form("histosCorrectionsUE_%s",name.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s",AliAnalysisManager::GetCommonFileName()));
   
   	if(jb.Length()==0)coutput2_UE = mgr->CreateContainer("containersCorrectionsUE", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,Form("%s",AliAnalysisManager::GetCommonFileName()));
   	else coutput2_UE = mgr->CreateContainer(Form("containersCorrectionsUE_%s",name.Data()), AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,Form("%s",AliAnalysisManager::GetCommonFileName()));



  }

   mgr->ConnectInput  (uecorr, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (uecorr,     0, coutput1_UE );
   mgr->ConnectOutput (uecorr,     1, coutput2_UE );
   
   return uecorr;
}

//------------------------------------------------------------------------
void ConfigTaskCorrectionsUE(AliAnalysisTaskCorrectionsUE * uecorr){
  // common config,
  uecorr->SetDebugLevel(10); 
  uecorr->SetPtRangeInHist(100, 0., 100.);
}

//------------------------------------------------------------------------
SetTrackCuts(AliAnalysisTaskCorrectionsUE * uecorr, Char_t *ct){
  
  uecorr->SetFilterBit(16); // Standard TPC+ITS cuts 2009

  switch (ct) {
       case "ALICE":
       uecorr->SetTrackPtCut(0.1);
       uecorr->SetTrackEtaCut(0.9);
       break;
  
       case "CDF":
       uecorr->SetTrackPtCut(0.5);
       uecorr->SetTrackEtaCut(0.9); 
       break;

       default:
       Printf("\n >>>>>>> AddTaskCorrectionsUE Error Wrong Set of Track Cuts selected\n");
       break;
  }
}

//------------------------------------------------------------------------
SetAnaTopology(AliAnalysisTaskCorrectionsUE * uecorr, Char_t *at){

  switch (at) {
       case "LJ":  // leading jet
        uecorr->SetAnaTopology(1);
        uecorr->SetJet1EtaCut(0.5);
       break;

       case "LJCC":  // leading jet with charged cluster
        uecorr->SetAnaTopology(1);
        uecorr->SetJet1EtaCut(0.9);  // CDF allow jet radius extend ouside |eta|<1
       break;
  
       case "BB":  // back-to-back
        uecorr->SetAnaTopology(2);
        uecorr->SetJet1EtaCut(0.5);
        uecorr->SetJet2DeltaPhiCut(2.616);
       break;

       case "BE":  // back-to-back exclusive
        uecorr->SetAnaTopology(3);
        uecorr->SetJet1EtaCut(0.5);
        uecorr->SetJet2DeltaPhiCut(2.616);
        uecorr->SetJet3PtCut(10.);
       break;

       case "MP_eta05":
        uecorr->SetAnaTopology(4);
        uecorr->SetJet1EtaCut(0.9);
       break;

       case "MP_eta09":  // max Pt charged particle
         uecorr->SetAnaTopology(4);
         uecorr->SetJet1EtaCut(0.9);
       break;


       default:
       Printf("\n >>>>>>> AddTaskCorrectionsUE Error Wrong Analysis Topology selected\n");
       break;
  }

}

//------------------------------------------------------------------------
SetRegionType(AliAnalysisTaskCorrectionsUE * uecorr, Char_t *rt){

  switch (rt) {
       case "TRANSV":  
        uecorr->SetRegionType(1);        
       break;
  
       case "CONE":
        uecorr->SetRegionType(2);        
        uecorr->SetConeRadius(0.4);   
       break;

       default:
       Printf("\n >>>>>>> AddTaskCorrectionsUE Error Wrong Region Type selected\n");
       break;
  }
}  

//------------------------------------------------------------------------
SetSorting(AliAnalysisTaskCorrectionsUE * uecorr, Char_t *sort){
   
   switch (sort){
      // Minimum zone is selected per quantity type
      case "MQA":
       uecorr->SetPtSumOrdering(1);             
      break;
      
      // Minimum Zone is selected at the region having lowest pt summation of all tracks 
      case "MSP":
       uecorr->SetPtSumOrdering(2);
      break;
      
      // Minimum Zone is selected at the region having lowest average pt per track
      case "MAP":
       uecorr->SetPtSumOrdering(3);
      break;
      
      default:
       Printf("\n >>>>>>> AddTaskCorrectionsUE Error Wrong Sorting selected\n");     
       Printf("\n >>>>>>> AddTaskCorrectionsUE Sort by Set to MSP \n");
       uecorr->SetPtSumOrdering(2);
      break;
   }
}
