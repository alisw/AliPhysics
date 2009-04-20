
void ConfigTaskUE(AliAnalysisTaskUE * ueana); // common config, extend with different cases
                  
AliAnalysisTaskUE *AddTaskUE()
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
   
   AliAnalysisTaskUE* ueana = new  AliAnalysisTaskUE("Underlying Event");
   ConfigTaskUE(ueana);
   mgr->AddTask(ueana);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_UE = mgr->CreateContainer("histosUE", TList::Class(),AliAnalysisManager::kOutputContainer,"pwg4UE.root");
   
   mgr->ConnectInput  (ueana, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (ueana,     0, coutput1_UE );
   
   return ueana;
}

AliAnalysisTaskUE *AddTaskUE(AliAnalysisManager* mgr,AliAnalysisDataContainer *cinput)
{
  // This is only for running on PROOF with the old root version 5-22-00 
  // and the older version of the AF

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJets", "No analysis manager to connect to.");
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
 
   AliAnalysisTaskUE* ueana = new  AliAnalysisTaskUE("Underlying Event");
   ConfigTaskUE(ueana);
   
   AliAnalysisDataContainer *coutput1_UE = mgr->CreateContainer("histosUE", TList::Class(),AliAnalysisManager::kOutputContainer,"pwg4UE.root");
   

   mgr->ConnectInput  (ueana,  0, cinput);    
   mgr->ConnectOutput (ueana,     0, coutput1_UE );
   
   return ueana;

}  

void ConfigTaskUE(AliAnalysisTaskUE * ueana){
  // common config, extend with different cases
  Int_t anaType =1; 
  Int_t regType =1;
  Double_t jetEtaCut=0.2;
  Double_t trackPtCut=0.5; 
  Double_t trackEtaCut= 0.9; 
  Double_t rad=0.7; 
  Double_t deltaPhiCut = 2.616;
  
  ueana->SetDebugLevel(10); 
  ueana->SetPtRangeInHist(25, 0., 250.);
  ueana->SetAnaTopology(anaType);      
  ueana->SetRegionType(regType);        
  ueana->SetJet1EtaCut(jetEtaCut);     
  ueana->SetTrackPtCut(trackPtCut); 
  ueana->SetPtSumOrdering(2);
  ueana->SetConeRadius(rad);   
  ueana->SetTrackEtaCut(trackEtaCut);
  ueana->SetJet2DeltaPhiCut(deltaPhiCut);
}
