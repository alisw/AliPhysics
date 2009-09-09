AliAnalysisTaskESDMCLabelAddition *AddTaskESDMCLabelAddition(Bool_t useKineFilter=kTRUE)
{
// Creates a filter task and adds it to the analysis manager.
// This file allows the creation of MC labels (based on the code of Philippe P.)

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskESDMCLabelAddition", "No analysis manager to connect to.");
      return NULL;
   }   
   
   TString type = mgr->GetInputEventHandler()->GetDataType();
   // Check if MC handler is connected in case kine filter requested
   AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
   if (!mcH && useKineFilter) {
      ::Error("AddTaskESDFilter", "No MC handler connected while kine filtering requested");
      return NULL;
   }   
   
   if (useKineFilter) {
      AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Kine Filter");
      mgr->AddTask(kinefilter);
   }   


   // Create the task, add it to the manager and configure it.
   //===========================================================================   
   // Barrel tracks filter
   AliAnalysisTaskESDMCLabelAddition *ESDMCLabeltask = new AliAnalysisTaskESDMCLabelAddition("ESD MC Labels addition");
   mgr->AddTask(ESDMCLabeltask);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (ESDMCLabeltask,  0, mgr->GetCommonInputContainer());

   return ESDMCLabeltask;
}   
