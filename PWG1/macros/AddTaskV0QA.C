///////////////////////////////////////////////////////////////////////////////
// Macro to setup AliAnalysisTaskV0QA for
// V0 performance QA to run on PWG1 QA train.
// ESD and MC input handlers must be attached to AliAnalysisManager
// 
// Output:
// V0QA.root file with V0QA for gammas,K0s,Lambdas, antilambdas is created.
//
// The file contains 4 THnSparse generic histograms which
// have to be analysed using the macros 
// PWG1/macros/plotSparseK0.C
// PWG1/macros/plotSparseGamma.C
// PWG1/macros/plotSparseL.C
// PWG1/macros/plotSparseAL.C
// 
// 14.10.09 A. Marin  a.marin@gsi.de
///////////////////////////////////////
 
AliAnalysisTaskV0QA *AddTaskV0QA(Bool_t bUseMCInfo=kTRUE)
{
// Creates a filter task and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskV0QA", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // This task requires an ESD input handler and an AOD output handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD")) {
      ::Error("AddTaskV0QA", "V0QA task needs the manager to have an ESD input handler.");
      return NULL;
   }   

   // Check if MC handler is connected in case kine filter requested
   AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
   if (!mcH && bUseMCInfo)) {
      ::Error("AddTaskV0QA", "No MC handler connected");
      return NULL;
   }   
   
   // Create the task, add it to the manager and configure it.
   //===========================================================================   
   // Barrel tracks filter
   AliAnalysisTaskV0QA * v0QA = new AliAnalysisTaskV0QA("V0QA");
   mgr->AddTask(v0QA);


   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (v0QA,  0, mgr->GetCommonInputContainer());

   // Create containers for output

   AliAnalysisDataContainer *coutput_v0QA =  
   mgr->CreateContainer("V0QA", TList::Class(), 
			AliAnalysisManager::kOutputContainer, Form("%s:%s", mgr->GetCommonFileName(),v0QA->GetName()));
   mgr->ConnectOutput(v0QA, 1, coutput_v0QA);
 	

   return v0QA;
}   
