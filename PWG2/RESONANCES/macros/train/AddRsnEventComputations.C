//
// Configuration script for monitor task with 2010 runs
//
// It contains a class definition where the cuts for each object
// are defined separately, functions are initialized and so on.
// This is used in the main function (named after the file name),
// which is called by the 'AddTask' function.
//

#if !defined(__CINT__) || defined(__MAKECINT__)

   #include "TString.h"
   #include "AliAnalysisManager.h"
   #include "AliRsnValue.h"
   #include "AliRsnFunction.h"
   #include "AliRsnCutValue.h"
   #include "AliRsnCutPrimaryVertex.h"
   #include "AliRsnAnalysisTask.h"

#endif

Bool_t AddRsnEventComputations(Bool_t isMC, const char *options = "", const char *taskName = "RSNtask")
{
   // ==================================================================================================================
   // == PRELIMINARY OPERATIONS ========================================================================================
   // ==================================================================================================================
   
   // retrieve task from manager, using its name
   AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
   AliRsnAnalysisTask *task = (AliRsnAnalysisTask*)mgr->GetTask(taskName);
   if (!task) {
      Error("RsnConfigMonitor", "Task not found");
      return kFALSE;
   }
   
   TString opt(options);
   opt.ToUpper();
   opt.ReplaceAll(" ", "");
   
   Bool_t central    = opt.Contains("CENT");
   Bool_t peripheral = opt.Contains("PERI");
   
   // ==================================================================================================================
   // == EVENT CUTS ====================================================================================================
   // ==================================================================================================================

   // event cuts are added directly to a cutSet in the task
   // we create all and then add thos which are needed
   
   // primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   // we switch on the check for pileup
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
   cutVertex->SetCheckPileUp(kTRUE);
      
   // centrality:
   // - 2nd      argument --> one of the centrality evaluation methods
   // - 3rd, 4th argument --> centrality ranges in percentile (0-10 for central, 60-70 for peripheral)
   AliRsnCutValue *cutCentrality = 0x0;
   if (central) 
      cutCentrality = new AliRsnCutValue("cutCentral", AliRsnValue::kEventCentralityV0,  0.0, 10.0);
   else if (peripheral)
      cutCentrality = new AliRsnCutValue("cutPeripheral", AliRsnValue::kEventCentralityV0, 60.0, 70.0);
   
   // primary vertex is always used
   task->GetEventCuts()->AddCut(cutVertex);
   
   // set cut scheme as AND of primary vertex and centrality, if initialized
   if (cutCentrality) {
      task->GetEventCuts()->AddCut(cutCentrality);
      task->GetEventCuts()->SetCutScheme(Form("%s & %s", cutVertex->GetName(), cutCentrality->GetName()));
   } else {
      task->GetEventCuts()->SetCutScheme(cutVertex->GetName());
   }
   ::Info("AddEventStuff", "Scheme for event cuts: %s", task->GetEventCuts()->GetCutScheme().Data());
   
   // ==================================================================================================================
   // == EVENT FUNCTIONS ===============================================================================================
   // ==================================================================================================================

   // we want to add an AliRsnFunction to compute multiplicity distribution
   // it is needed in order to know how many events we have in each multiplicity bin
   
   // axes
   AliRsnValue *axisEvMultSPD = new AliRsnValue("MultSPD", AliRsnValue::kEventMultSPD, 0.0, 150.0, 1.0);
   AliRsnValue *axisEvMultMC  = new AliRsnValue("MultMC" , AliRsnValue::kEventMultMC , 0.0, 150.0, 1.0);

   // create function and add axis
   AliRsnFunction *fcnEv = new AliRsnFunction;
   if (!fcnEv->AddAxis(axisEvMultSPD)) return kFALSE;
   if (isMC && !fcnEv->AddAxis(axisEvMultMC)) return kFALSE;

   // add functions to pairs
   task->GetInfo()->AddEventFunction(fcnEv);
   
   return kTRUE;
}
