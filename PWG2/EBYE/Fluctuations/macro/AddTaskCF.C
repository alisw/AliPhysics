/************************************************ 
 Charge Flatuation analysis task

 Auther: Satyajit Jena
 Email:  sjena@cern.ch
 Mon Oct 25 12:47:38 CEST 2010

*************************************************/

AliEbyECFAnalysisTask *AddTaskCF()
{
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTaskCF", "ERROR: No Analysis Manager");
      return NULL;
   }

   if (!mgr->GetInputEventHandler()) {
     Error("AddTaskCF", "ERROR: No input event handler");
     return NULL;
   }

   TString type = mgr->GetInputEventHandler()->GetDataType();
   TString OutName;
   OutName="CF."+type;
   TString outputFileName = AliAnalysisManager::GetCommonFileName(); 

   // getting default name

   outputFileName += ":PWG2EbyE_CF"; // adding directory type

   gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/macro/ConfigureCFAnalysis.C");

   AliEbyEChargeFluctuationAnalysis *analysis = GetAnalysisCFObject();

   AliEbyECFAnalysisTask *taskCF 
     = new AliEbyECFAnalysisTask("AliEbyECFAnalysisTask");
  
   taskCF->SetAnalysisObject(analysis);	
   
   mgr->AddTask(taskCF);
   
   AliAnalysisDataContainer *cout 
     = mgr->CreateContainer(OutName, TList::Class(), 
			    AliAnalysisManager::kOutputContainer,
			    outputFileName.Data());
   
   mgr->ConnectInput(taskCF, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskCF, 1, cout);
   return taskCF;

}

