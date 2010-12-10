/************************************************
 Multiplisity Flatuation analysis task

 Auther: Satyajit Jena
 Email:  sjena@cern.ch
 Mon Oct 25 12:47:38 CEST 2010

*************************************************/

AliEbyEMFAnalysisTask *AddTaskMF()
{
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTaskMF", "ERROR: No Analysis Manager");
      return NULL;
   }

   if (!mgr->GetInputEventHandler()) {
     Error("AddTaskMF", "ERROR: No input event handler");
     return NULL;
   }

   TString type = mgr->GetInputEventHandler()->GetDataType();
   TString OutName;
   OutName="MF."+type;
   TString outputFileName = AliAnalysisManager::GetCommonFileName(); 

   // getting default name

   outputFileName += ":PWG2EbyE_MF"; // adding directory type

   gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/macro/ConfigureMFAnalysis.C");

   AliEbyEMultiplicityFluctuationAnalysis *analysis = GetAnalysisMFObject();

   AliEbyEMFAnalysisTask *taskMF 
     = new AliEbyEMFAnalysisTask("AliEbyEMFAnalysisTask");
  
   taskMF->SetAnalysisObject(analysis);	
   
   mgr->AddTask(taskMF);
   
   AliAnalysisDataContainer *cout 
     = mgr->CreateContainer(OutName, TList::Class(), 
			    AliAnalysisManager::kOutputContainer,
			    outputFileName.Data());
   
   mgr->ConnectInput(taskMF, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskMF, 1, cout);
   return taskMF;

}

