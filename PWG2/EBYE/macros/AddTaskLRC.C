AliAnalysisTaskLRC *AddTaskLRC(double loForwardETA,double hiForwardETA,double loBackwardETA, double hiBackwardETA)
{
// This macro adds AliAnalysisTaskLRC to existing AnalysisManager
// Paramiters are : loForwardETA - lover ETA for Forward window 
// hiForwardETA - higer ETA for Forward window
// loBackwardETA, hiBackwardETA - same for Bakward window
// Ex: AddTaskLRC(-0.8, -0.6,0.6, 0.8); 



  
   // A. Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTaskLRC", "No analysis manager to connect to.");
      return NULL;
   }  

   // B. Check the analysis type using the event handlers connected to the analysis
   //    manager. The availability of MC handler cann also be checked here.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      Error("AddTaskLRC", "This task requires an input event handler");
      return NULL;
   }  
   TString type = mgr->GetInputEventHandler()->GetDataType(); 
   
   if(type!="ESD")
   {
   Error("AddTaskLRC","ESD is only supported now");
   return NULL;  
   
   }

   // Creating task and output container name 
   TString lF,hF,lB,hB;
   lF+=loForwardETA; lF.Remove(TString::kBoth,' ');
   hF+=hiForwardETA; hF.Remove(TString::kBoth,' ');
   lB+=loBackwardETA; lB.Remove(TString::kBoth,' ');
   hB+=hiBackwardETA; hB.Remove(TString::kBoth,' ');
   
   TString taskname="TaskLRCw"+lF+"to"+hF+"vs"+lB+"to"+hB;
    
   
   
   // C. Create the task, add it to manager.
   //===========================================================================
   AliAnalysisTaskLRC *taskLRC = new AliAnalysisTaskLRC(taskname);
   mgr->AddTask(taskLRC);

   // D. Configure the analysis task. Extra parameters can be used via optional
   // arguments of the AddTaskXXX() function.
   //===========================================================================
 
   taskLRC->SetETAWindows(loForwardETA,hiForwardETA,loBackwardETA,hiBackwardETA);
  
   // E. Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *cout_LRC = mgr->CreateContainer(taskname, TList::Class(),
                 AliAnalysisManager::kOutputContainer,"LRC.ESD.C.root");                             
   mgr->ConnectInput(taskLRC, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskLRC, 0, cout_LRC);

   // Return task pointer at the end
   return taskLRC;
};

TList *AddLRCTaskSet()
{// This routine uses AddTaskLRC to add a set of AliAnalysisTaskLRC 
//corresponding standary windows set  to existing analysis manager 

TList *TaskList = new TList();

//FB
TaskList->Add(AddTaskLRC(-0.2,-0.0,0.0,0.2));
TaskList->Add(AddTaskLRC(-0.4,-0.0,0.0,0.4));
TaskList->Add(AddTaskLRC(-0.6,-0.0,0.0,0.6));
TaskList->Add(AddTaskLRC(-0.8,-0.0,0.0,0.8));
//0.2 gap
TaskList->Add(AddTaskLRC(-0.4,-0.2,0.2,0.4));
TaskList->Add(AddTaskLRC(-0.6,-0.4,0.4,0.6));
TaskList->Add(AddTaskLRC(-0.8,-0.6,0.6,0.8));
//0.4 gap
TaskList->Add(AddTaskLRC(-0.6,-0.2,0.2,0.6));
TaskList->Add(AddTaskLRC(-0.8,-0.4,0.4,0.8));
//0.6 gap
TaskList->Add(AddTaskLRC(-0.8,-0.2,0.2,0.8));

//FULL
TaskList->Add(AddTaskLRC(-0.8,0.8,-0.8,0.8));

return TaskList;
}