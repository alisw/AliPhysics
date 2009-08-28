AliAnalysisTaskLRC *AddTaskLRC(Bool_t RunKine=kFALSE)
{
// This macro adds AliAnalysisTaskLRC to existing AnalysisManager
// RunKine paramiter switch task to kinematics analysis 

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

// Version line : 3.5
// Version 3.5.5


 
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
   
   cout<<" # TaskLRC - input :"<<type<<"\n";
   
  

     
   // C. Create the task, add it to manager.
   //===========================================================================
   AliAnalysisTaskLRC *taskLRC = new AliAnalysisTaskLRC("Tak_LRC",RunKine);
   mgr->AddTask(taskLRC);

   taskLRC->fMinPtLimit=0.2;    // 200MeV minimal Pt (5GeV max Pt by default)
      
   // D. Configure the analysis task. Extra parameters can be used via optional
   // arguments of the AddTaskXXX() function.
   //===========================================================================
 
   //Adding LRC processors to the task
   
 //FB
taskLRC->AddLRCProcess(new AliLRCProcess(-0.2,-0.0,0.0,0.2));
taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,0.0,0.4));
taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.0,0.0,0.6));
taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.0,0.0,0.8));
//0.2 gap
taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.2,0.2,0.4));
taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.4,0.4,0.6));
taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.6,0.6,0.8));
//0.4 gap
taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.2,0.2,0.6));
taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4,0.4,0.8));
//0.6 gap
taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.2,0.2,0.8));

//FULL
taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,0.8,-0.8,0.8));

TString OutName;
OutName="LRC_out_"+type;

   
   // E. Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *cout_LRC = mgr->CreateContainer(OutName, TList::Class(),
                 AliAnalysisManager::kOutputContainer,"LRC.C.root");                             
   mgr->ConnectInput(taskLRC, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskLRC, 1, cout_LRC);

   // Return task pointer at the end
   return taskLRC;
}

