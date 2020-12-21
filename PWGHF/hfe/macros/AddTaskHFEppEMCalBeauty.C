
///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskHFEppEMCalBeauty                                       //
// Author: Vivek Singh                                           //
//                                                               //
///////////////////////////////////////////////////////////////////

class AliAnalysisDataContainer;

AliAnalysisHFEppEMCalBeauty* AddTaskHFEppEMCalBeauty(

TString name = "", 
Bool_t isMC=kFALSE,
AliVEvent::EOfflineTriggerTypes trigger=AliVEvent::kINT7 ,
Bool_t isEG1=kFALSE,

//Bool_t PhysSelINT7 = kTRUE,
Bool_t useTender = kTRUE,
Bool_t ClsTypeEMC = kTRUE,
Bool_t ClsTypeDCAL = kTRUE

)

{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { ::Error("AliAnalysisHFEppEMCalBeauty", "No analysis manager to connect to.");
    return NULL;
    }

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file

    TString finDirname         = "_";
    TString outBasicname = "MyOutputContainer";
    TString profname       = "coutputProf";     
    
    finDirname        += name.Data();
    outBasicname      += finDirname.Data();

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    TString taskname="ElecAnalysis";
    AliAnalysisHFEppEMCalBeauty *HFeTask = new AliAnalysisHFEppEMCalBeauty(taskname.Data());
    //HFeTask->SetDebugLevel(2);

    HFeTask->SelectCollisionCandidates(trigger);
    HFeTask->SetMCAnalysis(isMC);
    HFeTask->SetTrigger(trigger);
    HFeTask->SetTenderSwitch(useTender);


  if(trigger==AliVEvent::kINT7){

    isEG1=kFALSE;
    HFeTask->SetEMCalTriggerEG1(kFALSE);
    HFeTask->SetEMCalTriggerDG1(kFALSE);
    HFeTask->SetEMCalTriggerEG2(kFALSE);
    HFeTask->SetEMCalTriggerDG2(kFALSE);
    HFeTask->SetClusterTypeEMC(ClsTypeEMC);
    HFeTask->SetClusterTypeDCAL(ClsTypeDCAL);
    //cout<<"1 trigger  "<<trigger<<"   "<< isEG1 <<endl;   	
    }
  
  if(trigger==AliVEvent::kEMCEGA ){

    HFeTask->SetEMCalTriggerEG1(kTRUE);
    HFeTask->SetEMCalTriggerDG1(kTRUE);
    //HFeTask->SetEMCalTriggerEG2(kFALSE);
    //HFeTask->SetEMCalTriggerDG2(kFALSE);
    HFeTask->SetClusterTypeEMC(ClsTypeEMC);
    HFeTask->SetClusterTypeDCAL(ClsTypeDCAL);
    //cout<<" 2 trigger  "<<trigger<<"   "<< isEG1 <<endl;
    }
 
 if(trigger==AliVEvent::kEMCEGA && isEG1==kFALSE){
        //HFeTask->SetEMCalTriggerEG1(kFALSE);
        //HFeTask->SetEMCalTriggerDG1(kFALSE);
        HFeTask->SetEMCalTriggerEG2(kTRUE);
        HFeTask->SetEMCalTriggerDG2(kTRUE);
        HFeTask->SetClusterTypeEMC(ClsTypeEMC);
        HFeTask->SetClusterTypeDCAL(ClsTypeDCAL);
        //cout<<"3 trigger  "<<trigger<<"   "<< isEG1 <<endl;    
    }

    mgr->AddTask(HFeTask);    

    //_________Structure of Task O/P
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outBasicname,TList::Class(),AliAnalysisManager::kOutputContainer, fileName.Data());

    // your HFeTask needs input: here we connect the manager to your HFeTask
    mgr->ConnectInput(HFeTask,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(HFeTask,1,coutput1);
  //  mgr->ConnectOutput(HFetaskINT7,2,mgr->CreateContainer(profname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    


  return HFeTask;
}
