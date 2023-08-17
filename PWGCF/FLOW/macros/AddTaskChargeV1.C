#include "TString.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskChargeV1.h"

AliAnalysisTaskChargeV1* AddTaskChargeV1(
    bool           doNUE=true,
    bool           doNUA=true,
    bool          ZDCcali=true
)
{
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskChargeV1* task = new AliAnalysisTaskChargeV1("MyTask");   //
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetNUEOn(doNUE);
    task->SetNUAOn(doNUA); 
    
    //=========================================================================
    //Read in Files
    TFile* fNUEFile = nullptr;
    TFile* fNUAFile = nullptr;
    TList* fListNUE = nullptr;
    TList* fListNUA = nullptr;
    TFile* fZDCCalibFile = nullptr;
    TList* fZDCCalibList = nullptr;

    if (!gGrid) TGrid::Connect("alien://");
    if (doNUE) {
   
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/j/jwan/CalibFile/18q/efficiencyBothpol18qnew.root","READ");
    //  fNUEFile = TFile::Open("./efficiencyBothpol18qnew.root","READ");
      fListNUE = dynamic_cast <TList*> (fNUEFile->Get("fMcEffiHij"));
   
  
    if(fListNUE) {
      task->SetListForNUE(fListNUE);
      std::cout<<"================  NUE List Set ================="<<std::endl;
    } else std::cout<<"!!!!!!!!!!!!!!!NUE List not Found!!!!!!!!!!!!!!!"<<std::endl;
  }
    if (doNUA) {
     
  
      fNUAFile = TFile::Open("alien:///alice/cern.ch/user/j/jwan/CalibFile/18q/WgtsNUAChargeAndPion_LHC18qPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root","READ");
     // fNUAFile = TFile::Open("./WgtsNUAChargeAndPion_LHC18qPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root","READ");
      fListNUA = dynamic_cast <TList*> (fNUAFile->Get("fNUA_ChPosChNeg"));
    
   
    if(fListNUA) {
      task->SetListForNUA(fListNUA);
      std::cout<<"================  NUA List Set ================="<<std::endl;
    } else std::cout<<"!!!!!!!!!!!!!!!NUA List not Found!!!!!!!!!!!!!!!"<<std::endl;
  }

  if(ZDCcali){
      fZDCCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/RecenteringResultFinal_2018q.root", "READ");
  //  fZDCCalibFile = TFile::Open("./RecenteringResultFinal_2018q.root");
    fZDCCalibList   = dynamic_cast <TList*> (fZDCCalibFile->Get("fOutputRecenter"));
    task->SetListForZDCCalib(fZDCCalibList);
    std::cout<<"================  ZDC List Set ================="<<std::endl;
  } else std::cout<<"!!!!!!!!!!!!!!!ZDC List not Found!!!!!!!!!!!!!!!"<<std::endl;

  


    // add your task to the managercd
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  //  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    AliAnalysisDataContainer* coutput_QA = mgr->CreateContainer("ListQA", TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           fileName.Data());
    AliAnalysisDataContainer* coutput_result = mgr->CreateContainer("ListResults", TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           fileName.Data());
    mgr->ConnectInput(task,0,cinput);
    mgr->ConnectOutput(task,1,coutput_result); 
    mgr->ConnectOutput(task,2,coutput_QA); 
   
       

  //  mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
