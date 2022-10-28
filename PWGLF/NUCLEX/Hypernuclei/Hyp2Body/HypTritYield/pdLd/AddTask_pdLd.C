//#include "AliAnalysisTask_pdLd.h"
//#include "AliAnalysisTask_pdLd.cxx"

//AliAnalysisTask* AddTask_pdLd(TString name = "name") {
AliAnalysisTask* AddTask_pdLd(
  TString Name = "AnalysisTask_pdLd",
  int CollisionSystem = 1,
  const char *Variation = "0") {

  // CollisionSystem:
  //
  // 0 = pp collisions
  // 1 = Pb-Pb collisions (central: 0-10%)
  // 2 = Pb-Pb collisions (semi-central: 30-50%)

  TString suffix = TString::Format("%s",Variation);



  // collision system
  bool ispp   = false;
  bool ispPb  = false;
  bool isPbPb = false;


  // trigger
  bool isMinBias      = false;
  bool isHighMult     = false;
  bool isCentral      = false;
  bool isSemiCentral  = false;

  if(CollisionSystem == 0)
  {
    ispp    = true;
  }

  if(CollisionSystem == 1)
  {
    isPbPb    = true;
    isCentral = true;
  }


  if(CollisionSystem == 2)
  {
    isPbPb	  = true;
    isSemiCentral = true;
  }

  if(CollisionSystem > 2)
  {
    std::cout << "AddTask_pdLd: CollisionSystem: " << CollisionSystem << " is not defined!" << std::endl;
  }




  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if(!manager)
  {
    std::cout << "AddTask_pdLd: No AliAnalysisManager found" << std::endl;
    return nullptr;
  }

  if(!manager->GetInputEventHandler())
  {
    std::cout << "AddTask_pdLd: No InputEventHandler found" << std::endl;
    return nullptr;
  }

  AliAnalysisTask_pdLd* task = new AliAnalysisTask_pdLd(Name.Data(),CollisionSystem);   
  if(!task)
  {
    std::cout << "AddTask_pdLd: No AliAnalysisTask found" << std::endl;
    return nullptr;
  }

  // add the task to the manager
  manager->AddTask(task);

  AliAnalysisDataContainer *container = manager->GetCommonInputContainer();
  if(!container)
  {
    std::cout << "AddTask_pdLd: No CommonInputContainer found" << std::endl;
    return nullptr;
  }

  // connect the manager to the task
  manager->ConnectInput(task,0,container);
  


  // pp collisions
  if(ispp == true && ispPb == false && isPbPb == false)
    {

      task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kHighMultV0);

    } // end of ispp



  // p-Pb collisions
  if(ispp == false && ispPb == true && isPbPb == false)
    {

      task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kHighMultV0);

    } // end of ispPb



  // Pb-Pb collisions
  if(ispp == false && ispPb == false && isPbPb == true)
    {

     task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kHighMultV0 | AliVEvent::kSemiCentral | AliVEvent::kCentral);
     //task->SelectCollisionCandidates(AliVEvent::kSemiCentral);

    } // end of isPbPb


  TString TaskName = TString::Format("%s:AnalysisTask_pdLd%s",AliAnalysisManager::GetCommonFileName(),suffix.Data());
  TString EventListName		= "Events";
  TString ProtonListName	= "Protons";
  TString DeuteronListName	= "Deuterons";
  TString pdListName		= "Proton-Deuteron-Pairs";
  TString AntiProtonListName	= "AntiProtons";
  TString AntiDeuteronListName	= "AntiDeuterons";
  TString apadListName		= "AntiProton-AntiDeuteron-Pairs";


  manager->ConnectOutput(task,1,manager->CreateContainer(EventListName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,TaskName.Data()));
  manager->ConnectOutput(task,2,manager->CreateContainer(ProtonListName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,TaskName.Data()));
  manager->ConnectOutput(task,3,manager->CreateContainer(DeuteronListName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,TaskName.Data()));
  manager->ConnectOutput(task,4,manager->CreateContainer(pdListName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,TaskName.Data()));
  manager->ConnectOutput(task,5,manager->CreateContainer(AntiProtonListName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,TaskName.Data()));
  manager->ConnectOutput(task,6,manager->CreateContainer(AntiDeuteronListName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,TaskName.Data()));
  manager->ConnectOutput(task,7,manager->CreateContainer(apadListName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,TaskName.Data()));



  // return a pointer to the task
  return task;

}
