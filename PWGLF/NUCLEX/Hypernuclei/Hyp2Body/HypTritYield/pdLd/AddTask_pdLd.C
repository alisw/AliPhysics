//#include "AliAnalysisTask_pdLd.h"
//#include "AliAnalysisTask_pdLd.cxx"



//AliAnalysisTask* AddTask_pdLd(TString name = "name") {
AliAnalysisTask_pdLd* AddTask_pdLd(
  TString Name = "AnalysisTask_pdLd",
  int CollisionSystem = 0,
  const char *Variation = "1") {


  bool DebugAddTask = true;

  // CollisionSystem:
  //
  // 0 = pp collisions
  // 1 = Pb-Pb collisions (central: 0-10%)
  // 2 = Pb-Pb collisions (semi-central: 30-50%)

  TString suffix = TString::Format("%s",Variation);
//  if(DebugAddTask) std::cout << "suffix: " << suffix.Data() << std::endl;



  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if(!manager)
  {
    std::cout << "AddTask_pdLd: x-x-x-x-> No AliAnalysisManager found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pdLd: AliAnalysisManager received" << std::endl;




  if(!manager->GetInputEventHandler())
  {
    std::cout << "AddTask_pdLd: x-x-x-x-> No InputEventHandler found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pdLd: InputEventHandler received" << std::endl;




  AliAnalysisTask_pdLd* task = new AliAnalysisTask_pdLd(Name.Data(),CollisionSystem);   
  if(!task)
  {
    std::cout << "AddTask_pdLd: x-x-x-x-> No AliAnalysisTask found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pdLd: AliAnalysisTask created" << std::endl;



  // add the task to the manager
  manager->AddTask(task);
  if(DebugAddTask) std::cout << "AddTask_pdLd: task added to AliAnalysisManager" << std::endl;



  AliAnalysisDataContainer *container = manager->GetCommonInputContainer();
  if(!container)
  {
    std::cout << "AddTask_pdLd: x-x-x-x-> No CommonInputContainer found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pdLd: CommonInputContainer received" << std::endl;



  // connect the manager to the task
  manager->ConnectInput(task,0,container);
  if(DebugAddTask) std::cout << "AddTask_pdLd: Input container connected" << std::endl;
  
  

  task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kHighMultV0 | AliVEvent::kSemiCentral | AliVEvent::kCentral);
  if(DebugAddTask) std::cout << "AddTask_pdLd: Collision candidates selected" << std::endl;


  TString TaskName = TString::Format("%s:AnalysisTask_pdLd%s",AliAnalysisManager::GetCommonFileName(),suffix.Data());

  TString EventListName	= Form("Events%s",suffix.Data());
  AliAnalysisDataContainer *EventContainer = manager->CreateContainer(
    EventListName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pdLd: EventContainer created" << std::endl;


  TString ProtonListName = Form("Protons%s",suffix.Data());
  AliAnalysisDataContainer *ProtonContainer = manager->CreateContainer(
    ProtonListName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pdLd: ProtonContainer created" << std::endl;


  TString DeuteronListName = Form("Deuterons%s",suffix.Data());
  AliAnalysisDataContainer *DeuteronContainer = manager->CreateContainer(
    DeuteronListName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pdLd: DeuteronContainer created" << std::endl;


  TString ProtonDeuteronListName = Form("Proton-Deuteron-Pairs%s",suffix.Data());
  AliAnalysisDataContainer *ProtonDeuteronContainer = manager->CreateContainer(
    ProtonDeuteronListName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pdLd: ProtonDeuteronContainer created" << std::endl;


  TString AntiProtonListName = Form("AntiProtons%s",suffix.Data());
  AliAnalysisDataContainer *AntiProtonContainer = manager->CreateContainer(
    AntiProtonListName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pdLd: AntiProtonContainer created" << std::endl;


  TString AntiDeuteronListName = Form("AntiDeuterons%s",suffix.Data());
  AliAnalysisDataContainer *AntiDeuteronContainer = manager->CreateContainer(
    AntiDeuteronListName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pdLd: AntiDeuteronContainer created" << std::endl;


  TString AntiProtonAntiDeuteronListName = Form("AntiProton-AntiDeuteron-Pairs%s",suffix.Data());
  AliAnalysisDataContainer *AntiProtonAntiDeuteronContainer = manager->CreateContainer(
    AntiProtonAntiDeuteronListName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pdLd: AntiProtonAntiDeuteronContainer created" << std::endl;





  manager->ConnectOutput(task,1,EventContainer);
  if(DebugAddTask) std::cout << "AddTask_pdLd: Output container 1 connected" << std::endl;

  manager->ConnectOutput(task,2,ProtonContainer);
  if(DebugAddTask) std::cout << "AddTask_pdLd: Output container 2 connected" << std::endl;

  manager->ConnectOutput(task,3,DeuteronContainer);
  if(DebugAddTask) std::cout << "AddTask_pdLd: Output container 3 connected" << std::endl;

  manager->ConnectOutput(task,4,ProtonDeuteronContainer);
  if(DebugAddTask) std::cout << "AddTask_pdLd: Output container 4 connected" << std::endl;

  manager->ConnectOutput(task,5,AntiProtonContainer);
  if(DebugAddTask) std::cout << "AddTask_pdLd: Output container 5 connected" << std::endl;

  manager->ConnectOutput(task,6,AntiDeuteronContainer);
  if(DebugAddTask) std::cout << "AddTask_pdLd: Output container 6 connected" << std::endl;

  manager->ConnectOutput(task,7,AntiProtonAntiDeuteronContainer);
  if(DebugAddTask) std::cout << "AddTask_pdLd: Output container 7 connected" << std::endl;


  // return a pointer to the task
  return task;

}
