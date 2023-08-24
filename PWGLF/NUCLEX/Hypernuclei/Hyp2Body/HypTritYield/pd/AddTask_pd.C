#include "AliAnalysisTask_pd.cxx"



//AliAnalysisTask* AddTask_pd(TString name = "name") {
AliAnalysisTask_pd* AddTask_pd(
  TString Name = "AnalysisTask_pd",
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
    std::cout << "AddTask_pd: x-x-x-x-> No AliAnalysisManager found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pd: AliAnalysisManager received" << std::endl;




  if(!manager->GetInputEventHandler())
  {
    std::cout << "AddTask_pd: x-x-x-x-> No InputEventHandler found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pd: InputEventHandler received" << std::endl;




  AliAnalysisTask_pd* task = new AliAnalysisTask_pd(Name.Data(),CollisionSystem);   
  if(!task)
  {
    std::cout << "AddTask_pd: x-x-x-x-> No AliAnalysisTask found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pd: AliAnalysisTask created" << std::endl;



  // add the task to the manager
  manager->AddTask(task);
  if(DebugAddTask) std::cout << "AddTask_pd: task added to AliAnalysisManager" << std::endl;



  AliAnalysisDataContainer *container = manager->GetCommonInputContainer();
  if(!container)
  {
    std::cout << "AddTask_pd: x-x-x-x-> No CommonInputContainer found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pd: CommonInputContainer received" << std::endl;



  // connect the manager to the task
  manager->ConnectInput(task,0,container);
  if(DebugAddTask) std::cout << "AddTask_pd: Input container connected" << std::endl;
  
  

  task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kHighMultV0 | AliVEvent::kSemiCentral | AliVEvent::kCentral);
  if(DebugAddTask) std::cout << "AddTask_pd: Collision candidates selected" << std::endl;


  TString TaskName = TString::Format("%s:AnalysisTask_pd%s",AliAnalysisManager::GetCommonFileName(),suffix.Data());


  TString ProtonTreeName = Form("ProtonTree%s",suffix.Data());
  AliAnalysisDataContainer *ProtonTreeContainer = manager->CreateContainer(
    ProtonTreeName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pd: ProtonContainer created" << std::endl;


  TString DeuteronTreeName = Form("DeuteronTree%s",suffix.Data());
  AliAnalysisDataContainer *DeuteronTreeContainer = manager->CreateContainer(
    DeuteronTreeName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pd: DeuteronContainer created" << std::endl;


  TString AntiProtonTreeName = Form("AntiProtonTree%s",suffix.Data());
  AliAnalysisDataContainer *AntiProtonTreeContainer = manager->CreateContainer(
    AntiProtonTreeName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pd: AntiProtonContainer created" << std::endl;


  TString AntiDeuteronTreeName = Form("AntiDeuteronTree%s",suffix.Data());
  AliAnalysisDataContainer *AntiDeuteronTreeContainer = manager->CreateContainer(
    AntiDeuteronTreeName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pd: AntiDeuteronContainer created" << std::endl;




  manager->ConnectOutput(task,1,ProtonTreeContainer);
  if(DebugAddTask) std::cout << "AddTask_pd: Output container 1 connected" << std::endl;

  manager->ConnectOutput(task,2,DeuteronTreeContainer);
  if(DebugAddTask) std::cout << "AddTask_pd: Output container 2 connected" << std::endl;

  manager->ConnectOutput(task,3,AntiProtonTreeContainer);
  if(DebugAddTask) std::cout << "AddTask_pd: Output container 3 connected" << std::endl;

  manager->ConnectOutput(task,4,AntiDeuteronTreeContainer);
  if(DebugAddTask) std::cout << "AddTask_pd: Output container 4 connected" << std::endl;


  // return a pointer to the task
  return task;

}
