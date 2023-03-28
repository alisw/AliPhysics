

AliAnalysisTask_pp_CreateTrees_PairsOnly* AddTask_pp_CreateTrees_PairsOnly(
  TString Name = "AnalysisTask_pp_CreateTrees_PairsOnly",
  int CollisionSystem = 1,
  const char *Variation = "0",
  bool UseOpenCuts = true) {

  bool DebugAddTask = true;

  // CollisionSystem:
  //
  // 1 = Pb-Pb collisions (central: 0-10%)
  // 2 = Pb-Pb collisions (semi-central: 30-50%)
  // 3 = pp collisions (MetaLHC16)
  // 4 = pp collisions (MetaLHC17)
  // 5 = pp collisions (MetaLHC18)



  TString suffix = TString::Format("%s",Variation);
  TString prefix = TString::Format("%i",CollisionSystem);



  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if(!manager)
  {
    std::cout << "AddTask_pp_CreateTrees_PairsOnly: x-x-x-x-> No AliAnalysisManager found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: AliAnalysisManager received" << std::endl;




  if(!manager->GetInputEventHandler())
  {
    std::cout << "AddTask_pp_CreateTrees_PairsOnly: x-x-x-x-> No InputEventHandler found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: InputEventHandler received" << std::endl;




  AliAnalysisTask_pp_CreateTrees_PairsOnly* task = new AliAnalysisTask_pp_CreateTrees_PairsOnly(Name.Data(),CollisionSystem,UseOpenCuts);   
  if(!task)
  {
    std::cout << "AddTask_pp_CreateTrees_PairsOnly: x-x-x-x-> No AliAnalysisTask found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: AliAnalysisTask created" << std::endl;



  // add the task to the manager
  manager->AddTask(task);
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: task added to AliAnalysisManager" << std::endl;



  AliAnalysisDataContainer *container = manager->GetCommonInputContainer();
  if(!container)
  {
    std::cout << "AddTask_pp_CreateTrees_PairsOnly: x-x-x-x-> No CommonInputContainer found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: CommonInputContainer received" << std::endl;



  // connect the manager to the task
  manager->ConnectInput(task,0,container);
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: Input container connected" << std::endl;
 

  if((CollisionSystem == 1) || (CollisionSystem == 2)){

    task->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kCentral);
    std::cout << "AddTask_pp_CreateTrees_PairsOnly: SelectCollisionCandidates(AliVEvent::kSemiCentral || AliVEvent::kCentral)" << std::endl;

  }


  if((CollisionSystem == 3) || (CollisionSystem == 4) || CollisionSystem == 5){
 
    task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kHighMultV0);
    std::cout << "AddTask_pp_CreateTrees_PairsOnly: SelectCollisionCandidates(AliVEvent::kHighMultV0)" << std::endl;

  }



  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: Collision candidates selected" << std::endl;


  TString TaskName = TString::Format("%s:%sAnalysisTask_pp_CreateTrees_PairsOnly%s",AliAnalysisManager::GetCommonFileName(),prefix.Data(),suffix.Data());


  TString ProtonTreeName = Form("%sProtonTree%s",prefix.Data(),suffix.Data());
  AliAnalysisDataContainer *ProtonTreeContainer = manager->CreateContainer(
    ProtonTreeName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: ProtonContainer created" << std::endl;



  TString AntiProtonTreeName = Form("%sAntiProtonTree%s",prefix.Data(),suffix.Data());
  AliAnalysisDataContainer *AntiProtonTreeContainer = manager->CreateContainer(
    AntiProtonTreeName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: AntiProtonContainer created" << std::endl;


  TString HistogramName = Form("%sHistograms%s",prefix.Data(),suffix.Data());
  AliAnalysisDataContainer *HistogramContainer = manager->CreateContainer(
    HistogramName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: HistogramContainer created" << std::endl;



  manager->ConnectOutput(task,1,ProtonTreeContainer);
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: Output container 1 connected" << std::endl;

  manager->ConnectOutput(task,2,AntiProtonTreeContainer);
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: Output container 2 connected" << std::endl;

  manager->ConnectOutput(task,3,HistogramContainer);
  if(DebugAddTask) std::cout << "AddTask_pp_CreateTrees_PairsOnly: Output container 3 connected" << std::endl;


  // return a pointer to the task
  return task;

}
