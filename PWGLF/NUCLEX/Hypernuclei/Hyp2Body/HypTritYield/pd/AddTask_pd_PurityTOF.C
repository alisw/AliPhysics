
AliAnalysisTask_pd_PurityTOF* AddTask_pd_PurityTOF(
  TString Name = "AnalysisTask_pd_PurityTOF",
  Int_t CollisionSystem = 1,
  const char *Variation = "0",
  Bool_t UseOpenCuts = true,
  Bool_t isMC = false,
  Bool_t SavePairsOnly = true) {

  Bool_t DebugAddTask = true;

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
    std::cout << "AddTask_pd_PurityTOF: x-x-x-x-> No AliAnalysisManager found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pd_PurityTOF: AliAnalysisManager received" << std::endl;




  if(!manager->GetInputEventHandler())
  {
    std::cout << "AddTask_pd_PurityTOF: x-x-x-x-> No InputEventHandler found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pd_PurityTOF: InputEventHandler received" << std::endl;




  AliAnalysisTask_pd_PurityTOF* task = new AliAnalysisTask_pd_PurityTOF(Name.Data(),CollisionSystem,UseOpenCuts,isMC,SavePairsOnly);   
  if(!task)
  {
    std::cout << "AddTask_pd_PurityTOF: x-x-x-x-> No AliAnalysisTask found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pd_PurityTOF: AliAnalysisTask created" << std::endl;



  // add the task to the manager
  manager->AddTask(task);
  if(DebugAddTask) std::cout << "AddTask_pd_PurityTOF: task added to AliAnalysisManager" << std::endl;



  AliAnalysisDataContainer *container = manager->GetCommonInputContainer();
  if(!container)
  {
    std::cout << "AddTask_pd_PurityTOF: x-x-x-x-> No CommonInputContainer found" << std::endl;
    return nullptr;
  }
  if(DebugAddTask) std::cout << "AddTask_pd_PurityTOF: CommonInputContainer received" << std::endl;



  // connect the manager to the task
  manager->ConnectInput(task,0,container);
  if(DebugAddTask) std::cout << "AddTask_pd_PurityTOF: Input container connected" << std::endl;
 

  if(CollisionSystem == 1){

    task->SelectCollisionCandidates(AliVEvent::kCentral);
    std::cout << "AddTask_pd_PurityTOF: SelectCollisionCandidates(AliVEvent::kCentral)" << std::endl;

  }


  if(CollisionSystem == 2){

    task->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    std::cout << "AddTask_pd_PurityTOF: SelectCollisionCandidates(AliVEvent::kSemiCentral)" << std::endl;

  }


  if((CollisionSystem == 3) || (CollisionSystem == 4) || CollisionSystem == 5){
 
    task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kHighMultV0);
    std::cout << "AddTask_pd_PurityTOF: SelectCollisionCandidates(AliVEvent::kHighMultV0)" << std::endl;

  }



  if(DebugAddTask) std::cout << "AddTask_pd_PurityTOF: Collision candidates selected" << std::endl;


  TString TaskName = TString::Format("%s:%sAnalysisTask_pd_PurityTOF%s",AliAnalysisManager::GetCommonFileName(),prefix.Data(),suffix.Data());


  TString OutputListName = Form("%sOutputList%s",prefix.Data(),suffix.Data());
  AliAnalysisDataContainer *OutputListContainer = manager->CreateContainer(
    OutputListName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    TaskName.Data()
  );
  if(DebugAddTask) std::cout << "AddTask_pd_PurityTOF: OutputListContainer created" << std::endl;



  manager->ConnectOutput(task,1,OutputListContainer);
  if(DebugAddTask) std::cout << "AddTask_pd_PurityTOF: Output container 1 connected" << std::endl;



  // return a pointer to the task
  return task;

}
