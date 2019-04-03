#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGGA/PHOSTasks/CaloCellQA/phys/macros/AddTaskCaloCellsPhysQA.C>
#endif

AliAnalysisTaskSE * AddTaskPhysPHOSQA()
{
  AliAnalysisTaskCaloCellsPhysQA *taskPHOSCellQA = dynamic_cast<AliAnalysisTaskCaloCellsPhysQA *> (
    AddTaskCaloCellsPhysQA(5, 1, "CaloCellsQA.root","PHOSCellsQA")
  );
  taskPHOSCellQA->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0); 

  AliCaloCellsPhysQA * qa = (AliCaloCellsPhysQA *) taskPHOSCellQA->GetCaloCellsQA();
  qa->SetPhysicsClusterCut(0.3, 3, 12.5e-9); 
  return taskPHOSCellQA;
}
