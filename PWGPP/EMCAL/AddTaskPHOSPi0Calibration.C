#include "AliLog.h"
#include "AliAnalysisTaskPHOSCalibSelection.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

AliAnalysisTaskPHOSPi0CalibSelection* AddTaskPHOSPi0Calibration(
    bool fSaveCells                     = true,
    bool fSaveFullTree                  = true,
    float minCellEnergy                 = 0.25,
    string outputFileName               = ""
){

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskEMCALTriggerQA", "No analysis manager to connect to.");
        return NULL;
    }  
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskEMCALPi0Calibration", "This task requires an input event handler");
        return NULL;
    }

    AliAnalysisTaskPHOSCalibSelection * pi0calib = new AliAnalysisTaskPHOSCalibSelection("PHOSPi0Calibration");
    pi0calib->SetCellMinEnergy(minCellEnergy);


    if( wagon.Length()==0 ){
    couthistos  = mgr->CreateContainer( "PHOSPi0Calibration",
                                        TList::Class(),
                                        AliAnalysisManager::kOutputContainer,
                                        outputFile.Data() );
    couttree    = mgr->CreateContainer("Pi0Calibration_Trig") ,
                                                                            TTree::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            ("Pi0CalibrationTree.root"); 
  } else {
    TString containerName = Form("Pi0Calibration") ;
    couthistos  = mgr->CreateContainer( wagon,
                                        TList::Class(),
                                        AliAnalysisManager::kOutputContainer,
    couttree    = mgr->CreateContainer( "Pi0Calibration_Trig"),
                                                                            TTree::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            ("Pi0CalibrationTree.root")); 
  }  
  
  mgr->AddTask(pi0calib);
  mgr->ConnectInput(pi0calib, 0, cinput1);                                
  mgr->ConnectOutput (pi0calib, 1, couthistos);
  if(fSaveCells || fSaveClusters) mgr->ConnectOutput( pi0calib, 2, couttree);
  

  return pi0calib;

    
}