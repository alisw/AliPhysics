#include "AliLog.h"
#include "AliAnalysisTaskPHOSCalibSelection.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

AliAnalysisTaskPHOSCalibSelection* AddTaskPHOSPi0Calibration(
    bool fSaveCells                     = true,
    bool fSaveFullTree                  = true,
    float minCellEnergy                 = 0.01,
    const char *trigSuffix              = "",
    const char * outputFileName               = ""
)
{

    TString wagon = trigSuffix;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskPHOSTriggerQA", "No analysis manager to connect to.");
        return NULL;
    }  
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskPHOSPi0Calibration", "This task requires an input event handler");
        return NULL;
    }

    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

    AliAnalysisDataContainer *couttree    = 0;

    AliAnalysisTaskPHOSCalibSelection * pi0calib = new AliAnalysisTaskPHOSCalibSelection("PHOSPi0Calibration");
    pi0calib->SetCellMinimumEnergy(minCellEnergy);
    pi0calib->SetSaveCells();
    pi0calib->SetSaveFullTree();

    TString containerName = Form("Pi0Calibration") ;


    couttree    = mgr->CreateContainer("Pi0Calibration", TTree::Class(), AliAnalysisManager::kOutputContainer, "Pi0CalibrationTree.root"); 
 
    mgr->AddTask(pi0calib);
    mgr->ConnectInput(pi0calib, 0, cinput1); 
    mgr->ConnectOutput(pi0calib, 1, couttree);
    
  return pi0calib;
 
}


