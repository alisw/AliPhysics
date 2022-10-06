/// \file AddTaskEMCALPi0Calibration.C
/// \ingroup EMCALPerfAddTaskMacros
/// \brief Configuration of task AliAnalysisTaskEMCALPi0CalibSelection.
///
/// Configuration of task AliAnalysisTaskEMCALPi0CalibSelection, which fills invariant mass
/// histograms for each of the EMCal channels. It has to be executed in several iterations.
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

// Root6
#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliAnalysisTaskEMCALPi0CalibSelectionV2.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)

#endif // CINT

/// Main method
///
/// The parameters for the analysis are:
/// \param calibPath : TString with full path and name of file with calibration factors from previous iteration.
/// \param trigger   : TString, event that triggered must contain this string. Leave for backward compatibility with old wagons
/// \param simu      : Int, 0: data, 1: MC, 2: JJ MC
/// \param fSaveCells: Bool, save the tree output for cells
/// \param fSaveClusters : Bool, save the tree output for clusters
/// \param fSaveFullTree: Bool, save the cell and cluster energy as float instead of short
/// \param minClusterEnergy : Double, minimum cluster energy used for pairing (GeV)
/// \param maxClusterEnergy : Double, maximum cluster energy used for pairing (GeV)
/// \param minNCells : Int, minimum number of cells for clusters
/// \param minClusterTime : Float, minimum cluster time (ns)
/// \param maxClusterTime : Float, maximum cluster time (ns)
/// \param maxDiffTimeClusterPair : Float, maximum cluster time difference between paired clusters (ns)
/// \param bSameSM : Bool, switch on/off paring only in same SM
/// \param outputFile: TString with name of output file (AnalysisResults.root).
/// \param trigSuffix :  A string with the trigger class, abbreviated, to run multiple triggers in same train
///
AliAnalysisTaskEMCALPi0CalibSelectionV2 * AddTaskEMCALPi0CalibrationV2(
  TString calibPath              = "", // "alienpath/RecalibrationFactors.root"
  TString trigger                = "",
  Int_t  simu                    = 0,
  Bool_t  fSaveHistos            = kTRUE,
  Bool_t  fSaveCells             = kFALSE,
  Bool_t  fSaveClusters          = kFALSE,
  Bool_t  fSaveFullTree          = kFALSE,
  Bool_t  isHeavyIon             = kFALSE,
  Bool_t  fNContCut              = kFALSE,
  Float_t minCellEnergy          = 0.05,
  Double_t minClusterEnergy      = 0.7,
  Double_t maxClusterEnergy      = 10,
  Int_t minNCells                = 2,
  Float_t minClusterTime         = 300,  // 560 ns in Run1
  Float_t maxClusterTime         = 800,  // 610 ns in Run1
  Float_t maxDiffTimeClusterPair = 100,  //  20 ns in Run1  
  Bool_t bSameSM                 = kTRUE,
  TString corrTaskSetting        = "",
  TString outputFile             = "", // AnalysisResults.root
  TString periodName             = "",
  const char *trigSuffix         = ""
) 
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
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
  
  TString wagon = trigSuffix;
  if ( wagon.Length() > 0 ) trigger = wagon;
  wagon+="_";
  wagon+=corrTaskSetting.Data();
  
  AliAnalysisTaskEMCALPi0CalibSelectionV2 * pi0calib = new AliAnalysisTaskEMCALPi0CalibSelectionV2(Form("EMCALPi0Calibration_%s_%s",trigger.Data(),corrTaskSetting.Data()));
  pi0calib->SetCellMinEnergy(minCellEnergy);
  pi0calib->SetClusterMinEnergy(minClusterEnergy);
  pi0calib->SetClusterMaxEnergy(maxClusterEnergy);
  pi0calib->SetClusterLambda0Cuts(0.1,0.5);
  
  pi0calib->SetAsymmetryCut(1.);
  pi0calib->SetClusterMinNCells(minNCells);
  pi0calib->SetNCellsGroup(0);
  if (bSameSM) pi0calib->SwitchOnSameSM();
    else pi0calib->SwitchOffSameSM();

  pi0calib->SetPeriodName(periodName);
  
  // Timing cuts
  pi0calib->SetPairDTimeCut(maxDiffTimeClusterPair);   
  pi0calib->SetClusterMinTime(minClusterTime); 
  pi0calib->SetClusterMaxTime(maxClusterTime); 

  pi0calib->SetTriggerName(trigger);
    
  //---------------------
  // Geometry alignment
  //---------------------
  
  pi0calib->SetGeometryName("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  pi0calib->SwitchOnLoadOwnGeometryMatrices();
  // pi0calib->SetOADBFilePath("/home/ewa/alice/AliPhysics/OADB/OADB/EMCAL");
  pi0calib->SetCorrectionTaskSetting(corrTaskSetting);

  if( simu == 1 ) {
    pi0calib->SetIsMC();
  } else if (simu == 2){
    pi0calib->SetJJMC();
  }

  if( isHeavyIon ) {
    pi0calib->SetHeavyIon();
  }

  if( fSaveHistos ){
    pi0calib->SetSaveHistos();
  }

  if( fSaveCells ){
    pi0calib->SetSaveCells();
  }

  if( fSaveClusters ){
    pi0calib->SetSaveClusters();
  }

  if( fSaveFullTree ){
    pi0calib->SetSaveFullTree();
  }

  if( fNContCut ){
    pi0calib->SetNContributorsCut();
  }
  
  //---------------------
  // Print info
  //---------------------
  
  pi0calib->PrintInfo();
  
  mgr->AddTask(pi0calib);
  
  if(outputFile.Length()==0) outputFile = AliAnalysisManager::GetCommonFileName(); 

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
                                                                                                                                                              
  AliAnalysisDataContainer *couthistos  = 0;
  AliAnalysisDataContainer *couttree    = 0;
  if( wagon.Length()==0 ){
    couthistos  = mgr->CreateContainer( Form( "Pi0Calibration_Trig%s",trigger.Data() ),
                                        TList::Class(),
                                        AliAnalysisManager::kOutputContainer,
                                        outputFile.Data() );
    if( fSaveCells || fSaveClusters)    couttree    = mgr->CreateContainer( Form("Pi0Calibration_Trig%s_%s",trigger.Data(),corrTaskSetting.Data()) ,
                                                                            TTree::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            Form("Pi0CalibrationTree_%s.root",corrTaskSetting.Data())); 
  } else {
    TString containerName = Form("Pi0Calibration_%s",corrTaskSetting.Data()) ;
    couthistos  = mgr->CreateContainer( wagon,
                                        TList::Class(),
                                        AliAnalysisManager::kOutputContainer,
                                        Form("%s:%s",outputFile.Data(),containerName.Data()) );
    if( (fSaveCells || fSaveClusters) )   couttree    = mgr->CreateContainer( Form("Pi0Calibration_Trig%s_%s",trigger.Data(), corrTaskSetting.Data()),
                                                                            TTree::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            Form("Pi0CalibrationTree_%s.root",corrTaskSetting.Data())); 
  }  
  
  mgr->AddTask(pi0calib);
  mgr->ConnectInput(pi0calib, 0, cinput1);                                
  mgr->ConnectOutput (pi0calib, 1, couthistos);
  if(fSaveCells || fSaveClusters) mgr->ConnectOutput( pi0calib, 2, couttree);

  return pi0calib;
}
