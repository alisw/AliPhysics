/// \file AddTaskEMCALTimeCalibration.C
/// \ingroup EMCALPerformanceMacros
/// \brief Configuration of task AliAnalysisTaskEMCALTimeCalib.
///
/// Configuration of task AliAnalysisTaskEMCALTimeCalib
///
/// The parameters for the analysis are:
/// \param outputFile: TString, output file name
/// \param geometryName: TString, geometry name
/// \param minClusterEne: Double_t, minimum cluster energy
/// \param maxClusterEne: Double_t, maximum cluster energy
/// \param minNcells: Int_t, minimum number of cells in cluster
/// \param maxNcells: Int_t, maximum number of cells in cluster
/// \param minLambda0LG: Double_t, minimum cluster lambda0 for cluster with low gain cell
/// \param maxLambda0LG: Double_t, maximum cluster lambda0 for cluster with low gain cell
/// \param minLambda0: Double_t, minimum cluster lambda0
/// \param maxLambda0: Double_t, maximum cluster lambda0
/// \param maxRtrack: Double_t, maximum cluster track distance
/// \param minCellEne: Double_t, minimum cell energy
/// \param referenceFileName: TString, name of reference file
/// \param referenceSMFileName: TString, name of reference file for SM by SM calib.
/// \param badReconstruction: Bool_t, flag to find L1 shift dur to errors in reconstruction
/// \param fillHeavyHistos: Bool_t, flag to fill heavy histograms with time per channel
/// \param badMapType: Int_t,settype of bad channel map acces
/// \param badMapFileName: TString, file with bad channels map (absID)
///
/// \author Adam Matyja <adam.tomasz.matyja@ifj.edu.pl>, INP PAN Cracow
///

AliAnalysisTaskEMCALTimeCalib  * AddTaskEMCALTimeCalibration(TString  outputFile = "", // timeResults.root
							     TString  geometryName = "",//EMCAL_COMPLETE12SMV1_DCAL_8SM
							     Double_t minClusterEne = 1.0,
							     Double_t maxClusterEne = 500,
							     Int_t    minNcells = 2,
							     Int_t    maxNcells = 200,
							     Double_t minLambda0LG = 0.1,
							     Double_t maxLambda0LG = 4.0,
							     Double_t minLambda0 = 0.1,
							     Double_t maxLambda0 = 0.4,
							     Double_t maxRtrack = 0.025,
							     Double_t minCellEne = 0.4,
							     Double_t minTime = -20.,
							     Double_t maxTime = 20.,
							     Bool_t   pileupFromSPDFlag = kFALSE,
							     TString  referenceFileName = "",//Reference.root
							     TString  referenceSMFileName = "",//ReferenceSM.root
							     Bool_t   badReconstruction = kFALSE,
							     Bool_t   fillHeavyHistos = kFALSE,
							     Int_t    badMapType = 0,
							     TString  badMapFileName = "")
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddTaskEMCALTimeCalibration", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
  {
    ::Error("AddTaskEMCALTimeCalibration", "This task requires an input event handler");
    return NULL;
  }
    
  AliAnalysisTaskEMCALTimeCalib *taskmbemcal = new AliAnalysisTaskEMCALTimeCalib("TimeCalibTask");
  taskmbemcal->SelectCollisionCandidates(AliVEvent::kEMC1|AliVEvent::kEMC7|AliVEvent::kEMC8|AliVEvent::kEMCEJE|AliVEvent::kEMCEGA);
  taskmbemcal->SetGeometryName(geometryName);
  taskmbemcal->SetMinClusterEnergy (minClusterEne);
  taskmbemcal->SetMaxClusterEnergy (maxClusterEne);
  taskmbemcal->SetMinNcells        (minNcells);  
  taskmbemcal->SetMaxNcells        (maxNcells);   
  taskmbemcal->SetMinLambda0       (minLambda0);	   
  taskmbemcal->SetMaxLambda0       (maxLambda0);	   
  taskmbemcal->SetMinLambda0LG     (minLambda0LG);	   
  taskmbemcal->SetMaxLambda0LG     (maxLambda0LG);	   
  taskmbemcal->SetMaxRtrack        (maxRtrack);	   
  taskmbemcal->SetMinCellEnergy    (minCellEne);	   
  taskmbemcal->SetMinTime          (minTime);	   
  taskmbemcal->SetMaxTime          (maxTime);

  if(fillHeavyHistos) taskmbemcal->SwithOnFillHeavyHisto();
  else taskmbemcal->SwithOffFillHeavyHisto();

  // pass1
  taskmbemcal->SetRawTimeHisto(200,400.,800.);
  taskmbemcal->SetPassTimeHisto (200,400.,800.);
  // pass2
  if(referenceSMFileName.Length()!=0){
    taskmbemcal->SetReferenceRunByRunFileName(referenceSMFileName);
    taskmbemcal->LoadReferenceRunByRunHistos();
    taskmbemcal->SetPassTimeHisto(800,400.,800.);
    if(badReconstruction) {  //add for runs before LHC15n muon_calo_pass1 in run2
      taskmbemcal->SwitchOnBadReco();
      taskmbemcal->SetPassTimeHisto(500,-100.,150.);
    }
  }

  //pass3
  if(referenceFileName.Length()!=0){   
    taskmbemcal->SetReferenceFileName(referenceFileName);
    taskmbemcal->LoadReferenceHistos();
    taskmbemcal->SetPassTimeHisto(1000,-250.,250.);
  }
  if(pileupFromSPDFlag==kTRUE) taskmbemcal->SwitchOnPileupFromSPD();
  else taskmbemcal->SwitchOffPileupFromSPD();

  //bad channel map
  taskmbemcal->SetBadChannelMapSource(badMapType);
  if(badMapType==2) taskmbemcal->SetBadChannelFileName(badMapFileName);


  //taskmbemcal->PrintInfo();

  if(outputFile.Length()==0) outputFile = AliAnalysisManager::GetCommonFileName(); 

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("chistolist", TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   outputFile.Data());

  mgr->AddTask(taskmbemcal);
  mgr->ConnectInput  (taskmbemcal, 0, cinput1);
  mgr->ConnectOutput (taskmbemcal, 1, coutput);

  return taskmbemcal;
}
