/**
 * This is the macro to include the Forward multiplicity in a train.  
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
AliAnalysisTask*
AddTaskFMD(Int_t nCutBins=1, Float_t correctionCut=0.1)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFMD", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task and add it to the manager ---------------------
  AliForwardMultiplicity* task = new AliForwardMultiplicity("FMD");
  mgr->AddTask(task);

  // --- Set parameters on the algorithms ----------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // Set the eta axis to use - note, this overrides whatever is used
  // by the rest of the algorithms - but only for the energy fitter
  // algorithm. 
  task->GetEnergyFitter().SetEtaAxis(200, -4, 6);
  // Set the low cut used for energy
  task->GetEnergyFitter().SetLowCut(0.4);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  task->GetEnergyFitter().SetBinsToSubtract(4);
  // Set the maximum number of landaus to try to fit (max 5)
  task->GetEnergyFitter().SetNLandau(5);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data
  task->GetEnergyFitter().SetMinEntries(1000);
  // Set maximum energy loss to consider 
  task->GetEnergyFitter().SetMaxE(10); 
  // Set number of energy loss bins 
  task->GetEnergyFitter().SetNEbins(300);
  // Set whether to use increasing bin sizes 
  task->GetEnergyFitter().SetUseIncreasingBins(true);
  // Set the low cut used for sharing 
  task->GetSharingFilter().SetLowCut(0.3);
  // Set the number of extra bins (beyond the secondary map border) 
  task->GetHistCollector().SetNCutBins(nCutBins);
  // Set the correction cut, that is, when bins in the secondary map 
  // is smaller than this, they are considered empty 
  task->GetHistCollector().SetCorrectionCut(correctionCut);
  // Set the overall debug level (1: some output, 3: a lot of output)
  task->SetDebug(0);
  // Set the debug level of a single algorithm 
  task->GetEnergyFitter().SetDebug(3);
  
  // --- Set up the parameer manager ---------------------------------
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  AliMCEventHandler* mcHandler = 
    dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  Info("AddTaskFMD", "MC handler %p", mcHandler);
  if(mcHandler) {
    pars->SetRealData(kFALSE);
    pars->SetProcessPrimary(kTRUE);
    pars->SetProcessHits(kFALSE);
  }
  else {
    pars->SetRealData(kTRUE);
    pars->SetProcessPrimary(kFALSE);
    pars->SetProcessHits(kFALSE);
  }
  pars->Init();
  
  // --- Makek the output container and connect it -------------------
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":%s",pars->GetDndetaAnalysisName());
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Forward", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histOut);

  return task;
}
