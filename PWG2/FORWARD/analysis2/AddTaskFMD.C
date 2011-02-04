/** 
 * @defgroup pwg2_forward_scripts Scripts used in the analysis
 *
 * @ingroup pwg2_forward
 */
/**
 * @file 
 * @ingroup pwg2_forward_scripts
 * 
 */
/**
 * This is the macro to include the Forward multiplicity in a train.  
 * 
 * @ingroup pwg2_forward_scripts
 */
AliAnalysisTask*
AddTaskFMD(Bool_t mc, UShort_t sys=0, UShort_t sNN=0, Short_t field)
{
  gSystem->Load("libPWG2forward2");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFMD", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task and add it to the manager ---------------------
  AliForwardMultiplicityBase* task = 0;
  if (mc) task = new AliForwardMCMultiplicityTask("FMD");
  else    task = new AliForwardMultiplicityTask("FMD");
  mgr->AddTask(task);
  
  // --- Do a local initialisation with assumed values ---------------
  if (sys > 0 && sNN > 0) 
    AliForwardCorrectionManager::Instance().Init(sys,sNN,field);
  
  // Whether to enable low flux specific code 
  task->SetEnableLowFlux(kFALSE);
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // Set the eta axis to use - note, this overrides whatever is used
  // by the rest of the algorithms - but only for the energy fitter
  // algorithm. 
  task->GetEnergyFitter().SetEtaAxis(200, -4, 6);
  // Set maximum energy loss to consider 
  task->GetEnergyFitter().SetMaxE(10); 
  // Set number of energy loss bins 
  task->GetEnergyFitter().SetNEbins(300);
  // Set whether to use increasing bin sizes 
  task->GetEnergyFitter().SetUseIncreasingBins(true);
  // Set whether to do fit the energy distributions 
  task->GetEnergyFitter().SetDoFits(kFALSE);
  // Set whether to make the correction object 
  task->GetEnergyFitter().SetDoMakeObject(kFALSE);
  // Set the low cut used for energy
  task->GetEnergyFitter().SetLowCut(0.4);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  task->GetEnergyFitter().SetFitRangeBinWidth(4);
  // Set the maximum number of landaus to try to fit (max 5)
  task->GetEnergyFitter().SetNParticles(5);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data
  task->GetEnergyFitter().SetMinEntries(1000);
  // Set the low cut used for sharing - overrides settings in eloss fits
    task->GetSharingFilter().SetLowCut(0.3);
  // Set the number of xi's (width of landau peak) to stop at 
  task->GetSharingFilter().SetNXi(1);
  // Set the maximum number of particle to try to reconstruct 
  task->GetDensityCalculator().SetMaxParticles(2);
  // Set the lower multiplicity cut.  Overrides setting in energy loss fits.
  task->GetDensityCalculator().SetMultCut(0.3); //was 0.3
  // Whether to use the secondary map correction
  task->GetCorrections().SetUseSecondaryMap(true);
  // Whether to use the vertex bias correction
  task->GetCorrections().SetUseVertexBias(false);
  // Whether to use the vertex bias correction
  task->GetCorrections().SetUseAcceptance(true);
  // Whether to use the merging efficiency correction 
  task->GetCorrections().SetUseMergingEfficiency(false);
  // Set the number of extra bins (beyond the secondary map border) 
  task->GetHistCollector().SetNCutBins(2);
  // Set the correction cut, that is, when bins in the secondary map 
  // is smaller than this, they are considered empty 
  task->GetHistCollector().SetCorrectionCut(0.5);
  // Set the overall debug level (1: some output, 3: a lot of output)
  task->SetDebug(0);
  // Set the debug level of a single algorithm 
  // task->GetEventInspector().SetDebug(4);
  // --- Set limits on fits the energy -------------------------------
  // Maximum relative error on parameters 
  AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
  // Least weight to use 
  AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
  // Maximum value of reduced chi^2 
  AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 5;

  
  // --- Make the output container and connect it --------------------
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  // outputfile += ":PWG2forwardDnDeta"; 
  // Form(":%s",pars->GetDndetaAnalysisName());
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Forward", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histOut);

  return task;
}
