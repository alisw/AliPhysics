/**
 * @file   AddTaskForwardQA.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:14:03 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 */
/**
 * @defgroup pwg2_forward_eloss Energy Loss Fits
 * @ingroup pwg2_forward_topical
 */

/**
 * This is the macro to include the FMD energy fitter in a train.  
 * 
 * @ingroup pwg2_forward_eloss
 */
AliAnalysisTask*
AddTaskForwardQA(Bool_t mc, Bool_t useCent)
{
  gSystem->Load("libPWG2forward2");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskForwardQA", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // --- Make the task and add it to the manager ---------------------
  AliForwardQATask* task = new AliForwardQATask("forwardQA");
  // task->SetBLow(blow);
  // task->SetBLow(bhigh);
  mgr->AddTask(task);

  Double_t nXi = mc ? 2 : 2;   //HHD test
  Bool_t   includeSigma = false; //true;

  AliFMDMultCuts cDensity;
  //c2.SetNXi(mc ? 1 : 1);
  //  c2.SetIncludeSigma(false);
  //c2.SetMPVFraction(0.5);
  //Double_t factor = 1.2;
  cDensity.SetMultCuts(0.3, 0.3, 0.3, 0.3, 0.3);
  
  // --- Set parameters on the algorithms ----------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);

  // --- Set parameters on energy loss fitter ------------------------
  // Set the eta axis to use - note, this overrides whatever is used
  // by the rest of the algorithms - but only for the energy fitter
  // algorithm. 
  task->GetEnergyFitter().SetEtaAxis(200, -4, 6);
  // Set maximum energy loss to consider 
  task->GetEnergyFitter().SetMaxE(15); 
  // Set number of energy loss bins 
  task->GetEnergyFitter().SetNEbins(500);
  // Set whether to use increasing bin sizes 
  task->GetEnergyFitter().SetUseIncreasingBins(true);
  // Set whether to do fit the energy distributions 
  task->GetEnergyFitter().SetDoFits(kTRUE);
  // Set whether to make the correction object 
  task->GetEnergyFitter().SetDoMakeObject(kFALSE);
  // Set the low cut used for energy
  task->GetEnergyFitter().SetLowCut(0.4);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  task->GetEnergyFitter().SetFitRangeBinWidth(4);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data
  task->GetEnergyFitter().SetMinEntries(1000);

  // --- Sharing filter ----------------------------------------------
  // Enable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetUseAngleCorrectedSignals(true);
  // Disable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetZeroSharedHitsBelowThreshold(false);
  // Enable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().GetHCuts().SetMultCuts(-1);
  // Set the number of xi's (width of landau peak) to stop at 
  task->GetSharingFilter().GetHCuts().SetNXi(nXi);
  // Set whether or not to include sigma in cut
  task->GetSharingFilter().GetHCuts().SetIncludeSigma(includeSigma);
  // Enable use of angle corrected signals in the algorithm 
  task->GetSharingFilter().SetLCuts(cDensity);

  // --- Density calculator ------------------------------------------
  // Set the maximum number of particle to try to reconstruct 
  task->GetDensityCalculator().SetMaxParticles(10);
  // Wet whether to use poisson statistics to estimate N_ch
  task->GetDensityCalculator().SetUsePoisson(true);
  // Set whether or not to include sigma in cut
  task->GetDensityCalculator().SetCuts(cDensity);
  // Set whether or not to use the phi acceptance
  //   AliFMDDensityCalculator::kPhiNoCorrect
  //   AliFMDDensityCalculator::kPhiCorrectNch
  //   AliFMDDensityCalculator::kPhiCorrectELoss
  task->GetDensityCalculator()
    .SetUsePhiAcceptance(AliFMDDensityCalculator::kPhiCorrectNch);
  // --- Set limits on fits the energy -------------------------------
  // Maximum relative error on parameters 
  AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
  // Least weight to use 
  AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
  // Maximum value of reduced chi^2 
  AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 20;
    
  // --- Make the output container and connect it --------------------
  // TString outputfile = ;
  // outputfile += ":PWG2forwardDnDeta"; 
  // Form(":%s",pars->GetDndetaAnalysisName());
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Forward", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("ForwardResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histOut);
  mgr->ConnectOutput(task, 2, output);


  return task;
}
