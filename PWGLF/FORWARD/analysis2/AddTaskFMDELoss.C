/**
 * @file   AddTaskFMDELoss.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:14:03 2011
 * 
 * @brief  Add energy loss task to train. 
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/**
 * @defgroup pwglf_forward_eloss Energy Loss Fits
 *
 * Fitting the energy loss @f$\Delta/\Delta_{mip}@f$ spectra 
 *
 * @ingroup pwglf_forward_topical
 */

/**
 * This is the macro to include the FMD energy fitter in a train.  
 * 
 * @param mc        Assume MC input 
 * @param onlyMB    Only collect statistics for MB (INEL) events
 * @param config    Configuration script 
 * @param corrs     Corrections
 * @param lowCut    (Optional) lower cut 
 * @param extraDead File with extra dead channels 
 *
 * @return Newly created task 
 *
 * @ingroup pwglf_forward_eloss
 */
AliAnalysisTask*
AddTaskFMDELoss(Bool_t        mc, 
		Bool_t        onlyMB=false,
		Double_t      lowCut=-1,
		const Char_t* config="elossFitConfig.C",
		const Char_t* corrs="",
		const Char_t* extraDead="")
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFMDELoss", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // --- Set alternative corrections path ----------------------------
  AliForwardCorrectionManager& cm = AliForwardCorrectionManager::Instance();
  if (corrs && corrs[0] != '\0') cm.SetPrefix(corrs);

  // --- Make the task and add it to the manager ---------------------
  AliFMDEnergyFitterTask* task = new AliFMDEnergyFitterTask("ForwardELoss");
  // --- Set parameters on the algorithms ----------------------------
  task->Configure(config);

  // For MC input we explicitly disable the noise correction 
  if (mc) task->GetESDFixer().SetRecoNoiseFactor(4);

  // If we get extra dead map from outside, add it here
  if (extraDead && extraDead[0] != '\0') task->GetESDFixer().AddDead(extraDead);
  // If we got a low cut value, set it 
  if (lowCut > 0)   task->GetEnergyFitter().SetLowCut(lowCut);
  
  // --- General -----------------------------------------------------
  // If set, only collect statistics for MB.  This is to prevent a
  // bias when looping over data where the MB trigger is downscaled.
  task->SetOnlyMB(onlyMB);

  // --- Make the output container and connect it --------------------
  task->Connect(0,0);

  Printf("Returning task %p", task);
  return task;
}
//
// EOF
//
