/** 
 * @defgroup pwglf_forward_scripts Scripts used in the analysis
 *
 * @ingroup pwglf_forward
 */
/** 
 * @defgroup pwglf_forward_scripts_tasks Scripts to add tasks to manager 
 * @ingroup pwglf_forward_scripts
 */
/**
 * @file   AddTaskForwardMult.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:13:54 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/**
 * This is the script to include the Forward multiplicity in a train.  
 * 
 * @param mc    Define as true for MC input. 
 * @param sys   Collision system (0: deduce, 1: pp, 2: pbpb, 3:pA)
 * @param sNN   Collision energy 
 * @param field L3 field setting. 
 *
 * @return newly allocated analysis task 
 *
 * @ingroup pwglf_forward_aod
 */
AliAnalysisTask*
AddTaskForwardMult(Bool_t mc, UShort_t sys=0, UShort_t sNN=0, Short_t field=0)
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskForwardMult", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task and add it to the manager ---------------------
  AliForwardMultiplicityBase* task = 0;
  
  if (mc)
    task = new AliForwardMCMultiplicityTask("FMD");
  else    
    task = new AliForwardMultiplicityTask("FMD");
  mgr->AddTask(task);
  
  // --- Do a local initialisation with assumed values ---------------
  if (sys > 0 && sNN > 0) 
    AliForwardCorrectionManager::Instance().Init(sys,sNN,field,mc);

  // --- Configure the task ------------------------------------------
  TString macroPath(gROOT->GetMacroPath());
  if (!macroPath.Contains("$(ALICE_ROOT)/PWGLF/FORWARD/analysis2")) { 
    macroPath.Append(":$(ALICE_ROOT)/PWGLF/FORWARD/analysis2");
    gROOT->SetMacroPath(macroPath);
  }
  const char* config = gSystem->Which(gROOT->GetMacroPath(),
				      "ForwardAODConfig.C");
  if (!config) 
    Warning("AddTaskForwardMult", "ForwardAODConfig.C not found in %s",
	    gROOT->GetMacroPath());
  else {
    Info("AddTaskForwardMult", 
	 "Loading configuration of '%s' from %s",
	 task->ClassName(), config);
    gROOT->Macro(Form("%s((AliForwardMultiplicityBase*)%p)", config, task));
    delete config;
  }
  
  // --- Make the output container and connect it --------------------
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  // outputfile += ":PWGLFforwardDnDeta"; 
  // Form(":%s",pars->GetDndetaAnalysisName());
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Forward", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histOut);

  return task;
}

/**
 * This is the script to include the Forward multiplicity in a train.  
 * 
 * @param type   Data type (if it contains MC, assume MC input): 
 *               - ppb, p-pb, pa, p-a:  proton-lead 
 *               - pp, p-p:             proton-proton
 *               - pbpb, pb-pb, a-a:    lead-lead
 *               
 * @param energy Collision energy in GeV
 * @param bfield L3 field setting in kG (-5, 0, 5)
 *
 * @return newly allocated analysis task 
 *
 * @ingroup pwglf_forward_aod
 */
AliAnalysisTask*
AddTaskForwardMult(const Char_t* type, 
		   Float_t       energy=0, 
		   Float_t       bfield=0)
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Deduce parameters -------------------------------------------
  TString  t(type);
  Bool_t   mc    = t.Contains("MC", TString::kIgnoreCase);
  UShort_t sys   = AliForwardUtil::ParseCollisionSystem(type);
  UShort_t sNN   = AliForwardUtil::ParseCenterOfMassEnergy(sys, energy);
  Short_t  field = AliForwardUtil::ParseMagneticField(field);

  return AddTaskForwardMult(mc, sys, sNN, field);
}
//
// EOF
//
