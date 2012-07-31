/**
 * @file   MakeMCCorrTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:54:47 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make Forward and Central MC corrections
 * 
 *
 * @ingroup pwglf_forward_mc
 * @ingroup pwglf_forward_trains
 */
class MakeMCCorrTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train (free form)
   */
  MakeMCCorrTrain(const  char* name) 
  : TrainSetup(name)
  {
    SetType(kESD);
  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param par  Whether to use par files 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(EMode /*mode*/, Bool_t par, AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_mccorr.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2", par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    if (!mgr->GetMCtruthEventHandler()) return;
    
    // --- Add the task ----------------------------------------------
    gROOT->Macro("AddTaskForwardMCCorr.C"); 

    // --- Add the task ----------------------------------------------
    gROOT->Macro("AddTaskCentralMCCorr.C");
  }
  //__________________________________________________________________
  /** 
   * Create physics selection , and add to manager
   * 
   * @param mc Whether this is for MC 
   * @param mgr Manager 
   */
  void CreatePhysicsSelection(Bool_t mc,
			      AliAnalysisManager* mgr)
  {
    TrainSetup::CreatePhysicsSelection(mc, mgr);

    // --- Get input event handler -----------------------------------
    AliInputEventHandler* ih =
      dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
    if (!ih) 
      Fatal("CreatePhysicsSelection", "Couldn't get input handler (%p)", ih);
    
    // --- Get Physics selection -------------------------------------
    AliPhysicsSelection* ps = 
      dynamic_cast<AliPhysicsSelection*>(ih->GetEventSelection());
    if (!ps) 
      Fatal("CreatePhysicsSelection", "Couldn't get PhysicsSelection (%p)", ps);

    // --- Ignore trigger class when selecting events.  This means ---
    // --- that we get offline+(A,C,E) events too --------------------
    // ps->SetSkipTriggerClassSelection(true);
  }
  /** 
   * Do not the centrality selection
   */
  void CreateCentralitySelection(Bool_t, AliAnalysisManager*) {}
  //__________________________________________________________________
  const char* ClassName() const { return "MakeMCCorrTrain"; }
};

//
// EOF
//
