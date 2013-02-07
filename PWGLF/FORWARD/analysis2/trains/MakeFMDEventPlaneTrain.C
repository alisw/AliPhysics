/**
 * @file   MakeFMDEventPlaneTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:52:39 2012
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_trains_specific
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ \Psi_R@f$
 * 
 *
 * @ingroup pwglf_forward_flow
 * @ingroup pwglf_forward_trains_specific
 */
class MakeFMDEventPlaneTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  
   * 
   * @param name     Name of train (free form)
   */
  MakeFMDEventPlaneTrain(const char* name) 
  : TrainSetup(name)
  {
    fOptions.Set("type", "AOD");
  }
protected:
  /** 
   * Create the tasks 
   * 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("AnalysisResults.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2)");
    
    Bool_t mc   = mgr->GetMCtruthEventHandler() != 0;

    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2:"
			     "$ALICE_ROOT/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskFMDEventPlane.C(%d)", mc));
  }
  /** 
   * Do not the centrality selection
   */
  void CreateCentralitySelection(Bool_t, AliAnalysisManager*) {}
  /** 
   * Crete output handler - we don't want one here. 
   * 
   * @return 0
   */
  AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeFMDEventPlaneTrain"; }
};
//
// EOF
//
