/**
 * @file   MakeFMDEventPlaneTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:52:39 2012
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_trains
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ \Psi_R@f$
 * 
 *
 * @ingroup pwglf_forward_flow
 * @ingroup pwglf_forward_trains
 */
class MakeFMDEventPlaneTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train (free form)
   */
  MakeFMDEventPlaneTrain(const char* name) 
  : TrainSetup(name)
  {
    SetType(kAOD);
  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param par  Whether to use par files 
   */
  void CreateTasks(EMode /*mode*/, Bool_t par, AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("AnalysisResults.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2", par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2:"
			     "$ALICE_ROOT/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskFMDEventPlane.C(%d)", fMC));
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
  AliVEventHandler* CreateOutputHandler(EType) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeFMDEventPlaneTrain"; }
};
//
// EOF
//
