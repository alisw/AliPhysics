/**
 * @file   MakeMultClassTrain.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:05:30 2011
 * 
 * @brief  Run first pass analysis - make AOD tree
 * 
 * @ingroup pwglf_forward_trains_specific
 */
#include "TrainSetup.C"
#include <sstream>

//====================================================================
/**
 * Analysis train to make Forward and Central multiplicity
 * 
 *
 * @ingroup pwglf_forward_aod
 * @ingroup pwglf_forward_trains_specific
 */
class MakeMultClassTrain : public TrainSetup
{
public:
  /** 
   * Constructor. 
   * 
   * @param name     Name of train (free form)
   */
  MakeMultClassTrain(const  TString& name) 
    : TrainSetup(name)
  {
    fOptions.Set("type", "ESD");
  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("multClass.root");

    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;

    gROOT->ProcessLine("AliMultEventClassifierTask* task = "
		       "new AliMultEventClassifierTask(\"MultEvent\"); "
		       "task->Connect(0)");
  }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeMultClassTrain"; }
};
//
// EOF
//
