/**
 * @file   MakeMultTrain.C
 * @author Valentina Zaccolo <Valentina.Zaccolo@cern.ch>
 * @date   Wed Nov  21 12:47:26 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains_specific
 * @ingroup pwglf_forward_multdist
 * 
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ Multiplicity Distributions@f$.
 * Valentina's modified task 
 *
 * @ingroup pwglf_forward_multdist
 * @ingroup pwglf_forward_trains_specific
 * @ingroup pwglf_forward_scripts_tasks_vz
 */
class MakeMultTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  
   * 
   * @param name     Name of train (free form)
   */
 MakeMultTrain(const char* name)
  : TrainSetup(name)
  {
    fOptions.Add("trig",    "TYPE",       "Trigger type",     "V0AND");
    fOptions.Add("ipz-min", "CENTIMETER", "Min Ip Z",         -4);
    fOptions.Add("ipz-max", "CENTIMETER", "Max Ip Z",         +4);
    // fOptions.Add("cent-min","%",          "Min Centrality",   0); 
    // fOptions.Add("cent-max","%",          "Max Centrality",   0);
    fOptions.Add("max-mult","N",          "Max Multiplicity", 500);
  }
protected:
  /** 
   * Create the tasks 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_multiplicity.root");

    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Add the task ----------------------------------------------
    AliAnalysisTaskSE* tsk =
      CoupleSECar("AddTaskMultDistributions.C", "", AliVEvent::kAny);
    FromOption(tsk, "TriggerMask",         "trig",     "V0AND");
    FromOption(tsk, "IpZMin",              "ipz-min",  -4.);
    FromOption(tsk, "IpZMax",              "ipz-max",  +4.);    
    FromOption(tsk, "MaxMult",             "max-mult", 500);
    // FromOption(tsk, "CentMin",             "cent-min", 0);
    // FromOption(tsk, "CentMax",             "cent-max", 100);
    
  }
  //__________________________________________________________________
  /** 
   * Do not the centrality selection
   */
  void CreateCentralitySelection(Bool_t) {}
  //__________________________________________________________________
  /** 
   * Crete output handler - we don't want one here. 
   * 
   * @return 0
   */
  AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeMultTrain"; }
};
//
// EOF
//
