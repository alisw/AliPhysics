/**
 * @file   MakeRespMatrTrain.C
 * @author Valentina Zaccolo <Valentina.Zaccolo@cern.ch>
 * @date   Fri Jan  11 14:47:26 2013
 * 
 * @brief  Make response matrices 
 * 
 * @ingroup pwglf_forward_trains_specific
 * @ingroup pwglf_forward_multdist
 * @ingroup pwglf_forward_scripts_tasks_vz
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make response matrices.
 * Valentina's modified task 
 *
 * @ingroup pwglf_forward_multdist
 * @ingroup pwglf_forward_trains_specific
 * @ingroup pwglf_forward_scripts_tasks_vz
 */
class MakeRespMatrTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  
   * 
   * @param name     Name of train (free form)
   */
  MakeRespMatrTrain(const char* name)
    : TrainSetup(name)
  {
    fOptions.Add("trig",    "TYPE",       "Trigger type", "V0AND");
    fOptions.Add("filter",  "TYPE",       "Filter type", "OUTLIER|PILEUP-BIN");
    fOptions.Add("ipz-min", "CENTIMETER", "Min Ip Z",     -4);
    fOptions.Add("ipz-max", "CENTIMETER", "Max Ip Z",     +4);
    fOptions.Add("max-mult","NUMBER",     "Max of histograms", 500);
    fOptions.Add("assume-nsd", "If set, assume MC events are NSD", false);
  }
protected:
  /** 
   * Make the tasks 
   * 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_response.root");

    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Add the task ----------------------------------------------
    AliAnalysisTaskSE* tsk = CoupleSECar("AddTaskCreateRespMatr.C",
					 AliVEvent::kAny);

    // --- Set options -----------------------------------------------
    FromOption(tsk, "TriggerMask",         "trig",     "V0AND");
    FromOption(tsk, "FilterMask",          "filter",   "OUTLIER|PILEUP-BIN");
    FromOption(tsk, "IpZMin",              "ipz-min",  -4.);
    FromOption(tsk, "IpZMax",              "ipz-max",  +4.);    
    FromOption(tsk, "MaxMult",             "max-mult", 500);
    FromOption(tsk, "MCIsNSD",             "assume-nsd", false);
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
  /** 
   * Never ever make an MC input handler! (we're working on AODs,
   * right?)
   */
  AliVEventHandler* CreateMCHandler(UShort_t type, bool mc) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeRespMatrTrain"; }
  //__________________________________________________________________
  
};
//
// EOF
//
