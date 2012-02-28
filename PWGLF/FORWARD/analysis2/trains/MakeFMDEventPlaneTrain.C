#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ flow@f$
 * 
 * To run, do 
 * @code 
 * gROOT->LoadMacro("TrainSetup.C");
 * // Make train 
 * MakeFMDEventPlaneTrain t("My Analysis");
 * // Set variaous parameters on the train 
 * t.SetDataDir("/home/of/data");
 * t.AddRun(118506)
 * // Run it 
 * t.Run("LOCAL", "FULL", -1, false, false);
 * @endcode 
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
   * @param dateTime Append date and time to name 
   * @param year     Year     - if not specified, current year
   * @param month    Month    - if not specified, current month
   * @param day      Day      - if not specified, current day
   * @param hour     Hour     - if not specified, current hour
   * @param min      Minutes  - if not specified, current minutes
   */
  MakeFMDEventPlaneTrain(const char* name, 
		Bool_t      mc = false,
		Bool_t      dateTime=false,
		UShort_t    year  = 0, 
		UShort_t    month = 0, 
		UShort_t    day   = 0, 
		UShort_t    hour  = 0, 
		UShort_t    min   = 0) 
    : TrainSetup(name, dateTime, year, month, day, hour, min),
      fMC(mc)
  {
  }
  /** 
   * Run this analysis 
   * 
   * @param mode     Mode - see TrainSetup::EMode
   * @param oper     Operation - see TrainSetup::EOperation
   * @param nEvents  Number of events (negative means all)
   * @param usePar   If true, use PARs 
   */
  void Run(const char* mode, const char* oper, 
	   Int_t nEvents=-1, Bool_t usePar=false)
  {
    Exec("AOD", mode, oper, nEvents, false, usePar);
  }
  /** 
   * Run this analysis 
   * 
   * @param mode     Mode - see TrainSetup::EMode
   * @param oper     Operation - see TrainSetup::EOperation
   * @param nEvents  Number of events (negative means all)
   * @param usePar   If true, use PARs 
   */
  void Run(EMode mode, EOper oper, Int_t nEvents=-1, 
	   Bool_t usePar=false)
  {
    Exec(kAOD, mode, oper, nEvents, false, usePar);
  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param mode Processing mode
   * @param par  Whether to use par files 
   */
  void CreateTasks(EMode mode, Bool_t par, AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("AnalysisResults.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2", mode, par, true);
    
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
  Bool_t fMC;
};
//
// EOF
//
