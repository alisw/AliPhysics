#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ flow@f$
 * 
 * To run, do 
 * @code 
 * gROOT->LoadMacro("TrainSetup.C");
 * // Make train 
 * MakeFlowTrain t("My Analysis");
 * // Set variaous parameters on the train 
 * t.SetDataDir("/home/of/data");
 * t.AddRun(118506)
 * // Run it 
 * t.Run("LOCAL", "FULL", -1, false, false);
 * @endcode 
 *
 * @ingroup pwg2_forward_flow
 * @ingroup pwg2_forward_trains
 */
class MakeFlowTrain : public TrainSetup
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
  MakeFlowTrain(const char* name, 
		char*       type="", 
		Bool_t      mc = false,
		char*       addFlow = "",
		Int_t       addFType = 0,
		Int_t       addFOrder = 0,
		Bool_t      dateTime=false,
		UShort_t    year  = 0, 
		UShort_t    month = 0, 
		UShort_t    day   = 0, 
		UShort_t    hour  = 0, 
		UShort_t    min   = 0) 
    : TrainSetup(name, dateTime, year, month, day, hour, min),
      fType(type), 
      fMC(mc),
      fAddFlow(addFlow),
      fAddFType(addFType),
      fAddFOrder(addFOrder)
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
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwardFlow.C(\"%s\",%d,\"%s\",%d,%d)",
		      fType, fMC, fAddFlow, fAddFType, fAddFOrder));
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
  char* fType;
  Bool_t fMC;
  char* fAddFlow;
  Int_t fAddFType;
  Int_t fAddFOrder;
};
//
// EOF
//
