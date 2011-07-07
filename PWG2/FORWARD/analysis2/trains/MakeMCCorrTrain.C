#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make Forward and Central MC corrections
 * 
 * To run, do 
 * @code 
 * gROOT->LoadMacro("TrainSetup.C");
 * // Make train 
 * MakeMCCorrTrain t("My Analysis");
 * // Set variaous parameters on the train 
 * t.SetDataDir("/home/of/data");
 * t.AddRun(118506)
 * // Run it 
 * t.Run("LOCAL", "FULL", -1, false, false);
 * @endcode 
 *
 * @ingroup pwg2_forward_mc
 * @ingroup pwg2_forward_trains
 */
class MakeMCCorrTrain : public TrainSetup
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
  MakeMCCorrTrain(const  char* name, 
		  Bool_t       dateTime = false,
		  UShort_t     year     = 0, 
		  UShort_t     month    = 0, 
		  UShort_t     day      = 0, 
		  UShort_t     hour     = 0, 
		  UShort_t     min      = 0) 
    : TrainSetup(name, dateTime, 
		 year, month, day, hour, min)
  {}
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
    Info("Run", "Running in mode=%s, oper=%s, events=%d, Par=%d", 
	 mode, oper, nEvents, usePar);
    Exec("ESD", mode, oper, nEvents, true, usePar);
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
	   Bool_t usePar = false)
  {
    Info("Run", "Running in mode=%d, oper=%d, events=%d, Par=%d", 
	 mode, oper, nEvents, usePar);
    Exec(kESD, mode, oper, nEvents, true, usePar);
  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param mode Processing mode
   * @param par  Whether to use par files 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(EMode mode, Bool_t par, AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_mccorr.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWG2forward2", mode, par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2",
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
  Double_t fVzMin;     // Least v_z
  Double_t fVzMax;     // Largest v_z
};

//
// EOF
//
