/**
 * @file   MakeAODTrain.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:05:30 2011
 * 
 * @brief  Run first pass analysis - make AOD tree
 * 
 * @ingroup pwg2_forward_trains
 */
//====================================================================
/**
 * Analysis train to make Forward and Central multiplicity
 * 
 * To run, do 
 * @code 
 * gROOT->LoadMacro("TrainSetup.C");
 * // Make train 
 * MakeAODTrain t("My Analysis");
 * // Set variaous parameters on the train 
 * t.SetDataDir("/home/of/data");
 * t.AddRun(118506)
 * // Run it 
 * t.Run("LOCAL", "FULL", -1, false, false);
 * @endcode 
 *
 * @ingroup pwg2_forward_aod
 * @ingroup pwg2_forward_trains
 */
class MakeAODTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train (free form)
   * @param sys      Collision system (1: pp, 2: PbPb)
   * @param sNN      Center of mass energy [GeV]
   * @param field    L3 magnetic field - one of {-5,0,+5} kG
   * @param useCent  Whether to use centrality 
   * @param dateTime Append date and time to name 
   * @param year     Year     - if not specified, current year
   * @param month    Month    - if not specified, current month
   * @param day      Day      - if not specified, current day
   * @param hour     Hour     - if not specified, current hour
   * @param min      Minutes  - if not specified, current minutes
   */
  MakeAODTrain(const  char* name, 
	       UShort_t     sys      = 0, 
	       UShort_t     sNN      = 0, 
	       Short_t      field    = 0, 
	       Bool_t       useCent  = false, 
	       Bool_t       dateTime = false, 
	       UShort_t     year     = 0, 
	       UShort_t     month    = 0, 
	       UShort_t     day      = 0, 
	       UShort_t     hour     = 0, 
	       UShort_t     min      = 0) 
    : TrainSetup(name, dateTime, 
		 year, month, day, hour, min),
      fSys(sys), 
      fSNN(sNN), 
      fField(field),
      fUseCent(useCent)
  {}
  /** 
   * Run this analysis 
   * 
   * @param mode     Mode - see TrainSetup::EMode
   * @param oper     Operation - see TrainSetup::EOperation
   * @param nEvents  Number of events (negative means all)
   * @param mc       If true, assume simulated events 
   * @param usePar   If true, use PARs 
   */
  void Run(const char* mode, const char* oper, 
	   Int_t nEvents=-1, Bool_t mc=false,
	   Bool_t usePar=false)
  {
    Info("Run", "Running in mode=%s, oper=%s, events=%d, MC=%d, Par=%d", 
	 mode, oper, nEvents, mc, usePar);
    Exec("ESD", mode, oper, nEvents, mc, usePar);
  }
  /** 
   * Run this analysis 
   * 
   * @param mode     Mode - see TrainSetup::EMode
   * @param oper     Operation - see TrainSetup::EOperation
   * @param nEvents  Number of events (negative means all)
   * @param mc       If true, assume simulated events 
   * @param usePar   If true, use PARs 
   */
  void Run(EMode mode, EOper oper, Int_t nEvents=-1, Bool_t mc=false, 
	   Bool_t usePar = false)
  {
    Info("Run", "Running in mode=%d, oper=%d, events=%d, MC=%d, Par=%d", 
	 mode, oper, nEvents, mc, usePar);
    Exec(kESD, mode, oper, nEvents, mc, usePar);
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
    AliAnalysisManager::SetCommonFileName("forward.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWG2forward2", mode, par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;
    
    // --- Task to copy header information ---------------------------
    gROOT->Macro("AddTaskCopyHeader.C");

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwardMult.C(%d,%d,%d,%d)", 
		      mc, fSys, fSNN, fField));
    AddExtraFile(gSystem->Which(gROOT->GetMacroPath(), "ForwardAODConfig.C"));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskCentralMult.C(%d,%d,%d,%d)", 
		      mc, fSys, fSNN, fField));
    AddExtraFile(gSystem->Which(gROOT->GetMacroPath(), "CentralAODConfig.C"));
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
      Fatal("CreatePhysicsSelection", "Couldn't get PhysicsSelection (%p)",ps);

    // --- Ignore trigger class when selecting events.  This means ---
    // --- that we get offline+(A,C,E) events too --------------------
    // ps->SetSkipTriggerClassSelection(true);
  }
  //__________________________________________________________________
  /** 
   * Create the centrality selection only if requested
   * 
   * @param mc  Monte-Carlo truth flag 
   * @param mgr Manager
   */
  void CreateCentralitySelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    if (!fUseCent) return;
    TrainSetup::CreateCentralitySelection(mc, mgr);
  }
  UShort_t fSys;
  UShort_t fSNN;
  Short_t  fField;
  Bool_t   fUseCent;
};
//
// EOF
//
