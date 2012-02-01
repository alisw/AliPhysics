//====================================================================
/**
 * Analysis train to do energy loss fits
 * 
 * @ingroup pwg2_forward_trains
 */
class MakeFMDELossTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train 
   * @param useCent  Whether to use centrality or not 
   * @param dateTime Append date and time to name
   * @param year     Year
   * @param month    Month 
   * @param day      Day
   * @param hour     Hour 
   * @param min      Minutes
   */
  MakeFMDELossTrain(const char* name  = "FMD Energy Loss",
		    Bool_t   useCent  = false,
		    Bool_t   dateTime = false, 
		    UShort_t year     = 0, 
		    UShort_t month    = 0, 
		    UShort_t day      = 0, 
		    UShort_t hour     = 0, 
		    UShort_t min      = 0) 
    : TrainSetup(name, dateTime, year, month, day, hour, min), 
      fUseCent(useCent)
  {}
  //__________________________________________________________________
  /** 
   * Run this analysis 
   * 
   * @param mode     Mode
   * @param oper     Operation
   * @param nEvents  Number of events (negative means all)
   * @param mc       If true, assume simulated events 
   * @param par      IF true, use par files 
   */
  void Run(const char* mode, const char* oper, 
	   Int_t nEvents=-1, Bool_t mc=false, Bool_t par=false)
  {
    Exec("ESD", mode, oper, nEvents, mc, par);
  }
  //__________________________________________________________________
  /** 
   * Run this analysis 
   * 
   * @param mode     Mode
   * @param oper     Operation
   * @param nEvents  Number of events (negative means all)
   * @param mc       If true, assume simulated events 
   * @param par      IF true, use par files 
   */
  void Run(EMode mode, EOper oper, Int_t nEvents=-1, Bool_t mc=false,
	   Bool_t par=false)
  {
    Exec(kESD, mode, oper, nEvents, mc, par);
  }
protected:
  //__________________________________________________________________
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
    AliAnalysisManager::SetCommonFileName("forward_eloss.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWG2forward2", mode, par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskFMDELoss.C(%d,%d)", mc, fUseCent));
  }
  /** 
   * Create entrality selection if enabled 
   * 
   * @param mc   Whether this is MC or not
   * @param mgr  Analysis manager 
   */
  virtual void CreateCentralitySelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    if (!fUseCent) return;

    gROOT->Macro("AddTaskCentrality.C");
    AliCentralitySelectionTask* ctask = 
      dynamic_cast<AliCentralitySelectionTask*>(mgr->GetTask("CentralitySelection"));
    if (!ctask) return;
    ctask->SetPass(fESDPass);
    if (mc) ctask->SetMCInput();
  }
  /** 
   * Crete output handler - we don't want one here. 
   * 
   * @return 0
   */
  AliVEventHandler* CreateOutputHandler(EType) { return 0; }
  Bool_t fUseCent; // Whether to use centrality or not 
};

//
// EOF
//
