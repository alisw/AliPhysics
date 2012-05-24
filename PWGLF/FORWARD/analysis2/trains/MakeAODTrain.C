#include "TrainSetup.C"

/**
 * @file   MakeAODTrain.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:05:30 2011
 * 
 * @brief  Run first pass analysis - make AOD tree
 * 
 * @ingroup pwglf_forward_trains
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
 * @ingroup pwglf_forward_aod
 * @ingroup pwglf_forward_trains
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
      fUseCent(useCent), 
      fDebugLvl(0), 
      fUseTPCEventPlane(false)
  {}
  /** 
   * Set the debug level on Forward objects
   * 
   * @param lvl Debug level
   */
  void SetDebugLevel(Int_t lvl) { fDebugLvl = lvl; }
  /** 
   * If set to true, add TPC event plane task. 
   * 
   * @param use Wheter to include TPC event plane task 
   */
  void SetUseTPCEventPlance(Bool_t use) { fUseTPCEventPlane = use; }
  /** 
   * Run this analysis 
   * 
   * @param mode     Mode - see TrainSetup::EMode
   * @param oper     Operation - see TrainSetup::EOperation
   * @param nEvents  Number of events (negative means all)
   * @param mc       If true, assume simulated events 
   * @param usePar   If true, use PARs 
   */
  void Run(const char* mode, 
	   const char* oper, 
	   Int_t       nEvents=-1, 
	   Bool_t      mc=false,
	   Bool_t      usePar=false, 
	   Int_t       dbg=0)
  {
    Info("Run", "Running in mode=%s, oper=%s, events=%d, MC=%d, Par=%d", 
	 mode, oper, nEvents, mc, usePar);
    Exec("ESD", mode, oper, nEvents, mc, usePar, dbg);
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
  void Run(EMode  mode, 
	   EOper  oper, 
	   Int_t  nEvents=-1, 
	   Bool_t mc=false, 
	   Bool_t usePar = false,
	   Int_t  dbg = 0)
  {
    Info("Run", "Running in mode=%d, oper=%d, events=%d, MC=%d, Par=%d", 
	 mode, oper, nEvents, mc, usePar);
    Exec(kESD, mode, oper, nEvents, mc, usePar, dbg);
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
    LoadLibrary("PWGLFforward2", mode, par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;
    
    // --- Add TPC eventplane task
    if (fUseTPCEventPlane) gROOT->Macro("AddTaskEventplane.C");

    // --- Task to copy header information ---------------------------
    gROOT->Macro("AddTaskCopyHeader.C");

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwardMult.C(%d,%d,%d,%d,%d)", 
		      mc, fSys, fSNN, fField, fDebugLvl));
    AddExtraFile(gSystem->Which(gROOT->GetMacroPath(), "ForwardAODConfig.C"));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskCentralMult.C(%d,%d,%d,%d,%d)", 
		      mc, fSys, fSNN, fField, fDebugLvl));
    AddExtraFile(gSystem->Which(gROOT->GetMacroPath(), "CentralAODConfig.C"));

    // --- Add MC particle task --------------------------------------
    if (mc) gROOT->Macro("AddTaskMCParticleFilter.C");

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
  //__________________________________________________________________
  /** 
   * Create the setup script.  This is overloaded here, so that we can
   * create our dN/deta job script here.
   * 
   * @param type   Type of analysis
   * @param mode   Mode of the analysis 
   * @param mc     Whether this is MC or not 
   * @param usePar Whether to use par files or not 
   * @param dbg    Debug level
   */
  virtual void CreateSetupScript(EType  type, 
				 EMode  mode, 
				 Bool_t mc, 
				 Bool_t usePar, 
				 Int_t  dbg) const
  {
    TrainSetup::CreateSetupScript(type, mode, mc, usePar, dbg);
    
    Info("CreateSetupScript", "Creating dNdeta train setup script");
    TString base(Form("%s_dndeta", fName.Data()));
    TString escaped = EscapeName(base, fDatime);
    std::ofstream o(Form("%s.C", escaped.Data()));
    if (!o) { 
      Error("CreateSetupScript", "Failed to make dNdeta script %s.C", 
	    escaped.Data());
      return;
    }
    
    o << std::boolalpha 
      << "// Script to analyse AOD pass " << EscapedName() << " for dN/deta\n"
      << "// Automatically generated by MakeAODTrain\n"
      << "void " << escaped << "(bool terminate=false,"
      << "const char* trig=\"INEL\",Double_t vzMin=-10, Double_t vzMax=+10,"
      << "const char* scheme=\"FULL\")\n"
      << "{\n";
    WriteBuild(o, "MakedNdetaTrain");

    o << "  MakedNdetaTrain t(\"" << base << "\",trig,vzMin,vzMax,scheme,"
      << fUseCent << ");\n";
    TrainSetup::WriteSettings(o, "t");
    o << "  t.SetDataDir(\"" << GetOutputDirectory(mode) << "\");\n";
    WriteRuns(o, "t");
      
    const char* cmode = (mode == kLocal ? "\"LOCAL\"" : 
			 mode == kProof ? "\"PROOF\"" : "\"GRID\"");
    o << "  t.Run(" << cmode << ",(terminate ? \"TERMINATE\" : \"FULL\"),-1,"
      << usePar << ',' << dbg << ");\n"
      << "}\n"
      << "// EOF" << std::endl;
    
    o.close();
  }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeAODTrain"; }
  //__________________________________________________________________
  void WriteConstruction(std::ostream& o, const char* obj) const
  {
    o << "  UShort_t sys = " << fSys     << "; // 1:pp, 2:PbPb, 3:pPb\n"
      << "  UShort_t sNN = " << fSNN     << "; // sqrt(sNN) in GeV\n"
      << "  Short_t  fld = " << fField   << "; // L3 field in kG\n"
      << "  Bool_t   cen = " << fUseCent << "; // enable centrality\n" 
      << "  MakeAODTrain "<< obj << "(\"" << fName << "\",sys,sNN,fld,cen);\n"; 
  }
  //__________________________________________________________________
  void WriteSettings(std::ostream& o, const char* obj) const
  {
    TrainSetup::WriteSettings(o, obj);
    o << "  " << obj << ".SetDebugLvl(" << fDebugLvl << ");\n"
      << "  " << obj << ".SetUseTPCEventPlane(" << fUseTPCEventPlane << ");\n"
      << std::endl;
  }
  //__________________________________________________________________
  void WriteRun(std::ostream& o, 
		const char* obj, 
		const char* /*type*/, 
		const char* mode, 
		const char* oper, 
		Bool_t      mc, 
		Bool_t      usePar, 
		Int_t       dbg) const
  {
    o << "  " << obj << ".Run(" << mode << ',' << oper << ",-1," << mc << ','
      << usePar << ',' << dbg << ");" << std::endl;
  }
  UShort_t fSys;
  UShort_t fSNN;
  Short_t  fField;
  Bool_t   fUseCent;
  Int_t    fDebugLvl;
  Bool_t   fUseTPCEventPlane;
};
//
// EOF
//
