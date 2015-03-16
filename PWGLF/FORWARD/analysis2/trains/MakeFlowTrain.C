/**
 * @file   MakeFlowTrain.C
 * @author Alexander Hansen                                      
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains
 * 
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make flow
 * 
 *
 * @ingroup pwglf_forward_flow
 */
class MakeFlowTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  
   * 
   * @param name     Name of train (free form)
   */
  MakeFlowTrain(const char* name)
  : TrainSetup(name)
  {
    // General options
    fOptions.Add("sys",   "SYSTEM",  "1:pp, 2:PbPb, 3:pPb", "");
    fOptions.Add("mc", "Do MC analysis");
    fOptions.Add("max-mom", "2|3|4|5", "Max flow moment to analyse", "5");
    fOptions.Add("use-cent", "Whether to use the impact parameter "
		 "for centrality");
    fOptions.Add("mc-vtx", "Whether to get the vertex from the MC header");
    fOptions.Add("fwd-dets", "[fmd,vzero]", "Forward detectors", "fmd+vzero");
    fOptions.Add("afterburner", "Whether to afterburn");
    fOptions.Set("type", "AOD");
    fOptions.Show(std::cout);
    // QC specific options
    fOptions.Add("QC", "Add QC tasks (requires AODs with FMD and SPD objects)");
    fOptions.Add("qc-type", "[std,eta-gap,3cor,all]", 
		 "Which types of QC's to do", "all");
    fOptions.Add("eg-value", "EGVALUE", 
		 "Set value in |eta| of the eta gap", "2.0");
    fOptions.Add("tracks", "[tpc,hybrid,only]", 
		 "Whether or only to use tracks for reference flow", 
		 "tpc+hybrid");
    fOptions.Add("sat-vtx", "Whether to use satellite interactions");
    fOptions.Add("outlier-fmd", "NSIGMA", "Outlier cut for FMD", "4.0");
    fOptions.Add("outlier-spd", "NSIGMA", "Outlier cut for SPD", "4.0");
    // EP specific options
    fOptions.Add("EP", "Add EP tasks (requires VZERO AOD objects)");
  }
protected:
  //__________________________________________________________________
  /** 
   * Create the tasks 
   * 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_flow.root");

    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
                             gROOT->GetMacroPath()));

    // --- Add the tasks ---------------------------------------------
    Bool_t doQC = fOptions.AsBool("QC");
    Bool_t doEP = fOptions.AsBool("EP");
    if (doQC) AddQCTasks();
    if (doEP) AddEPTasks();
    if (!doQC && !doEP) Fatal("CreateTasks", "No flow tasks were added");

  }
  //__________________________________________________________________
  /**
   * Add the QC task(s)
   */
  void AddQCTasks()
  {
    // --- Get the parameters ----------------------------------------
    UShort_t sys      = fOptions.AsInt("sys", 0);
    Int_t    moment   = fOptions.AsInt("max-mom");
    TString  fwdDets  = fOptions.Get("fwd-dets");
    TString  types    = fOptions.Get("qc-type");
    Double_t egValue  = fOptions.AsDouble("eg-value");
    TString  tracks   = fOptions.Get("tracks");
    Bool_t   useCent  = fOptions.AsBool("use-cent");
    Bool_t   useMCVtx = fOptions.AsBool("mc-vtx");
    Bool_t   addFlow  = fOptions.AsBool("afterburner");
    Bool_t   satVtx   = fOptions.AsBool("sat-vtx");
    Double_t fmdCut   = fOptions.AsDouble("outlier-fmd");
    Double_t spdCut   = fOptions.AsDouble("outlier-spd");
    Bool_t   mc       = fOptions.AsBool("mc"); 

    types.ToLower();
    fwdDets.ToUpper();
    Bool_t doFMD      = fwdDets.Contains("FMD");
    Bool_t doVZERO    = fwdDets.Contains("VZERO");
    tracks.ToLower();
    Bool_t onlyTr     = tracks.Contains("only");
    Bool_t tpcTr      = tracks.Contains("tpc");
    Bool_t hybridTr   = tracks.Contains("hybrid");
    Bool_t ispA = (sys == 3 ? kTRUE : kFALSE);

    // Notice the place holders at arg=2,3,4, and 9, These are
    // 
    //   2: Detector to use (FMD/VZERO)
    //   3: Whether to use eta gap (true/false)
    //   4: Do 3-subevent correlations (true/false)
    //   9: Use tracks for referernce flow (true/false)
    TString args;
    args = TString::Format("%d,\"%%s\",%%d,%%d,%d,%f,%f,%f,%%d,%d,%d,%d,%d,%d",
			   moment,
			   mc, 
			   fmdCut, 
			   spdCut,
			   egValue,
			   useCent,
			   ispA,
			   useMCVtx,
			   satVtx, 
			   addFlow);
    
    // --- Add the task ----------------------------------------------
    const char* mac = "AddTaskForwardFlowQC.C";
    if (doFMD) {
      if (types.Contains("std") || types.Contains("all")) {
        if (!onlyTr)
	  CoupleCar(mac, Form(args.Data(), "FMD", false, false, 0));
      	if (tpcTr)
	  CoupleCar(mac, Form(args.Data(), "FMD", false, false, 1));
      	if (hybridTr)
	  CoupleCar(mac, Form(args.Data(), "FMD", false, false, 2));
      }
      if (types.Contains("eta-gap") || types.Contains("all")) {
        if (!onlyTr)
	  CoupleCar(mac, Form(args.Data(), "FMD", true, false, 0));
      	if (tpcTr)
      	  CoupleCar(mac, Form(args.Data(), "FMD", true, false, 1));
      	if (hybridTr)
      	  CoupleCar(mac, Form(args.Data(), "FMD", true, false, 2));
      }
      if (types.Contains("3cor") || types.Contains("all")) {
        if (!onlyTr)
	  CoupleCar(mac, Form(args.Data(), "FMD", false, true, 0));
      }
    }
    if (doVZERO) {
      if (types.Contains("std") || types.Contains("all")) {
        if (!onlyTr)
	  CoupleCar(mac, Form(args.Data(), "VZERO", false, false, 0));
      	if (tpcTr)
	  CoupleCar(mac, Form(args.Data(), "VZERO", false, false, 1));
      	if (hybridTr)
	  CoupleCar(mac, Form(args.Data(), "VZERO", false, false, 2));
      }
      if (types.Contains("eta-gap") || types.Contains("all")) {
        if (!onlyTr)
	  CoupleCar(mac, Form(args.Data(), "VZERO", true, false, 0));
      	if (tpcTr)
      	  CoupleCar(mac, Form(args.Data(), "VZERO", true, false, 1));
      	if (hybridTr)
      	  CoupleCar(mac, Form(args.Data(), "VZERO", true, false, 2));
      }
      if (types.Contains("3cor") || types.Contains("all")) {
        if (!onlyTr)
	  CoupleCar(mac, Form(args.Data(), "VZERO", false, true, 0));
      }
    }
  }
  //__________________________________________________________________
  /**
   * Add the EP task(s)
   */
  void AddEPTasks()
  {
    // --- Get the parameters ----------------------------------------
    Int_t    moment   = fOptions.AsInt("max-mom");
    TString  etaGap   = fOptions.Get("eta-gap");
    TString  addFlow  = fOptions.Get("afterburner");
    Bool_t   mc       = fOptions.AsBool("mc"); 
    TString  fwdDets  = fOptions.Get("fwd-dets");

    // --- Add the task ----------------------------------------------
    CoupleCar("AddTaskEventplane.C","");
    CoupleCar("AddTaskForwardFlowEP.C", 
	      Form("%d, %d, \"%s\"", mc, moment, fwdDets.Data()));
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
  /** 
   * Create MC input handler.  Since this train does not use the MC
   * information from @c galice.root, @c Kinematics.root, and @c
   * TrackRefs.root directly, we define this member function to return
   * null even when doing MC analysis.  This train _only_ looks at the
   * AOD object of the _processed_ MC data.
   * 
   * @return Always 0
   */
  AliVEventHandler* CreateMCHandler(UShort_t, bool)
  { 
    return 0;
  }
  //__________________________________________________________________
  /** 
   * Get the class name of the train setup 
   * 
   * @return Class name 
   */
  const char* ClassName() const { return "MakeFlowTrain"; }
  //__________________________________________________________________
};
//
// EOF
//
