#include "TrainSetup.C"
/**
 * @file   MakeFullTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:53:43 2012
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_trains_specific
 */

//====================================================================
/**
 * Analysis train to make Forward and Central multiplicity, @f$
 * dN/d\eta@f$, flow and @f$\Psi_R@f$ in one loop over the ESDs 
 *
 * @ingroup pwglf_forward_aod
 * @ingroup pwglf_forward_dndete
 * @ingroup pwglf_forward_flow
 * @ingroup pwglf_forward_trains_specific
 */
class MakeFullTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train (free form)
   */
  MakeFullTrain(const  char* name) 
    : TrainSetup(name)
  {
    // General
    fOptions.Add("cent",          "ESTIMATOR", "Use centrality", "none");
    fOptions.Add("cent-bins",     "EDGES",     "Centrality bin edges", "");
    fOptions.Add("ipz-min",       "CENTIMETER","Min Ip Z",       "-10");
    fOptions.Add("ipz-max",       "CENTIMETER","Max Ip Z",        "+10");
    fOptions.Add("trigger",       "TYPE",      "Trigger to use", "INEL");
    fOptions.Add("filter",        "FILTER",    "Vetos",          "");
    fOptions.Add("satellite",                  "Use satellite interactions");
    fOptions.Add("sys",           "SYSTEM",    "1:pp, 2:PbPb, 3:pPb", "");
    // ESD settings
    fOptions.Add("aod-run",        "NUMBER",  "Run number", 0);
    fOptions.Add("aod-snn",        "ENERGY",  "Center of mass energy GeV","");
    fOptions.Add("aod-field",      "STRENGTH","L3 field strength in kG", "");
    fOptions.Add("aod-forward-config", "FILE", "Forward configuration", 
		 "ForwardAODConfig.C");
    fOptions.Add("aod-central-config", "FILE", "Forward configuration", 
		 "CentralAODConfig.C");
    fOptions.Add("aod-corr", "DIR",    "Corrections dir", "");
    // dNdeta AOD settings
    fOptions.Add("dndeta",                     "Add dN/deta tasks");
    fOptions.Add("dndeta-config",  "FILE",     "dN/deta configuration", 
		 "dNdetaConfig.C");
    fOptions.Add("dndeta-trig",    "TYPE",     "Trigger type",         "INEL");
    fOptions.Add("dndeta-scheme",  "SCHEME",   "Normalization scheme",     "");
    fOptions.Add("dndeta-trigEff", "EFFICENCY","Trigger effeciency",       1);
    fOptions.Add("dndeta-trigEff0","EFFICENCY","0-bin trigger effeciency", 1);
    // Multiplicity setttings
    fOptions.Add("pnch", "Run the P(Nch) task", false);
    fOptions.Add("pnch-resp", "Run the response matrix task for P(Nch)", false);
    fOptions.Add("pnch-bias", "Run the trigger bias task for P(Nch)", false);
    // Flow AOD settings
    fOptions.Add("flow",                     "Add flow tasks");
    fOptions.Add("flow-max-mom",    "MOMENT","Max flow moment to analyse", "5");
    fOptions.Add("flow-detectors",  "NAMES", "Forward detectors", "fmd+vzero");
    fOptions.Add("flow-qc-types",   "TYPES", 
		 "Which types of QC's to do [std,eta-gap,3cor,all]", "all");
    fOptions.Add("flow-eta-gap",    "DISTANCE","Size of eta gap",     "2.0");
    fOptions.Add("flow-use-b",      "Use impact param. for centrality",false);
    fOptions.Add("flow-afterburner","WHAT",
		 "What to afterburn [eta,phi,b,pid]", "");
    fOptions.Add("flow-outlier-fmd","NSIGMA", "Outlier cut for FMD", "4.0");
    fOptions.Add("flow-outlier-spd","NSIGMA", "Outlier cut for SPD", "0.0");
    fOptions.Add("flow-mc-vtx",
		 "Whether to get the vertex from the MC header");
    fOptions.Add("flow-ref-tracks", "TYPE", 
		 "Whether or only to use tracks for reference flow "
		 "[tpc,hybrid,only]", "tpc+hybrid");
    fOptions.Add("flow-ep", "Add Event Plane tasks (need VZERO AOD objects)");

    fOptions.Set("type", "ESD"); 
    fOptions.Show(std::cout);
  }
protected:
  /**
   * Create tasks to make Nch per (eta,phi) per event in AOD 
   */
  Bool_t CreateAODTasks()
  {
    // --- Get options -----------------------------------------------
    ULong_t  run    = fOptions.AsInt   ("aod-run",   0);
    UShort_t sNN    = fOptions.AsInt   ("aod-snn",   0);
    UShort_t fld    = fOptions.AsInt   ("aod-field", 0);
    UShort_t sys    = fOptions.AsInt   ("sys",       0);
    Bool_t   satVtx = fOptions.AsBool  ("satellite");
    Bool_t   mc     = HasMCHandler();
    TString  cor    = "";
    if (fOptions.Has("aod-corr")) cor = fOptions.Get("aod-corr");
    
    // --- Add forward task ------------------------------------------
    TString fwdConfig = fOptions.Get("aod-forward-config");
    AliAnalysisTask* fwd = CoupleCar("AddTaskForwardMult.C",
				     Form("%d,%lu,%hu,%hu,%hd,\"%s\",\"%s\"", 
					  mc, run, sys, sNN, fld, 
					  fwdConfig.Data(), cor.Data()));
    if (!fwd) return false;
    fRailway->LoadAux(gSystem->Which(gROOT->GetMacroPath(), fwdConfig));

    // --- Add central task ------------------------------------------
    TString cenConfig = fOptions.Get("aod-central-config");
    AliAnalysisTask* cen = CoupleCar("AddTaskCentralMult.C",
				     Form("%d,%lu,%hu,%hu,%hd,\"%s\",\"%s\"", 
					  mc, run, sys, sNN, fld, 
					  cenConfig.Data(),cor.Data()));
    if (!cen) return false;
    fRailway->LoadAux(gSystem->Which(gROOT->GetMacroPath(), cenConfig));
    if (!cor.IsNull()) {
      if (fwd) 
	fRailway->LoadAux(Form("%s/fmd_corrections.root",cor.Data()), true);
      if (cen) 
	fRailway->LoadAux(Form("%s/spd_corrections.root",cor.Data()), true);
    }
    return true;
  }
  /** 
   * Create a dN/deta task 
   * 
   * @param which Which kind 
   * 
   * @return true on success
   */
  Bool_t CreatedNdetaTask(const char* which)
  {
    // --- Get parameters ------------------------------------------
    TString  config = fOptions.Get     ("dndeta-config");
    TString  args   = Form("\"%s\",\"%s\"", which, config.Data());
    UInt_t   mask   = AliVEvent::kAny;
    
    AliAnalysisTaskSE* tsk = CoupleSECar("AddTaskdNdeta.C",args,mask);
    if (!tsk) {
      Printf("Failed to add task via AddTaskdNdeta.C(%s,%s)",
	     which, config.Data());
      return false;
    }
    FromOption(tsk, "TriggerMask",         "trigger",        "INEL");
    FromOption(tsk, "FilterMask",          "filter",      "OUTLIER|PILEUP-BIN");
    FromOption(tsk, "CentralityMethod",    "cent",           "");
    FromOption(tsk, "CentralityAxis",      "cent-bins",      "default");
    FromOption(tsk, "IpZMin",              "ipz-min",        -10.);
    FromOption(tsk, "IpZMax",              "ipz-max",        +10.);
    FromOption(tsk, "SatelliteVertices",   "satellite",      false);
    FromOption(tsk, "NormalizationScheme", "dndeta-scheme",  "EVENT,TRIGGER");
    FromOption(tsk, "TriggerEff",          "dndeta-trigEff", 1.);
    FromOption(tsk, "TriggerEff0",         "dndeta-trigEff0",1.);
    return true;
  }
  /** 
   * Create a flow task 
   * 
   * @param fmt        Arguments format 
   * @param det        Detector to make task for 
   * @param useEtaGap  Whether to do eta-gaps 
   * @param use3cor    Whether to do 3 particle correlations 
   * @param tracks     Which kinds of tracks to do
   * 
   * @return true on success
   */
  Bool_t CreateFlowTask(const char* fmt,
			const char* det,
			Bool_t      useEtaGap,
			Bool_t      use3cor,
			Short_t     tracks)
  {
    const char* mac = "AddTaskForwardFlowQC.C";
    AliAnalysisTaskSE* task =
      CoupleSECar(mac,Form(fmt, det, useEtaGap, use3cor, tracks));
    if (!task) return false;
    FromOption(task, "CentralityAxis",      "cent-bins",      "default");
    return true;
  }
  /** 
   * Create all the flow tasks 
   * 
   * @return true on success 
   */
  Bool_t CreateFlowTasks()
  {
    // --- Get the parameters ----------------------------------------
    Bool_t   mc       = HasMCHandler();
    Int_t    sys      = fOptions.AsInt   ("sys");
    Int_t    moment   = fOptions.AsInt   ("flow-max-mom");
    TString  fwdDets  = fOptions.Get     ("flow-detectors");
    TString  types    = fOptions.Get     ("flow-qc-types");
    Double_t egValue  = fOptions.AsDouble("flow-eta-gap");
    TString  tracks   = fOptions.Get     ("flow-ref-tracks");
    Bool_t   useB     = fOptions.AsBool  ("flow-use-b");
    Bool_t   useMCVtx = fOptions.AsBool  ("flow-mc-vtx");
    Bool_t   addFlow  = fOptions.AsBool  ("flow-afterburner");
    Double_t fmdCut   = fOptions.AsDouble("flow-outlier-fmd");
    Double_t spdCut   = fOptions.AsDouble("flow-outlier-spd");
    Bool_t   satVtx   = fOptions.AsBool  ("satellite");
    
    fwdDets.ToUpper();
    Bool_t doFMD      = fwdDets.Contains("FMD");
    Bool_t doVZERO    = fwdDets.Contains("VZERO");
    tracks.ToLower();
    Bool_t onlyTr     = tracks.Contains("only");
    Bool_t tpcTr      = tracks.Contains("tpc");
    Bool_t hybridTr   = tracks.Contains("hybrid");
    Bool_t ispA       = sys == 3 || sys == 4;
    types.ToLower();
    Bool_t std        = types.Contains("std")     || types.Contains("all");
    Bool_t etaGap     = types.Contains("eta-gap") || types.Contains("all");
    Bool_t cor3       = types.Contains("3cor")    || types.Contains("all");
    
    // Notice the place holders at arg=2,3,4, and 9, These are
    // 
    //   2: Detector to use (FMD/VZERO)
    //   3: Whether to use eta gap (true/false)
    //   4: Do 3-subevent correlations (true/false)
      //   9: Use tracks for referernce flow (true/false)
    TString args;
    args=TString::Format("%d,\"%%s\",%%d,%%d,%d,%f,%f,%f,%%d,%d,%d,%d,%d,%d",
			 moment,
			 mc, 
			 fmdCut, 
			 spdCut,
			 egValue,
			 !useB,
			 ispA,
			 useMCVtx,
			 satVtx, 
			 addFlow);
    
    // --- Add the task ----------------------------------------------
    if (doFMD) {
      if (std) {
	if (!onlyTr &&!CreateFlowTask(args,  "FMD",false,false,0)) return false;
	if (tpcTr   &&!CreateFlowTask(args,  "FMD",false,false,1)) return false;
	if (hybridTr&&!CreateFlowTask(args,  "FMD",false,false,2)) return false;
      }
      if (etaGap) {
	if (!onlyTr &&!CreateFlowTask(args,  "FMD", true,false,0)) return false;
	if (tpcTr   &&!CreateFlowTask(args,  "FMD", true,false,1)) return false;
	if (hybridTr&&!CreateFlowTask(args,  "FMD", true,false,2)) return false;
      }
      if (cor3) {
	if (!onlyTr &&!CreateFlowTask(args,  "FMD",false,true, 0)) return false;
      }
    }
    if (doVZERO) {
      if (std) {
	if (!onlyTr &&!CreateFlowTask(args,"VZERO",false,false,0)) return false;
	if (tpcTr   &&!CreateFlowTask(args,"VZERO",false,false,1)) return false;
	if (hybridTr&&!CreateFlowTask(args,"VZERO",false,false,2)) return false;
      }
      if (etaGap) {
	if (!onlyTr &&!CreateFlowTask(args,"VZERO", true,false,0)) return false;
	if (tpcTr   &&!CreateFlowTask(args,"VZERO", true,false,1)) return false;
	if (hybridTr&&!CreateFlowTask(args,"VZERO", true,false,2)) return false;
      }
      if (cor3) {
	if (!onlyTr &&!CreateFlowTask(args,"VZERO",false, true,0)) return false;
      }
    }
    
    // --- Add the task ----------------------------------------------
    if (fOptions.Has("flow-ep")) {
      CoupleCar("AddTaskEventplane.C", "");
      CoupleCar("AddTaskForwardFlowEP.C", 
		  Form("%d, %d, \"%s\"", mc, moment, fwdDets.Data()));
    }
  }    
  /** 
   * Create P(Nch) task 
   * 
   * @return true on success
   */
  Bool_t CreatePNchTask(const char* mac)
  {
    Info("CreatePNchTask","\n"
	 "*****************************************************\n"
	 "Creating a task via %s\n"
	 "*****************************************************", mac);
    AliAnalysisTaskSE* tsk =  CoupleSECar(mac);
    if (!tsk) return false;
    FromOption(tsk, "TriggerMask",      "trigger",     "INEL");
    FromOption(tsk, "FilterMask",       "filter",      "OUTLIER|PILEUP-BIN");
    // FromOption(tsk, "CentralityMethod", "cent",        "");
    FromOption(tsk, "CentralityAxis",   "cent-bins",   "default");
    FromOption(tsk, "IpZMin",           "ipz-min",     -10.);
    FromOption(tsk, "IpZMax",           "ipz-max",     +10.);
    return true;
  }
  /** 
   * Create the tasks 
   * 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward.root");

    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = HasMCHandler();

    // --- Task to copy header information ---------------------------
    CoupleCar("AddTaskCopyHeader.C", "");


    // --- Add MC particle task --------------------------------------
    if (mc) CoupleCar("AddTaskMCParticleFilter.C","");

    // --- Create AOD tasks ------------------------------------------
    if (!CreateAODTasks()) return;
    
    // --- Add dN/deta tasks -----------------------------------------
    if (fOptions.Has("dndeta")) {
      if (!CreatedNdetaTask("Forward")) return;
      if (!CreatedNdetaTask("Central")) return;
      if (mc && !CreatedNdetaTask("MCTruth")) return;
    }
    
    // --- P(Nch) task -----------------------------------------------
    if (fOptions.Has("pnch") &&
	!CreatePNchTask("AddTaskMultDistributions.C")) return;

    // --- P(Nch) response matrix task -------------------------------
    if (mc && fOptions.Has("pnch-resp") &&
	!CreatePNchTask("AddTaskCreateRespMatr.C")) return;

    // --- P(Nch) trigger bias task ----------------------------------
    if (mc && fOptions.Has("pnch-bias") &&
	!CreatePNchTask("AddTaskTriggerCorrection.C")) return;
    
    // --- Add the flow task -----------------------------------------
    if (fOptions.Has("flow") && !CreateFlowTasks()) return;
  }
  // The code below is redundant and shouldn't be used. We keep it
  // here for reference.  If one needs a special physics selection,
  // one can use the option "ps" defined by the TrainSetup.C class.
#if 0
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

    // --- Special for pPb pilot run Sep. 2012 -----------------------
    UShort_t sys = fOptions.AsInt("sys", 0);
    if (sys == 3) { 
      Warning("CreatePhysicsSelection", 
	      "Special setup for pPb pilot run September, 2012");
      gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/ANALYSIS/macros",
			       gROOT->GetMacroPath()));
      gROOT->LoadMacro("PhysicsSelectionOADB_CINT5_pA.C");
      gROOT->ProcessLine(Form("((AliPhysicsSelection*)%p)"
			      "->SetCustomOADBObjects("
			      "OADBSelection_CINT5_V0A(),0);", ps));
      ps->SetSkipTriggerClassSelection(true);
    }
    // --- Ignore trigger class when selecting events.  This means ---
    // --- that we get offline+(A,C,E) events too --------------------
    // ps->SetSkipTriggerClassSelection(true);
  }
  //__________________________________________________________________
  /** 
   * Create the centrality selection only if requested
   * 
   * @param mc  Monte-Carlo truth flag 
   */
  void CreateCentralitySelection(Bool_t mc)
  {
    if (!fOptions.Has("cent")) return;
    TrainSetup::CreateCentralitySelection(mc);
  }
#endif
  //__________________________________________________________________
  const char* ClassName() const { return "MakeFullTrain"; }
};
//
// EOF
//
