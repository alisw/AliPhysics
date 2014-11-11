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
    fOptions.Add("dndeta",  "Add dN/deta tasks");
    fOptions.Add("flow",    "Add flow tasks");
    fOptions.Add("cent",    "Use centrality");
    fOptions.Add("sat-vtx", "Use satellite interactions");
    fOptions.Add("vzMin",    "CENTIMETER", "Min Ip Z",                 "-10");
    fOptions.Add("vzMax",    "CENTIMETER", "Max Ip Z",                 "+10");
    // ESD settings
    fOptions.Add("run",            "NUMBER",  "Run number", 0);
    fOptions.Add("sys",            "SYSTEM",  "1:pp, 2:PbPb, 3:pPb", "");
    fOptions.Add("snn",            "ENERGY",  "Center of mass energy GeV","");
    fOptions.Add("field",          "STRENGTH","L3 field strength in kG", "");
    fOptions.Add("aod-forward-config", "FILE", "Forward configuration", 
		 "ForwardAODConfig.C");
    fOptions.Add("aod-central-config", "FILE", "Forward configuration", 
		 "CentralAODConfig.C");
    fOptions.Add("aod-corr", "DIR",    "Corrections dir", "");
    // dNdeta AOD settings
    fOptions.Add("dndeta-config",  "FILE", "dN/deta configuration", 
		 "dNdetaConfig.C");
    fOptions.Add("dndeta-trig",    "TYPE",     "Trigger type",         "INEL");
    fOptions.Add("dndeta-scheme",  "SCHEME",   "Normalization scheme",     "");
    fOptions.Add("dndeta-trigEff", "EFFICENCY","Trigger effeciency",       1);
    fOptions.Add("dndeta-trigEff0","EFFICENCY","0-bin trigger effeciency", 1);
    // Flow AOD settings
    fOptions.Add("flow-max-mom", "2|3|4|5", "Max flow moment to analyse", "5");
    fOptions.Add("flow-detectors", "[fmd,vzero]", 
		 "Forward detectors", "fmd+vzero");
    fOptions.Add("flow-qc-types", "[std,eta-gap,3cor,all]", 
		 "Which types of QC's to do", "all");
    fOptions.Add("flow-eta-gap","DISTANCE",  "Size of eta gap",          "2.0");
    fOptions.Add("flow-b-cent", "Use the impact parameter for centrality");
    fOptions.Add("flow-afterburner","[eta,phi,b,pid]", "What to afterburn", "");
    fOptions.Add("flow-outlier-fmd", "NSIGMA", "Outlier cut for FMD", "4.0");
    fOptions.Add("flow-outlier-spd", "NSIGMA", "Outlier cut for SPD", "0.0");
    fOptions.Add("flow-mc-vtx", "Whether to get the vertex from the MC header");
    fOptions.Add("flow-ref-tracks", "[tpc,hybrid,only]", 
		 "Whether or only to use tracks for reference flow", 
		 "tpc+hybrid");
    fOptions.Add("flow-ep", "Add Event Plane tasks (need VZERO AOD objects)");

    fOptions.Set("type", "ESD"); 
    fOptions.Show(std::cout);

  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward.root");

    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = HasMCHandler();
   
    // --- Task to copy header information ---------------------------
    CoupleCar("AddTaskCopyHeader.C", "");

    // --- Get options -----------------------------------------------
    ULong_t  run    = fOptions.AsInt   ("run", 0);
    UShort_t sys    = fOptions.AsInt   ("sys", 0);
    UShort_t sNN    = fOptions.AsInt   ("snn", 0);
    UShort_t fld    = fOptions.AsInt   ("field", 0);
    Bool_t   cent   = fOptions.Has     ("cent");
    Bool_t   satVtx = fOptions.AsBool  ("sat-vtx");
    Double_t vzMin  = fOptions.AsDouble("vzmin", -10);
    Double_t vzMax  = fOptions.AsDouble("vzmax", +10);
    TString  cor    = "";
    if (fOptions.Has("aod-corr")) cor = fOptions.Get("aod-corr");
    
    // --- Add the task ----------------------------------------------
    TString fwdConfig = fOptions.Get("aod-forward-config");
    AliAnalysisTask* fwd = CoupleCar("AddTaskForwardMult.C",
				     Form("%d,%lu,%hu,%hu,%hd,\"%s\",\"%s\"", 
					  mc, run, sys, sNN, fld, 
					  fwdConfig.Data(), cor.Data()));
    fRailway->LoadAux(gSystem->Which(gROOT->GetMacroPath(), fwdConfig));

    // --- Add the task ----------------------------------------------
    TString cenConfig = fOptions.Get("aod-central-config");
    AliAnalysisTask* cen = CoupleCar("AddTaskCentralMult.C",
				     Form("%d,%lu,%hu,%hu,%hd,\"%s\",\"%s\"", 
					  mc, run, sys, sNN, fld, 
					  cenConfig.Data(),cor.Data()));
    fRailway->LoadAux(gSystem->Which(gROOT->GetMacroPath(), cenConfig));
    if (!cor.IsNull()) {
      if (fwd) 
	fRailway->LoadAux(Form("%s/fmd_corrections.root",cor.Data()), true);
      if (cen) 
	fRailway->LoadAux(Form("%s/spd_corrections.root",cor.Data()), true);
    }

    // --- Add MC particle task --------------------------------------
    if (mc) CoupleCar("AddTaskMCParticleFilter.C","");

    
    // --- Add dN/deta tasks -----------------------------------------
    if (fOptions.Has("dNdeta")) {
      // --- Get parameters --------------------------------------------
      TString  trig   = fOptions.Get     ("dndeta-trig");
      TString  scheme = fOptions.Get     ("dndeta-scheme");
      Double_t effT   = fOptions.AsDouble("dndeta-trigEff", 1);
      Double_t effT0  = fOptions.AsDouble("dndeta-trigEff0", 1);
      TString  config = fOptions.Get     ("dndeta-config");

      // --- Form arguments --------------------------------------------
      TString args;
      args.Form("\"%s\",\"%s\",%f,%f,%d,\"%s\",%g,%g",
		config.Data(),trig.Data(), vzMin, vzMax, cent, scheme.Data(),
		effT, effT0);
      // --- Add the task ----------------------------------------------
      CoupleCar("AddTaskForwarddNdeta.C", args);
      CoupleCar("AddTaskCentraldNdeta.C", args);
      CoupleCar("AddTaskMCTruthdNdeta.C", args);
    }
    
    // --- Add the flow task -----------------------------------------
    if (fOptions.Has("flow")) {
      // --- Get the parameters ----------------------------------------
      Int_t    moment   = fOptions.AsInt   ("flow-max-mom");
      TString  fwdDets  = fOptions.Get     ("flow-detectors");
      TString  types    = fOptions.Get     ("flow-qc-types");
      Double_t egValue  = fOptions.AsDouble("flow-eta-gap");
      TString  tracks   = fOptions.Get     ("flow-ref-tracks");
      Bool_t   useCent  = fOptions.AsBool  ("flow-b-cent");
      Bool_t   useMCVtx = fOptions.AsBool  ("flow-mc-vtx");
      Bool_t   addFlow  = fOptions.AsBool  ("flow-afterburner");
      Double_t fmdCut   = fOptions.AsDouble("flow-outlier-fmd");
      Double_t spdCut   = fOptions.AsDouble("flow-outlier-spd");

      types.ToLower();
      fwdDets.ToUpper();
      Bool_t doFMD      = fwdDets.Contains("FMD");
      Bool_t doVZERO    = fwdDets.Contains("VZERO");
      tracks.ToLower();
      Bool_t onlyTr     = tracks.Contains("only");
      Bool_t tpcTr      = tracks.Contains("tpc");
      Bool_t hybridTr   = tracks.Contains("hybrid");
      Bool_t ispA       = (sys == 3 ? kTRUE : kFALSE);
      
      // Notice the place holders at arg=2,3,4, and 9, These are
      // 
      //   2: Detector to use (FMD/VZERO)
      //   3: Whether to use eta gap (true/false)
      //   4: Do 3-particle correlations (true/false)
      //   9: Use tracks for referernce flow (true/false)
      TString args;
      args=TString::Format("%d,\"%%s\",%%d,%%d,%d,%f,%f,%f,%%d,%d,%d,%d,%d,%d",
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

      // --- Add the task ----------------------------------------------
      if (fOptions.Has("flow-ep")) {
	CoupleCar("AddTaskEventplane.C", "");
	CoupleCar("AddTaskForwardFlowEP.C", 
		  Form("%d, %d, \"%s\"", mc, moment, fwdDets.Data()));
      }
    }
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

    // --- Special for pPb pilot run Sep. 2012 -----------------------
    UShort_t sys = fOptions.AsInt("sys", 0);
    if (sys == 3) { 
      Warning("CreatePhysicsSelection", 
	      "Special setup for pPb pilot run September, 2012");
      gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
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
   * @param mgr Manager
   */
  void CreateCentralitySelection(Bool_t mc)
  {
    if (!fOptions.Has("cent")) return;
    TrainSetup::CreateCentralitySelection(mc);
  }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeFullTrain"; }
};
//
// EOF
//
