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
    fOptions.Add("dNdeta",  "Add dN/deta tasks");
    fOptions.Add("flow",    "Add flow tasks");
    fOptions.Add("cent",    "Use centrality");
    fOptions.Add("sat-vtx", "Use satellite interactions");
    // ESD settings
    fOptions.Add("run",   "NUMBER",  "Run number", 0);
    fOptions.Add("sys",   "SYSTEM",  "1:pp, 2:PbPb, 3:pPb", "");
    fOptions.Add("snn",   "ENERGY",  "Center of mass energy in GeV", "");
    fOptions.Add("field", "STRENGTH","L3 field strength in kG", "");
    fOptions.Add("tpc-ep", "Use TPC event plane");
    fOptions.Add("forward-config", "FILE", "Forward configuration", 
		 "ForwardAODConfig.C");
    fOptions.Add("central-config", "FILE", "Forward configuration", 
		 "CentralAODConfig.C");
    fOptions.Add("corr", "DIR", "Corrections dir", "");
    // dNdeta AOD settings
    fOptions.Add("dndeta-config", "FILE", "dN/deta configuration", 
		 "dNdetaConfig.C");
    fOptions.Add("trig",     "TYPE",       "Trigger type",             "INEL");
    fOptions.Add("vzMin",    "CENTIMETER", "Min Ip Z",                 "-10");
    fOptions.Add("vzMax",    "CENTIMETER", "Max Ip Z",                 "+10");
    fOptions.Add("scheme",   "SCHEME",     "Normalization scheme",     "");
    fOptions.Add("trigEff",  "EFFICENCY",  "Trigger effeciency",       "1");
    fOptions.Add("trigEff0", "EFFICENCY",  "0-bin trigger effeciency", "1");
    // Flow AOD settings
    fOptions.Add("mom",     "MOMENTS",      "Flow moments to analyse", "234");
    fOptions.Add("eta-gap", "[no,yes,both]","Whether to use an eta-gap","both");
    fOptions.Add("eg-value",   "DISTANCE",  "Size of eta gap",          "2.0");
    fOptions.Add("b-cent",     "Use the impact parameter for centrality");
    fOptions.Add("afterburner","[eta,phi,b,pid]", "What to afterburn", "");
    fOptions.Add("ab-type",    "[1-4]",     "Type of afterburner",     "");
    fOptions.Add("ab-order",    "[2-6]",     "Order of afterburner",    "");
    fOptions.Add("outlier-fmd", "NSIGMA", "Outlier cut for FMD", "4.0");
    fOptions.Add("outlier-spd", "NSIGMA", "Outlier cut for SPD", "0.0");

    fOptions.Set("type", "ESD"); 
    fOptions.Show(std::cout);

  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward.root");

    // --- Load libraries/pars ---------------------------------------
    fHelper->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;
   
    // --- Add TPC eventplane task -----------------------------------
    if (fOptions.Has("tpc-ep")) gROOT->Macro("AddTaskEventplane.C");

    // --- Task to copy header information ---------------------------
    gROOT->Macro("AddTaskCopyHeader.C");

    // --- Get options -----------------------------------------------
    ULong_t  run    = fOptions.AsInt("run", 0);
    UShort_t sys    = fOptions.AsInt("sys", 0);
    UShort_t sNN    = fOptions.AsInt("snn", 0);
    UShort_t fld    = fOptions.AsInt("field", 0);
    Bool_t   cent   = fOptions.Has("cent");
    Bool_t   satVtx = fOptions.AsBool("sat-vtx");
    TString  cor = "";
    if (fOptions.Has("corr")) cor = fOptions.Get("corr");
    
    // --- Add the task ----------------------------------------------
    TString fwdConfig = fOptions.Get("forward-config");
    gROOT->Macro(Form("AddTaskForwardMult.C(%d,%d,%d,%d,%d,\"%s\")", 
		      mc, run, sys, sNN, fld, fwdConfig.Data()));
    fHelper->LoadAux(gSystem->Which(gROOT->GetMacroPath(), fwdConfig));

    // --- Add the task ----------------------------------------------
    TString cenConfig = fOptions.Get("central-config");
    gROOT->Macro(Form("AddTaskCentralMult.C(%d,%d,%d,%d,%d,\"%s\")", 
		      mc, run, sys, sNN, fld, cenConfig.Data()));
    fHelper->LoadAux(gSystem->Which(gROOT->GetMacroPath(), cenConfig));

    // --- Add MC particle task --------------------------------------
    if (mc) gROOT->Macro("AddTaskMCParticleFilter.C");

    
    // --- Add dN/deta tasks -----------------------------------------
    if (fOptions.Has("dNdeta")) {
      // --- Get parameters --------------------------------------------
      TString  trig   = fOptions.Get("trig");
      TString  scheme = fOptions.Get("scheme");
      Double_t vzMin  = fOptions.AsDouble("vzmin", -10);
      Double_t vzMax  = fOptions.AsDouble("vzmax", +10);
      Double_t effT   = fOptions.AsDouble("trigEff", 1);
      Double_t effT0  = fOptions.AsDouble("trigEff0", 1);
      TString  config = fOptions.Get("dndeta-config");

      // --- Form arguments --------------------------------------------
      TString args;
      args.Form("\"%s\",\"%s\",%f,%f,%d,\"%s\",%g,%g",
		config.Data(),trig.Data(), vzMin, vzMax, cent, scheme.Data(),
		effT, effT0);
      // --- Add the task ----------------------------------------------
      gROOT->Macro(Form("AddTaskForwarddNdeta.C(%s);", args.Data()));
      gROOT->Macro(Form("AddTaskCentraldNdeta.C(%s);", args.Data()));
      gROOT->Macro(Form("AddTaskMCTruthdNdeta.C(%s);", args.Data()));
    }
    
    // --- Add the flow task -----------------------------------------
    if (fOptions.Has("flow")) {
      // --- Get the parameters ----------------------------------------
      TString  type    = fOptions.Get("mom");
      TString  etaGap  = fOptions.Get("eta-gap");
      Double_t egValue = fOptions.AsDouble("eg-value");
      Bool_t   useCent = fOptions.AsBool("b-cent");
      TString  addFlow = fOptions.Get("afterburner");
      Int_t    abType  = fOptions.AsInt("ab-type");
      Int_t    abOrder = fOptions.AsInt("ab-order");
      Double_t fmdCut  = fOptions.AsDouble("outlier-fmd");
      Double_t spdCut  = fOptions.AsDouble("outlier-spd");

      // --- Add the task ----------------------------------------------
      // Notice the '%%d' for second argument.  We substitute that later
      TString arg(Form("AddTaskForwardFlow.C(\"%s\",%%d,%d,%f,%f,%f,%d,%d,\"%s\",%d,%d)",
		       type.Data(),
		       mc, 
		       fmdCut, 
		       spdCut,
		       egValue,
		       useCent,
		       satVtx, 
		       addFlow.Data(), 
		       abType, 
		       abOrder));
      
      if (etaGap.Contains("no") || etaGap.Contains("false") ||
	  etaGap.Contains("both"))
	gROOT->Macro(Form(arg.Data(), false));
      
      if (etaGap.Contains("yes") || etaGap.Contains("true") ||
	  etaGap.Contains("both"))
	gROOT->Macro(Form(arg.Data(), true));
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
  void CreateCentralitySelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    if (!fOptions.Has("cent")) return;
    TrainSetup::CreateCentralitySelection(mc, mgr);
  }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeFullTrain"; }
};
//
// EOF
//
