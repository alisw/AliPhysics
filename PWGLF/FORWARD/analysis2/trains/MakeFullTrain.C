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
    fOptions.Add("sys",   "SYSTEM",  "1:pp, 2:PbPb, 3:pPb", "");
    fOptions.Add("snn",   "ENERGY",  "Center of mass energy in GeV", "");
    fOptions.Add("field", "STRENGTH","L3 field strength in kG", "");
    fOptions.Add("cent",  "Use centrality");
    fOptions.Add("tpc-ep", "Use TPC event plane");
    fOptions.Add("forward-config", "FILE", "Forward configuration", 
		 "ForwardAODConfig.C");
    fOptions.Add("central-config", "FILE", "Forward configuration", 
		 "CentralAODConfig.C");
    fOptions.Add("satelitte", "Use satelitte interactions");
    fOptions.Add("trig",     "TYPE", "Trigger type", "INEL");
    fOptions.Add("vzMin",    "CENTIMETER", "Min Ip Z", "-10");
    fOptions.Add("vzMax",    "CENTIMETER", "Max Ip Z", "+10");
    fOptions.Add("scheme",   "SCHEME", "Normalization scheme", "");
    fOptions.Add("trigEff",  "EFFICENCY", "Trigger effeciency", "1");
    fOptions.Add("trigEff0", "EFFICENCY", "0-bin trigger effeciency", "1");
    fOptions.Add("cut-edges", "Cut acceptance edges");
    fOptions.Add("mom", "MOMENTS", "Flow moments to determine", "123456");
    fOptions.Add("outlier-fmd", "NSIGMA", "Cut on outliers in FMD", "4.1");
    fOptions.Add("outlier-spd", "NSIGMA", "Cut on outliers in SPD", "4.1");
    fOptions.Add("afterburner", "WHAT", "What to afterburn", "eta phi b pid");
    fOptions.Add("ab-type", "1|2|3|4", "Type of afterburner", "");
    fOptions.Add("ab-order", "1|2|3|4|5|6", "Order used by afterburner", "");

    fOptions.Set("type", "ESD");
  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param par  Whether to use par files 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2");
    
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
    UShort_t sys = fOptions.AsInt("sys", 0);
    UShort_t sNN = fOptions.AsInt("snn", 0);
    UShort_t fld = fOptions.AsInt("field", 0);
    
    // --- Add the task ----------------------------------------------
    TString fwdConfig = fOptions.Get("forward-config");
    gROOT->Macro(Form("AddTaskForwardMult.C(%d,%d,%d,%d,\"%s\")", 
		      mc, sys, sNN, fld, fwdConfig.Data()));
    fHelper->LoadAux(gSystem->Which(gROOT->GetMacroPath(), fwdConfig));

    // --- Add the task ----------------------------------------------
    TString cenConfig = fOptions.Get("central-config");
    gROOT->Macro(Form("AddTaskCentralMult.C(%d,%d,%d,%d,\"%s\")", 
		      mc, sys, sNN, fld, cenConfig.Data()));
    fHelper->LoadAux(gSystem->Which(gROOT->GetMacroPath(), cenConfig));

    // --- Get parameters --------------------------------------------
    TString  trig   = fOptions.Get("trig");
    TString  scheme = fOptions.Get("scheme");
    Double_t vzMin  = fOptions.AsDouble("vzmin", -10);
    Double_t vzMax  = fOptions.AsDouble("vzmax", +10);
    Double_t effT   = fOptions.AsDouble("trigEff", 1);
    Double_t effT0  = fOptions.AsDouble("trigEff0", 1);
    Bool_t   cent   = fOptions.Has("cent");
    Bool_t   edges  = fOptions.Has("cut-edges");

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwarddNdeta.C(\"%s\",%f,%f,%d,\"%s\")",
		      trig.Data(), vzMin, vzMax, cent, scheme.Data(),
		      edges, trigEff, trigEff0));

    gROOT->Macro(Form("AddTaskCentraldNdeta.C(\"%s\",%f,%f,%d,\"%s\")",
		      trig.Data(), vzMin, vzMax, cent, scheme.Data()),
		      edges, trigEff, trigEff0);

    gROOT->Macro(Form("AddTaskMCTruthdNdeta.C(\"%s\",%f,%f,%d,\"%s\")",
		      trig.Data(), vzMin, vzMax, cent, scheme.Data()),
		      edges, trigEff, trigEff0);

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskFMDEventPlane.C(%d)", mc));

    // --- Add MC particle task --------------------------------------
    if (mc) gROOT->Macro("AddTaskMCParticleFilter.C");

    // --- Get the parameters ----------------------------------------
    TString  type        = fOptions.Get("mom");
    Bool_t   satelitte   = fOptions.Has("satelitte");
    Double_t outlierFMD  = fOptions.AsDouble("outlier-fmd");
    Double_t outlierSPD  = fOptions.AsDouble("outlier-spd");
    TString  afterburner = fOptions.Get("afterburner");
    Int_t    abType      = fOptions.Get("ab-type");
    Int_t    abOrder     = fOptions.Get("ab-order");
    
    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwardFlow.C(\"%s\",%d,%d,%f,%f,\"%s\",%d,%d)",
		      type.Data(), mc, satelitte, outlierFMD, outlierSPD,
		      afterburner.Data(), abType, abOrder));
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
