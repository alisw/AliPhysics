/**
 * @file   MakeAODTrain.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:05:30 2011
 * 
 * @brief  Run first pass analysis - make AOD tree
 * 
 * @ingroup pwglf_forward_trains_specific
 */
#include "TrainSetup.C"
#include <sstream>

//====================================================================
/**
 * Analysis train to make Forward and Central multiplicity
 * 
 *
 * @ingroup pwglf_forward_aod
 * @ingroup pwglf_forward_trains_specific
 */
class MakeAODTrain : public TrainSetup
{
public:
  /** 
   * Constructor. 
   * 
   * @param name     Name of train (free form)
   */
  MakeAODTrain(const  TString& name) 
    : TrainSetup(name)
  {
    fOptions.Add("run",  "NUMBER",  "Run number for corrs", 0);
    fOptions.Add("sys",  "SYSTEM",  "1:pp, 2:PbPb, 3:pPb", 0);
    fOptions.Add("snn",  "ENERGY",  "Center of mass energy in GeV", 0);
    fOptions.Add("field","STRENGTH","L3 field strength in kG", 0);
    fOptions.Add("corr", "DIR",     "Corrections dir", "");
    fOptions.Add("dead", "FILE",    "Additional dead-map script", "");
    fOptions.Add("copy", "LIST",    "',' separated list to copy","cent");
    fOptions.Add("max-strips", "NUMBER", 
                 "Maximum consequtive strips (MC)", 2);
    fOptions.Add("forward-config", "FILE", "Forward configuration", 
		 "ForwardAODConfig.C");
    fOptions.Add("central-config", "FILE", "Central configuration", 
		 "CentralAODConfig.C");
    fOptions.Add("no-central","Do not at SPD cluster task",false);
    fOptions.Add("cent",      "Use centrality");
    fOptions.Add("tpc-ep",    "Use TPC event plane");
    fOptions.Add("satelitte", "Use satelitte interactions");
    fOptions.Add("secmap",    "Use secondary maps to correct", false);
    fOptions.Add("mc-tracks", "Enable MC track filter", false);
    fOptions.Add("hit-threshold", "CUT", "Threshold for hits", 0.9);
    fOptions.Set("type", "ESD");
  }
protected:
  /** 
   * Set a parameter on the forward task 
   * 
   * @param task Task to set on
   * @param call The call 
   * @param arg  The argument(s) for the call
   */
  void SetOnFwd(AliAnalysisTask* task,
		const TString& call,
		const TString& arg)
  {
    TString cmd;
    cmd.Form("((AliForwardMultiplicityBase*)%p)->%s(%s)",
	     task, call.Data(), arg.Data());
    gROOT->ProcessLine(cmd);
  }
  /** 
   * Set a parameter on the forward task 
   * 
   * @param task Task to set on
   * @param call The call 
   * @param arg  The argument(s) for the call
   */
  void SetOnMCFwd(AliAnalysisTask* task,
		  const TString& call,
		  const TString& arg)
  {
    TString cmd;
    cmd.Form("((AliForwardMCMultiplicityTask*)%p)->%s(%s)",
	     task, call.Data(), arg.Data());
    gROOT->ProcessLine(cmd);
  }
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
    fRailway->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;
    
    UInt_t mask = AliVEvent::kAny;
    // --- Add TPC eventplane task
    if (fOptions.Has("tpc-ep")) CoupleSECar("AddTaskEventplane.C","", mask);

    // --- Task to copy header information ---------------------------
    TString cpy = fOptions.Get("copy");
    Info("", "What to copy: %s", cpy.Data());
    CoupleSECar("AddTaskCopyHeader.C", Form("\"%s\"", cpy.Data()), mask);

    // --- Get options -----------------------------------------------
    ULong_t  run  = fOptions.AsInt("run", 0);
    UShort_t sys  = fOptions.AsInt("sys", 0);
    UShort_t sNN  = fOptions.AsInt("snn", 0);
    UShort_t fld  = fOptions.AsInt("field", 0);
    UShort_t mSt  = fOptions.AsInt("max-strips", 2);
    Double_t thr  = fOptions.AsDouble("hit-threshold", .9);
    Bool_t   sec  = fOptions.Has("secmap");
    TString  corr = "";
    TString  dead = "";
    if (fOptions.Has("corr")) corr = fOptions.Get("corr"); 
    if (fOptions.Has("dead")) dead = fOptions.Get("dead"); 
    
    // --- Add the task ----------------------------------------------
    TString fwdConfig = fOptions.Get("forward-config");
    AliAnalysisTask* fwd =
      CoupleSECar("AddTaskForwardMult.C",
		  Form("%d,%ld,%d,%d,%d,\"%s\",\"%s\",\"%s\"", 
		       mc, run, sys, sNN, fld, 
		       fwdConfig.Data(), corr.Data(),
		       dead.Data()), mask);
    if (!fwd)
      Fatal("CoupleCars", "Failed to add forward task");
    SetOnFwd(fwd, "GetCorrections().SetUseSecondaryMap", Form("%d",sec));
    SetOnFwd(fwd, "GetDensityCalculator()->SetHitThreshold", Form("%f", thr));
    if (mc)
      SetOnMCFwd(fwd, "GetTrackDensity().SetMaxConsequtiveStrips",
		 Form("%d",mSt));
    fRailway->LoadAux(gSystem->Which(gROOT->GetMacroPath(), fwdConfig), true);
    if (!corr.IsNull()) {
      // fRailway->LoadAux(Form("%s/fmd_corrections.root",corr.Data()), true);
      const char* names[] = {
	"secondary",	
	"elossfits",	
	// "vertexbias",	
	// "merging",	
	// "doublehit",	
	"acceptance",	
	"noisegain",
	0 };
      const char** ptr = names;
      while (*ptr) {
	fRailway->LoadAux(Form("%s/%s",corr.Data(),*ptr),true);
	ptr++;
      }
    }
    if (!dead.IsNull())
      fRailway->LoadAux(Form("%s",dead.Data()), true);

    // --- Add the task ----------------------------------------------
    Bool_t noCentral = fOptions.Has("no-central");
    AliAnalysisTask* cen = 0;
    TString cenConfig = "";
    if (!noCentral) { 
      cenConfig = fOptions.Get("central-config");
      cen = CoupleSECar("AddTaskCentralMult.C",
			Form("%d,%ld,%d,%d,%d,\"%s\",\"%s\"", 
			     mc, run, sys, sNN, fld, 
			     cenConfig.Data(), corr.Data()), mask);
      fRailway->LoadAux(gSystem->Which(gROOT->GetMacroPath(),cenConfig),true);
      if (!corr.IsNull())
	fRailway->LoadAux(Form("%s/spd_corrections.root",corr.Data()), true);
      gROOT->ProcessLine(Form("((AliCentralMultiplicityTask*)%p)"
			      "->SetUseSecondary(%d)",
			      cen, sec));
    }

    // --- Add MC particle task --------------------------------------
    if (mc && fOptions.Has("mc-tracks")) 
      CoupleSECar("AddTaskMCParticleFilter.C","", mask);
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
  //__________________________________________________________________
  const char* ClassName() const { return "MakeAODTrain"; }
  //__________________________________________________________________
  /** 
   * Overloaded to create new dNdeta.C and dndeta.sh in the output 
   * directory
   * 
   * @param asShellScript 
   */
  void SaveSetup(Bool_t asShellScript)
  {
    TrainSetup::SaveSetup(asShellScript);
    SaveSummarize();
    SavedNdeta(asShellScript);
  }
  void SavedNdeta(Bool_t asShellScript)
  {
    if (!fRailway) { 
      Warning("MakeAODTrain::SaveSetup", 
	      "Cannot make dNdeta.C script without helper");
      return;
    }
    
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    Bool_t              mc  = mgr && (mgr->GetMCtruthEventHandler() != 0);
    OptionList          uopts(fRailway->Options());
    OptionList          opts(fOptions);
    TString             cls("MakedNdetaTrain");
    TString             name(fName);
    Int_t               sys = fOptions.AsInt("sys", 0);
    if (name.Contains("aod")) name.ReplaceAll("aod", "dndeta");
    else                      name.Append("_dndeta");
    opts.Remove("run");
    opts.Remove("sys");
    opts.Remove("snn");
    opts.Remove("field");
    opts.Remove("tpc-ep");
    opts.Remove("corr");
    opts.Remove("max-strips");
    opts.Remove("mc-tracks");
    opts.Remove("secmap");
    opts.Remove("dead");
    opts.Add("abs-min-cent","PERCENT",   "Absolute least centrality",-1.);
    opts.Add("cent-bins",   "BINS",      "Centrality bins",          "");
    opts.Add("satellite",   "Restrict analysis to satellite events", false);
    opts.Add("trig",        "TRIGGER",   "Trigger type",             "INEL");
    opts.Add("filter",      "FILTER",    "Filter type", "OUTLIER|PILEUP-BIN");
    opts.Add("ipz-bins",    "BINS",      "Lower bound on Ip Z",      "u10");
    opts.Add("scheme",      "FLAGS",     "Normalization scheme", "TRIGGER,EVENT");
    opts.Add("trigEff",     "EFFICIENCY","Trigger efficiency",       1.);
    opts.Add("trigEff0",    "EFFICIENCY","0-bin trigger efficiency", 1.);
    opts.Add("mc",          "Also analyse MC truth",           fRailway->IsMC());
    opts.Add("truth-config","FILE",      "MC-Truth configuration", "");
    opts.Add("mean-ipz",    "MU",        "Mean of IPz dist.",        0);
    opts.Add("var-ipz",     "SIGMA",     "Variance of IPz dist.",    -1);
    
    // Rewrite our URL 
    TString outString = fRailway->OutputLocation();
    if (outString.IsNull()) outString = fEscapedName;
    TUrl    outUrl(outString);
    
    if (uopts.Find("pattern")) // && outString.EndsWith("AliAOD.root")) 
      uopts.Set("pattern", "*/AliAOD.root");
    // if (uopts.Find("concat")) uopts.Set("concat", true);
    if (uopts.Find("par")) uopts.Set("par", "task");

    std::stringstream s;
    uopts.Store(s, "", "&", false, true);
    outUrl.SetOptions(s.str().c_str());
      
    const char* defConfig="$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dNdetaConfig.C";
    opts.Set("url", outUrl.GetUrl());
    opts.Set("type",          "AOD");
    opts.Set("forward-config",defConfig);
    opts.Set("central-config",defConfig);
    opts.Set("truth-config",defConfig);
    if (!fDatimeString.IsNull()) opts.Set("date", fDatimeString);

    if (sys != 1) {
      opts.Set("cent",     "default");
      opts.Set("trig",     "V0AND");
      opts.Set("scheme",   "default");
      opts.Set("centBins", "default");
      SaveSetupROOT("dNdeta", cls, name, opts, &uopts);
      if (asShellScript) 
	SaveSetupShell("dndeta", cls, name, opts, &uopts);
    }
    else {
      name.ReplaceAll("dndeta", "dndeta_inel");
      SaveSetupROOT("dNdetaINEL", cls, name, opts, &uopts);
      if (asShellScript) 
	SaveSetupShell("dndeta_inel", cls, name, opts, &uopts);
      
      name.ReplaceAll("inel", "nsd");
      opts.Set("trig", "V0AND");
      SaveSetupROOT("dNdetaNSD", cls, name, opts, &uopts);
      if (asShellScript) 
	SaveSetupShell("dndeta_nsd", cls, name, opts, &uopts);
      
      name.ReplaceAll("nsd", "inelgt0");
      opts.Set("trig", "INELGT0");
      SaveSetupROOT("dNdetaINELGt0", cls, name, opts, &uopts);
      if (asShellScript) 
	SaveSetupShell("dndeta_inelgt0", cls, name, opts, &uopts);
    }
  }
  /** 
   * Write a ROOT script to draw summary 
   * 
   */
  void SaveSummarize()
  {
    std::ofstream f("Summarize.C");
    if (!f) { 
      Error("SaveSummarize", "Failed to open Summarize.C script");
      return;
    }
    f << "// Generated by " << ClassName() << "\n"
      << "// WHAT is a bit mask of\n"
      << "//   0x001     Event inspector\n"
      << "//   0x002     Sharing filter\n"
      << "//   0x004     Density calculator\n"
      << "//   0x008     Corrector\n"
      << "//   0x010     Histogram collector\n"
      << "//   0x020     Analysis step cartoon\n"
      << "//   0x040     Results\n"
      << "//   0x080     Central\n"
      << "//   0x100     Landscape\n"
      << "//   0x200     Pause\n"
      << "//\n"
      << "void Summarize(const char* filename=\"forward.root\",\n"
      << "               UShort_t what=0x1FF)\n"
      << "{\n"
      << "  const char* fwd=\"$ALICE_PHYSICS/PWGLF/FORWARD/analysis2\";\n"
      << "  gROOT->LoadMacro(Form(\"%s/DrawAODSummary.C\",fwd));\n"
      << "  DrawAODSummary(filename,what);\n"
      << "}\n"
      << "// EOF" << std::endl;
    f.close();
  }
  void PostShellCode(std::ostream& f)
  {
    f << "  echo \"=== Summarizing results ...\"\n"
      << "  aliroot -l -b -q ${prefix}Summarize.C\n"
      << std::endl;
  }
};
//
// EOF
//
