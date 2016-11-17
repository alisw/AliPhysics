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
class BothAODTrain : public TrainSetup
{
public:
  /** 
   * Constructor. 
   * 
   * @param name     Name of train (free form)
   */
  BothAODTrain(const  TString& name) 
    : TrainSetup(name)
  {
    // Forward options 
    fOptions.Add("run",  "NUMBER",  "Run number for corrs", 0);
    fOptions.Add("sys",  "SYSTEM",  "1:pp, 2:PbPb, 3:pPb", 0);
    fOptions.Add("snn",  "ENERGY",  "Center of mass energy in GeV", 0);
    fOptions.Add("field","STRENGTH","L3 field strength in kG", 0);
    fOptions.Add("corr", "DIR",     "Corrections dir", "");
    fOptions.Add("dead", "FILE",    "Additional dead-map script", "");
    fOptions.Add("max-strips", "NUMBER", 
                 "Maximum consequtive strips (MC)", 2);
    fOptions.Add("forward-config", "FILE", "Forward configuration", 
		 "ForwardAODConfig.C");
    fOptions.Add("central-config", "FILE", "Central configuration", 
		 "CentralAODConfig.C");
    fOptions.Add("no-central","Do not at SPD cluster task",false);
    fOptions.Add("cent",      "Use centrality");
    fOptions.Add("satelitte", "Use satelitte interactions");
    fOptions.Add("secmap",    "Use secondary maps to correct", false);
    fOptions.Add("hit-threshold", "CUT", "Threshold for hits", 0.9);

    // Other options 
    fOptions.Add("copy", "LIST",    "',' separated list to copy","cent");
    fOptions.Add("mc-tracks", "Enable MC track filter", false);
    fOptions.Add("tpc-ep",    "Use TPC event plane");

    // SPD tracklet options 
    fOptions.Add("max-delta",        "X","Cut on weighted distance",25.);
    fOptions.Add("scale-dtheta",         "Scale dTheta" ,           true);
    fOptions.Add("dphi-window",      "X","dPhi window",             0.06);
    fOptions.Add("dtheta-window",    "X","dTheta window",           0.025);
    fOptions.Add("dphi-shift",       "X","Bending shift",           0.0045);
    fOptions.Add("phi-overlap-cut",  "X","Phi overlap cut",         0.005);
    fOptions.Add("z-eta-overlap-cut","X","Z-Eta overlap cut",       0.05);
    fOptions.Add("filter-mode",   "MODE","Filter strange clusters",  0);
    fOptions.Add("filter-weight", "FILE","File with filter weights", "");
    fOptions.Add("need-clusters",        "If set, insist on RecPoints",false);
    fOptions.Set("type", "ESD");
  }
protected:
  /** 
   * Create the input handler.  This is overwritten from the base
   * class to allow using AliESDInputHandlerRP for rec. points., and
   * AliMixInputEventHandler if requested.
   * 
   * @param type Type of analysis 
   * @param needRec Whether we need rec-points (clusters)
   * 
   * @return The input handler 
   */
  AliVEventHandler* CreateInputHandler(UShort_t type, Bool_t needRec=false)
  {
    Bool_t needRP = fOptions.AsBool("need-clusters");
    return TrainSetup::CreateInputHandler(type, needRP);
  }
  /** 
   * Create the MC input handler.  Overwritten here to allow setting
   * the pre-read mode.
   * 
   * @param type Input type 
   * @param mc   True for MC 
   * 
   * @return The MC input handler 
   */
  AliVEventHandler* CreateMCHandler(UShort_t type, bool mc)
  {
    AliMCEventHandler* ret =
      static_cast<AliMCEventHandler*>(TrainSetup::CreateMCHandler(type,mc));
    if (ret) ret->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    return ret;
  }
  
  /** 
   * Couple tasks to run before main body 
   * 
   * @param mgr 
   * @param mc 
   * @param mask 
   */
  void CouplePreCars(AliAnalysisManager* mgr, Bool_t mc, UInt_t mask)
  {
    // --- Add TPC eventplane task
    if (fOptions.Has("tpc-ep")) CoupleSECar("AddTaskEventplane.C","", mask);

    // --- Task to copy header information ---------------------------
    TString cpy = fOptions.Get("copy");
    Info("", "What to copy: %s", cpy.Data());
    CoupleSECar("AddTaskCopyHeader.C", Form("\"%s\"", cpy.Data()), mask);
  }
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
   * Couple the forward multiplicity task 
   * 
   * @param mgr 
   * @param mc  
   * @param mask 
   */
  void CoupleForwardCar(AliAnalysisManager* mgr, Bool_t mc, UInt_t mask)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward.root");

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
      Fatal("CoupleForwardCar", "Failed to add forward task");

    SetOnFwd(fwd, "GetCorrections().SetUseSecondaryMap", Form("%d",sec));
    //SetOnFwd(fwd, "GetDensityCalculator()->SetHitThreshold", Form("%f", thr));
    if (mc)
      SetOnMCFwd(fwd, "GetTrackDensity().SetMaxConsequtiveStrips",
		 Form("%d",mSt));
    fRailway->LoadAux(gSystem->Which(gROOT->GetMacroPath(), fwdConfig), true);
    if (!corr.IsNull()) 
      fRailway->LoadAux(Form("%s/fmd_corrections.root",corr.Data()), true);
    if (!dead.IsNull())
      fRailway->LoadAux(Form("%s",dead.Data()), true);
  }
  /** 
   * Couple the central multiplicity task  (clusters)
   * 
   * @param mgr 
   * @param mc  
   * @param mask 
   */
  void CoupleClusterCar(AliAnalysisManager* mgr, Bool_t mc, UInt_t mask)
  {
    Bool_t noCentral = fOptions.Has("no-central");
    if (noCentral) return;
    
    // --- Get options -----------------------------------------------
    ULong_t  run  = fOptions.AsInt("run", 0);
    UShort_t sys  = fOptions.AsInt("sys", 0);
    UShort_t sNN  = fOptions.AsInt("snn", 0);
    UShort_t fld  = fOptions.AsInt("field", 0);
    UShort_t mSt  = fOptions.AsInt("max-strips", 2);
    Bool_t   sec  = fOptions.Has("secmap");
    TString  corr = "";
    if (fOptions.Has("corr")) corr = fOptions.Get("corr"); 

    // --- Add the task ----------------------------------------------
    TString cenConfig = fOptions.Get("central-config");
    AliAnalysisTask* cen = CoupleSECar("AddTaskCentralMult.C",
				       Form("%d,%ld,%d,%d,%d,\"%s\",\"%s\"", 
					    mc, run, sys, sNN, fld, 
					    cenConfig.Data(), corr.Data()),
				       mask);
    if (!cen)
      Fatal("CoupleClusterCar", "Failed to add central (cluster) task");
    fRailway->LoadAux(gSystem->Which(gROOT->GetMacroPath(),cenConfig),true);
    if (!corr.IsNull())
      fRailway->LoadAux(Form("%s/spd_corrections.root",corr.Data()), true);
    gROOT->ProcessLine(Form("((AliCentralMultiplicityTask*)%p)"
			    "->SetUseSecondary(%d)",
			    cen, sec));
  }
  /** 
   * Couple central multiplicity task (tracklets)
   * 
   * @param mgr 
   * @param mc 
   * @param mask 
   */
  void CoupleTrackletCar(AliAnalysisManager* mgr, Bool_t mc, UInt_t mask)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("tracklet.root");
    
    // Load our sources 
    fRailway->LoadSource("FixPaths.C");
    fRailway->LoadSource("AliAODSimpleHeader.C");
    fRailway->LoadSource("AliSimpleHeaderTask.C");    
    fRailway->LoadSource("AliAODTracklet.C");
    fRailway->LoadSource("AliTrackletWeights.C");
    fRailway->LoadSource("AliTrackletAODUtils.C");
    fRailway->LoadSource("AliTrackletAODTask.C");

    // --- Task to create simple header ------------------------------
    if (!gROOT->ProcessLine("AliSimpleHeaderTask::Create()"))
      Fatal("CoupleTrackletTask", "Failed to make simple header task");

    // --- Create the task using interpreter -------------------------
    Long_t             ret  =
      gROOT->ProcessLine(Form("AliTrackletAODTask::Create(\"%s\")",
			      fOptions.AsString("filter-weight")));
    AliAnalysisTaskSE* task = reinterpret_cast<AliAnalysisTaskSE*>(ret);
    if (!task) 
      Fatal("CoupleTrackletCar", "Failed to make tracklet task");
    
    // --- Set various options on task -------------------------------
    FromOption(task, "MaxDelta",	"max-delta",	     25.);
    FromOption(task, "ScaleDTheta",	"scale-dtheta",	     true);
    FromOption(task, "DPhiWindow",	"dphi-window",	     0.06);
    FromOption(task, "DThetaWindow",	"dtheta-window",     0.025);
    FromOption(task, "DPhiShift",	"dphi-shift",	     0.0045);
    FromOption(task, "PhiOverlapCut",	"phi-overlap-cut"  , 0.005);
    FromOption(task, "ZEtaOverlapCut",	"z-eta-overlap-cut", 0.05);
    FromOption(task, "FilterMode",      "filter-mode",       0);
  }
  /** 
   * Couple tasks to run after main body 
   * 
   * @param mgr 
   * @param mc 
   * @param mask 
   */
  void CouplePostCars(AliAnalysisManager* mgr, Bool_t mc, UInt_t mask)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("AnalysisResults.root");
    // --- Add MC particle task --------------------------------------
    if (!mc || !fOptions.Has("mc-tracks")) return;
    CoupleSECar("AddTaskMCParticleFilter.C","", mask);
  }    
  /** 
   * Create the tasks 
   * 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {  
    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    TString mp(gROOT->GetMacroPath());
    mp.Append(":$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2");
    mp.Append(":$(ALICE_ROOT)/ANALYSIS/macros");
    mp.Append(":$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2/dndeta/tracklets3");
    gROOT->SetMacroPath(mp);

    // --- Check if this is MC ---------------------------------------
    Bool_t mc   = mgr->GetMCtruthEventHandler() != 0;    
    UInt_t mask = AliVEvent::kAny;

    // --- Couple our tasks ------------------------------------------
    CouplePreCars    (mgr, mc, mask);
    CoupleForwardCar (mgr, mc, mask);
    CoupleClusterCar (mgr, mc, mask);
    CoupleTrackletCar(mgr, mc, mask);
    CouplePostCars   (mgr, mc, mask);      
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
  const char* ClassName() const { return "BothAODTrain"; }
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
  /** 
   * Write script to do dN/deta analysis 
   * 
   * @param asShellScript 
   */
  void SavedNdeta(Bool_t asShellScript)
  {
    if (!fRailway) { 
      Warning("BothAODTrain::SaveSetup", 
	      "Cannot make dNdeta.C script without helper");
      return;
    }
    
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    Bool_t              mc  = mgr && (mgr->GetMCtruthEventHandler() != 0);
    OptionList          uopts(fRailway->Options());
    OptionList          opts(fOptions);
    TString             cls("BothdNdetaTrain");
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
    opts.Remove("copy");
    opts.Remove("scale-dtheta");
    opts.Remove("dphi-window");
    opts.Remove("dtheta-window");
    opts.Remove("phi-overlap-cut");
    opts.Remove("z-eta-overlap-cut");
    opts.Remove("filter-k0s");
    // Various
    opts.Add("cent-bins",    "BINS",      "Centrality bins",         "");
    opts.Add("ipz-bins",     "BINS",      "IPz bins",               "u15");
    opts.Add("satellite",    "Restrict analysis to satellite events",false);
    opts.Add("trig",         "TRIGGER",   "Trigger type",            "INEL");
    opts.Add("filter",       "FILTER",    "Filter type", "OUTLIER|PILEUP-BIN");
    opts.Add("mc",           "Also analyse MC truth", fRailway->IsMC());
    opts.Add("multsel",      "Enable MultSelection",                 false);
    opts.Add("abs-min-cent", "PERCENT",   "Absolute least cent.",    -1);
    // For forward 
    opts.Add("scheme",       "FLAGS",   "Normalization scheme","TRIGGER,EVENT");
    opts.Add("trigEff",      "EFFICIENCY","Trigger efficiency",      1.);
    opts.Add("trigEff0",     "EFFICIENCY","0-bin trigger efficiency",1.);
    opts.Add("truth-config", "FILE",      "MC-Truth configuration", "");
    opts.Add("mean-ipz",     "MU",        "Mean of IPz dist.",       0);
    opts.Add("var-ipz",      "SIGMA",     "Variance of IPz dist.",   -1);
    // For tracklets 
    opts.Add("eta-bins",     "BINS",      "Eta bins",              "r16:-2:+2");
    opts.Add("tail-delta",   "X",         "Tail cut on distance",    5.);
    opts.Add("tail-max",     "X",         "Tail cut on distance",    -1.);
    opts.Add("delta-cut",    "X",         "Cut on weighted distance",1.5);
    opts.Add("shifted-dphi-cut","RADIANS","Cut on dPhi-phiBent",     -1.);
    opts.Add("reweigh",      "FILE",      "File with weights",       "");
    opts.Add("reweigh-calc", "MODE",      "prod,square,sum",         "prod");
    opts.Add("reweigh-mask", "MASK",      "Tracklet mask for weighing", 0xFF);
    opts.Add("reweigh-veto", "VETO",      "Tracklet veto for weighing", 0x0);
    opts.Add("reweigh-inv",               "W -> 1/W",                 false);
    
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
      
    const char* defConfig=
      "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dNdetaConfig.C";
    opts.Set("url", outUrl.GetUrl());
    opts.Set("type", "AOD");
    opts.Set("forward-config",defConfig);
    opts.Set("central-config",defConfig);
    opts.Set("truth-config",defConfig);
    if (!fDatimeString.IsNull()) opts.Set("date", fDatimeString);

    if (sys != 1) {
      opts.Set("cent", "default");
      opts.Set("trig", "V0AND");
      opts.Set("scheme", "default");
      opts.Set("cent-bins", "default");
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
  /** 
   * Code to run from the @c post.sh shell script 
   * 
   * @param f Output stream 
   */
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
