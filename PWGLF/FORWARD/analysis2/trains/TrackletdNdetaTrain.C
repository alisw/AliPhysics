/**
 * @file   TrackletdNdetaTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:51:26 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains_specific
 * 
 */

#include "TrainSetup.C"
#include <TF1.h>

//====================================================================
/**
 * Analysis train to make @f$ dN/d\eta@f$
 * 
 *
 * @ingroup pwglf_forward_dndeta
 * @ingroup pwglf_forward_trains_specific
 */
class TrackletdNdetaTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  
   * 
   * @param name     Name of train (free form)
   */
  TrackletdNdetaTrain(const char* name)
    : TrainSetup(name)
  {
    fOptions.Add("trig",          "TYPE",      "Trigger type",       "INEL");
    fOptions.Add("trigEff",       "EFFICIENCY","Trigger efficiency", 1.);
    fOptions.Add("filter",        "TYPE",   "Filter type","OUTLIER|PILEUP-BIN");
    fOptions.Add("cent",          "ESTIMATOR", "Use centrality",     "none");
    fOptions.Add("cent-bins",     "BINS",      "Centrality bins",    "default");
    fOptions.Add("abs-min-cent",  "PERCENT", "Absolute least centrality", -1.);
    fOptions.Add("ipz-bins",      "BINS",      "IPz bins",           "u15");
    fOptions.Add("eta-bins",      "BINS",      "Eta bins",         "r16:-2:+2");
    // For simulations 
    fOptions.Add("mc",            "Also make dN/deta for MC truth",  false);
    // For SPD tracklets 
    fOptions.Add("max-delta",  "X",      "Cut on weighted distance",25.);
    fOptions.Add("tail-delta", "X",      "Tail cut on distance",    5.);
    fOptions.Add("tail-max",   "X",      "Tail cut on distance",    -1.);
    fOptions.Add("delta-cut",  "X",      "Cut on weighted distance",1.5);
    fOptions.Add("dphi-shift", "RADIANS","Bending shift",           0.0045);
    fOptions.Add("shifted-dphi-cut", "RADIANS", "Cut on dPhi-phiBent",  -1.);
    fOptions.Add("reweigh",    "FILE",   "File with weights",      "");
    fOptions.Add("reweigh-calc", "MODE", "prod,square,sum",        "prod");
    fOptions.Add("reweigh-mask", "MASK", "Tracklet mask for weighing", 0xFF);
    fOptions.Add("reweigh-veto", "VETO", "Tracklet veto for weighing", 0x0);
    fOptions.Add("reweigh-inv",          "W -> 1/W",                 false);
    
  }
protected:
  /** 
   * Create our tasks.  This uses the interpreter to make the object.
   * 
   * @param mgr 
   */
  void CoupleTrackletCar(AliAnalysisManager* mgr)
  {
    AliAnalysisManager::SetCommonFileName("tracklet_dndeta.root");
    TString fwd(gSystem->Getenv("ANA_SRC"));
    if (fwd.IsNull()) fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
    gROOT->SetMacroPath(Form("%s:%s/dndeta/tracklets3",
			     gROOT->GetMacroPath(), fwd.Data()));
    gSystem->AddIncludePath(Form("-I%s/dndeta/tracklets3", fwd.Data()));
    
    Info("CreateTasks", "Loading code");
    fRailway->LoadSource("FixPaths.C");
    fRailway->LoadSource("AliAODSimpleHeader.C");
    fRailway->LoadSource("AliAODTracklet.C");
    fRailway->LoadSource("AliTrackletWeights.C");
    fRailway->LoadSource("AliTrackletAODUtils.C");
    fRailway->LoadSource("AliTrackletAODdNdeta.C");

    // --- Create the task using interpreter -------------------------
    Bool_t   mc = fOptions.Has("mc");
    if (!mc) mc = fRailway->IsMC();     
    Long_t ret  =
      gROOT->ProcessLine(Form("AliTrackletAODdNdeta::Create(%d,\"%s\")",mc,
			      fOptions.AsString("reweigh")));
    AliAnalysisTaskSE* task = reinterpret_cast<AliAnalysisTaskSE*>(ret);
    if (!task) return;
    
    // --- Figure out the trigger options ----------------------------
    TString trg = fOptions.Get("trig"); trg.ToUpper();
    UInt_t  sel = AliVEvent::kINT7;
    UInt_t  mb  = AliVEvent::kINT1; // AliVEvent::kMB;
    if      (trg.EqualTo("MB"))      sel = mb;
    else if (trg.EqualTo("MBOR"))    sel = mb;
    else if (trg.EqualTo("INEL"))    sel = mb;
    else if (trg.EqualTo("INELGT0")) sel = mb;
    else if (trg.EqualTo("V0AND"))   sel = AliVEvent::kINT7;
    else if (trg.EqualTo("NSD"))     sel = AliVEvent::kINT7;
    else if (trg.EqualTo("V0OR"))    sel = AliVEvent::kCINT5;
    else if (trg.EqualTo("ANY"))     sel = AliVEvent::kAny;
    else if (trg.EqualTo("NONE"))    sel = 0;
    task->SelectCollisionCandidates(sel);
    Int_t minTrk = trg.EqualTo("INELGT0") ? 1 : 0;

    // --- Figure out calculation mode -------------------------------
    TString calc = fOptions.Get("reweigh-calc"); calc.ToUpper();
    UChar_t mcal = 0;
    if      (calc.BeginsWith("PROD")) mcal = 0;
    else if (calc.BeginsWith("SQ"))   mcal = 1;
    else if (calc.BeginsWith("SUM"))  mcal = 2;
    else if (calc.BeginsWith("AV"))   mcal = 3;
    
    // --- Set various options on task -------------------------------
    const char* defCent = "0-5-10-20-30-40-50-60-70-80-90";
    FromOption(task, "CentralityMethod", "cent", 	    "V0M");
    FromOption(task, "CentralityAxis",   "cent-bins",        defCent);
    FromOption(task, "EtaAxis",          "eta-bins",         "r16:2");
    FromOption(task, "IPzAxis",          "ipz-bins",         "u15");
    FromOption(task, "DeltaCut",	 "delta-cut",	     1.5);
    FromOption(task, "TailDelta",	 "tail-delta",	     5.);
    FromOption(task, "TailMaximum",	 "tail-max",	     -1);
    FromOption(task, "MaxDelta",	 "max-delta",	     25.);
    FromOption(task, "DPhiShift",	 "dphi-shift",	     0.0045);
    FromOption(task, "ShiftedDPhiCut",	 "shifted-dphi-cut", -1.);
    FromOption(task, "AbsMinCent",       "abs-min-cent",     -1.);
    FromOption(task, "WeightMask",       "reweigh-mask",     0xFF);
    FromOption(task, "WeightVeto",       "reweigh-veto",     0x0);
    SetOnTask (task, "WeightCalc",                           mcal);
    FromOption(task, "WeightInverse",    "reweigh-inv",      false);
    FromOption(task, "TriggerEfficiency","trigEff",          0.);
    SetOnTask (task, "MinEta1",                              minTrk);
    // if (mc && we) {
    //   TUrl wurl(fOptions.AsString("reweight"));
    //   TFile* wfile = TFile::Open(wurl.GetFile());
    //   if (!wfile) {
    // 	Warning("CreateTasks", "Failed to open weights file: %s",
    // 		wurl.GetUrl());
    // 	return;
    //   }
    //   TString wnam(wurl.GetAnchor());
    //   if (wnam.IsNull()) wnam = "weights";
    //   TObject* wobj = wfile->Get(wnam);
    //   if (!wobj) {
    // 	Warning("CreateTasks", "Failed to get weights %s from file %s",
    // 		wnam.Data(), wfile->GetName());
    // 	return;
    //   }
    //   if (!wobj->IsA()->InheritsFrom("AliTrackletWeights")) {
    // 	Warning("CreateTasks", "Object %s from file %s not an "
    // 		"AliTrackletWeights but a %s",
    // 		wnam.Data(), wfile->GetName(), wobj->ClassName());
    // 	return;
    //   }
    //   SetOnTaskGeneric(task, "Weights",
    // 		       Form("((AliTrackletWeights*)%p)", wobj));
    // }
    Printf("Print the generated task");
    task->Print("");    
  }
  /** 
   * Create the tasks 
   * 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    // fRailway->LoadLibrary("AOD");
    // gSystem->ListLibraries();
    // gSystem->Load("libAOD");

    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Add the task ----------------------------------------------
    CoupleTrackletCar(mgr);

  }
  //__________________________________________________________________
  /** 
   * Do not the centrality selection
   */
  //__________________________________________________________________
  // void CreateCentralitySelection(Bool_t) {}
  /** 
   * Do not create MC input handler 
   * 
   * @return Always null
   */
  AliVEventHandler* CreateMCHandler(UShort_t, bool) { return 0; }
  //__________________________________________________________________
  /** 
   * Crete output handler - we don't want one here. 
   * 
   * @return 0
   */
  AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "TrackletdNdetaTrain"; }
  //__________________________________________________________________
  /** 
   * Overloaded to create new draw.C 
   * 
   * @param asShellScript 
   */
  void SaveSetup(Bool_t asShellScript)
  {
    TrainSetup::SaveSetup(asShellScript);

    SaveCollect();
  }
  void SaveCollect()
  {
    std::ofstream o("Collect.C");
    if (!o) { 
      Error("TrackletdNdetaTrain::SavePost", "Failed to open Collect.C");
      return;
    }
    o << "// Created by " << ClassName() << "\n"
      << "// Will draw dN/deta results from produced files\n"
      << "// \n"
      << "void Collect(const char* other,\n"
      << "             const char* output,\n"
      << "             UInt_t      proc=0x22,\n"
      << "             UInt_t      viz=0x32f,\n"
      << "             UInt_t      n=10)\n"
      << "{\n"
      << "  TString fwd=\n"
      << "    \"$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta/tracklets3\";\n"
      << "  if (gSystem->Getenv(\"ANA_SRC\"))\n"
      << "    fwd = \"$ANA_SRC/dndeta/tracklets3\";\n"
      << "  gROOT->LoadMacro(Form(\"%s/Post.C\",fwd.Data()));\n"
      << "  const char* thisDir = \"" << fEscapedName << "\";\n";
    if (fRailway->IsMC())
      o << "  Post(thisDir,other,output,proc,viz,n);\n";
    else
      o << "  Post(other,thisDir,output,proc,viz,n);\n";
    o << "}\n"
      << "// EOF\n"
      << std::endl;
    o.close();
  }
  void PostShellCode(std::ostream& f)
  {
    f << "  echo \"=== Summarizing results ...\"\n"
      << "  aliroot -l -b -q ${prefix}Summarize.C\n"
      << "  echo \"=== Draw results ...\"\n"
      << "  aliroot -l -b -q ${prefix}Draw.C\n"
      << std::endl;
  }
};
//
// EOF
//
