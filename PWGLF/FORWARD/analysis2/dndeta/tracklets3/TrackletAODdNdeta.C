/**
 * @file   TrackletAODdNdeta.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:51:47 2016
 * 
 * @brief  A tracklet dNdeta train
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
#include "TrainSetup.C"
#ifndef __CINT__
#include <AliESDInputHandlerRP.h>
// #include <AliMixInputEventHandler.h>
// #include <AliMixEventPool.h>
// #include <AliMixEventCutObj.h>
#include <TMacro.h>
#else
class AliAnalysisTaskSE;
#endif

/** 
 * Train definition to make custom AOD with tracklets 
 * 
 * Run as 
 @verbatim 
 runTrain --class=TrackletAODTrain --name=NAME --url=URL 
 @endverbatim 
 *
 * See also 
 *
 * - https://twiki.cern.ch/twiki/bin/view/ALICE/FMDTrains
 * - http://hehi00.nbi.dk:8888/pwglfforward/train_setup_doc.html
 * 
 * This train uses a custom PAR file (RubensCode.par) of the classes
 *
 * - AliTrackletTaskMulti
 * - AliITSMultRecBg
 *
 * If these classes lived in a compiled AliPhysics library (Say
 * libPWGUD.so), the we wouldn't need that PAR file.
 *
 * @ingroup pwglf_forward_tracklets
 */
struct TrackletAODdNdeta : public TrainSetup
{
  /** 
   * Constructor.  This sets up the available options 
   * 
   * @param name Free form name 
   */
  TrackletAODdNdeta(const char* name)
    : TrainSetup(name)
  {
    // Define all our options here
    const char* defCent = DefaultCentBins();
    fOptions.Add("trig",       "NAME",   "Trigger to use",         "V0AND");
    fOptions.Add("cent",       "METHOD", "Centrality selector",    "V0M");
    fOptions.Add("cent-bins",  "BINS",   "Centrality bins",        defCent); 
    fOptions.Add("eta-bins",   "BINS",   "Eta bins",               "r16:-2:+2");
    fOptions.Add("ipz-bins",   "BINS",   "IPz bins",               "u15");
    fOptions.Add("max-delta",  "X",      "Cut on weighted distance",25.);
    fOptions.Add("tail-delta", "X",      "Tail cut on distance",    5.);
    fOptions.Add("tail-max",   "X",      "Tail cut on distance",    -1.);
    fOptions.Add("delta-cut",  "X",      "Cut on weighted distance",1.5);
    fOptions.Add("dphi-shift", "RADIANS","Bending shift",           0.0045);
    fOptions.Add("shifted-dphi-cut", "RADIANS", "Cut on dPhi-phiBent",  -1.);
    fOptions.Add("mc",                   "For MC data",             false);
    fOptions.Add("multsel",              "Enable MultSelection",    false);
    fOptions.Add("abs-min-cent","PERCENT","Absolute least cent.",   -1);
    fOptions.Add("max-ntracklet", "N",   "For FAKE centrality",    6000);
    fOptions.Add("reweigh",      "FILE", "File with weights",      "");
    fOptions.Add("reweigh-calc", "MODE", "prod,square,sum",        "prod");
    fOptions.Add("reweigh-mask", "MASK", "Tracklet mask for weighing", 0xFF);
    fOptions.Add("reweigh-veto", "VETO", "Tracklet veto for weighing", 0x0);
    fOptions.Add("reweigh-inv",          "W -> 1/W",                 false);
    fOptions.SetDescription("Analyse AOD for dN/deta from tracklets");
    fOptions.Set("type", "AOD");
  }
  /** 
   * The default centrality bins 
   * 
   * @return The list of centrality bins 
   */
  const char* DefaultCentBins() const
  {
    return "0-5-10-20-30-40-50-60-70-80-90";
  }
  /** 
   * Create output handler.  Overloaded here to set no output handler,
   * as this train does not define any AOD output.
   * 
   * @return Always null
   */
  virtual AliVEventHandler* CreateOutputHandler(UShort_t)
  {
    return 0;
  }
  /** 
   * Create centrality selection, and add to manager
   * 
   * @param mc Whether this is for MC 
   */
  virtual void CreateCentralitySelection(Bool_t mc)
  {
    mc               = fOptions.Has("mc");
    if (!mc) mc      = fRailway->IsMC();     
    if (fOptions.AsBool("multsel"))
      TrainSetup::CreateCentralitySelection(mc);
  }
  /** 
   * Create our tasks.  This uses the interpreter to make the object.
   * 
   * @param mgr 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
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
    if      (trg.EqualTo("MB"))    sel = AliVEvent::kMB;
    else if (trg.EqualTo("V0AND")) sel = AliVEvent::kINT7;
    else if (trg.EqualTo("V0OR"))  sel = AliVEvent::kCINT5;
    else if (trg.EqualTo("ANY"))   sel = AliVEvent::kAny;
    task->SelectCollisionCandidates(sel);

    // --- Figure out calculation mode -------------------------------
    TString calc = fOptions.Get("reweigh-calc"); calc.ToUpper();
    UChar_t mcal = 0;
    if      (calc.BeginsWith("PROD")) mcal = 0;
    else if (calc.BeginsWith("SQ"))   mcal = 1;
    else if (calc.BeginsWith("SUM"))  mcal = 2;
    else if (calc.BeginsWith("AV"))   mcal = 3;
    
    // --- Set various options on task -------------------------------
    const char* defCent = DefaultCentBins();
    FromOption(task, "CentralityMethod","cent", 	    "V0M");
    FromOption(task, "CentralityAxis",  "cent-bins",        defCent);
    FromOption(task, "AbsMinCent",      "abs-min-cent",    -1.);
    FromOption(task, "MaxNTracklet",    "max-ntracklet",    6000.);
    FromOption(task, "EtaAxis",         "eta-bins",         "r16:2");
    FromOption(task, "IPzAxis",         "ipz-bins",         "u15");
    FromOption(task, "DeltaCut",	"delta-cut",	    1.5);
    FromOption(task, "TailDelta",	"tail-delta",	    5.);
    FromOption(task, "TailMaximum",	"tail-max",	    -1);
    FromOption(task, "MaxDelta",	"max-delta",	    25.);
    FromOption(task, "DPhiShift",	"dphi-shift",	    0.0045);
    FromOption(task, "ShiftedDPhiCut",	"shifted-dphi-cut",-1.);
    FromOption(task, "WeightMask",      "reweigh-mask",     0xFF);
    FromOption(task, "WeightVeto",      "reweigh-veto",     0x0);
    SetOnTask (task, "WeightCalc",                          mcal);
    FromOption(task, "WeightInverse",   "reweigh-inv",      false);
    task->Print("");    
  }
  //__________________________________________________________________
  /** 
   * Overloaded to create new Collect.C 
   * 
   * @param asShellScript 
   */
  void SaveSetup(Bool_t asShellScript)
  {
    TrainSetup::SaveSetup(asShellScript);

    SaveCollect();
  }
  /** 
   * Make script to collect final
   * @f$\mathrm{d}N_{\mathrm{ch}}/\mathrm{d}\eta@f$ from this and
   * other pass.
   * 
   */
  void SaveCollect()
  {
    std::ofstream o("Collect.C");
    if (!o) { 
      Error("TrackletAODdNdeta::SavePost", "Failed to open Collect.C");
      return;
    }
    o << "// Created by " << ClassName() << "\n"
      << "// Will draw dN/deta results from produced files\n"
      << "// \n"
      << "// Arguments:\n"
      << "//   other     Directory of other result (MC or real)\n"
      << "//   output    Optional output directory name\n"
      << "//   proc      Bit mask of processing options\n"
      << "//             - 0x1  Unit normalisation of C\n"
      << "//             - 0x2  Constant normalisation of C\n"
      << "//             - 0x4  eta-dependent normalisation of C\n"
      << "//             - 0x8  (eta,IPz)-dependent normalisation of C\n"
      << "//   viz       Bit mask of visualisation options\n"
      << "//             - 0x0001  General information\n"
      << "//             - 0x0002  Parameters used\n"
      << "//             - 0x0004  Show weights if used\n"
      << "//             - 0x0008  Show dN/deta \n"
      << "//             - 0x0010  Show species\n"
      << "//             - 0x0020  Show delta distribtions\n"
      << "//             - 0x0040  Show details of the calculation\n"
      << "//             - 0x0080  Reserved\n"
      << "//             - 0x0100  Generate PDF and PNGs\n"
      << "//             - 0x0200  Pause after each page\n"
      << "//             - 0x0400  Generate plots in landscape orientation\n"
      << "//             - 0x0800  Use an alternate marker\n"
      << "//   n         Number of centrality bins to process\n" 
      << "// \n"
      << "// Output is generated in the sub-directory\n"
      << "// \n";
    if (fRailway->IsMC()) 
      o << "//    " << fEscapedName << "_<other>\n";
    else
      o << "//    <other>_" << fEscapedName << "\n";
    o << "//\n"
      << "// where <other> is the first argument given\n"
      << "void Collect(const char* other,\n"
      << "             const char* output=0\n"
      << "             UInt_t      proc=0x2,\n"
      << "             UInt_t      viz=0x32f,\n"
      << "             UInt_t      n=10)\n"
      << "{\n"
      << "  TString fwd=\n"
      << "    \"$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta/tracklets3\";\n"
      << "  gSystem->AddIncludePath(Form(\"-I%s\",fwd.Data()));\n"
      << "  gROOT->LoadMacro(Form(\"%s/AliTrackletAODUtils.C+g\",fwd.Data()));"
      << "\n"
      << "  gROOT->LoadMacro(Form(\"%s/AliTrackletdNdeta2.C+g\",fwd.Data()));"
      << "\n";
    TString oName(fEscapedName); oName.Append("/AnalysisResults.root");
    if (fRailway->IsMC()) 
      o << "  TString realFile = Form(\"%s/AnalysisResults.root\",other);\n"
	<< "  TString simFile  = \"" << oName << "\";\n"
	<< "  TString outFile  = Form(\""<<fEscapedName <<"_%s\",other);\n";
    else
      o << "  TString realFile = \"" << oName << "\";\n"
	<< "  TString simFile  = Form(\"%s/AnalysisResults.root\",other);\n"
	<< "  TString outFile  = Form(\"%s_"<<fEscapedName <<"\",other);\n";
    o << "  if (output && output[0] != '\0') outFile = output;\n"
      << "  AliTrackletdNdeta2* p = new AliTrackletdNdeta2;\n"
      << "  if (proc & 0x1) outFile.Append(\"_unit\");\n"
      << "  if (proc & 0x2) outFile.Append(\"_const\");\n"
      << "  if (proc & 0x4) outFile.Append(\"_eta\");\n"
      << "  if (proc & 0x8) outFile.Append(\"_etaipz\");\n"
      << "  \n"
      << "  p->Run(proc,viz,n,realFile,simFile,outFile);\n"
      << "  \n"
      << "  gROOT->LoadMacro(Form(\"%s/ExtractGSE2.C\",fwd.Data()));\n"
      << "  ExtractGSE2(outFile);\n"
      << "}\n"
      << "// EOF\n"
      << std::endl;
    o.close();
  }
  /** 
   * Get the train setup name 
   * 
   * @return The train setup name 
   */
  const char* ClassName() const { return "TrackletAODdNdeta"; }
};
//
// EOF
//
