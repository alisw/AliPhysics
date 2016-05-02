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
    fOptions.Add("delta-cut",  "X",      "Cut on weighted distance",1.5);
    fOptions.Add("dphi-shift", "RADIANS","Bending shift",           0.0045);
    fOptions.Add("shifted-dphi-cut", "RADIANS", "Cut on dPhi-phiBent",  -1.);
    fOptions.Add("mc",                   "For MC data",             false);
    fOptions.Add("multsel",              "Enable MultSelection",    false);
    fOptions.Add("abs-min-cent","PERCENT","Absolute least cent.",   -1);
    fOptions.Add("reweight",   "FILE",    "File with weights",      "");
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
			      fOptions.AsString("reweight")));
    AliAnalysisTaskSE* task = reinterpret_cast<AliAnalysisTaskSE*>(ret);
    if (!task) return;
    
    // --- Figure out the trigger options ----------------------------
    TString trg = fOptions.Get("trig");
    trg.ToUpper();
    UInt_t  sel = AliVEvent::kINT7;
    if      (trg.EqualTo("MB"))    sel = AliVEvent::kMB;
    else if (trg.EqualTo("V0AND")) sel = AliVEvent::kINT7;
    else if (trg.EqualTo("V0OR"))  sel = AliVEvent::kCINT5;
    else if (trg.EqualTo("ANY"))   sel = AliVEvent::kAny;
    task->SelectCollisionCandidates(sel);

    // --- Set various options on task -------------------------------
    const char* defCent = DefaultCentBins();
    FromOption(task, "CentralityMethod","cent", 	    "V0M");
    FromOption(task, "CentralityAxis",  "cent-bins",        defCent);
    FromOption(task, "EtaAxis",         "eta-bins",         "r16:2");
    FromOption(task, "IPzAxis",         "ipz-bins",         "u15");
    FromOption(task, "DeltaCut",	"delta-cut",	    1.5);
    FromOption(task, "TailDelta",	"tail-delta",	    5.);
    FromOption(task, "MaxDelta",	"max-delta",	    25.);
    FromOption(task, "DPhiShift",	"dphi-shift",	    0.0045);
    FromOption(task, "ShiftedDPhiCut",	"shifted-dphi-cut",-1.);
    FromOption(task, "AbsMinCent",      "abs-min-cent",    -1.);

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
	
    task->Print("");    
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
