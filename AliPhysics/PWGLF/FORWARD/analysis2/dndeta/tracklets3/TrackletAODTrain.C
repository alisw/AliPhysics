/**
 * @file   TrackletAODTrain.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:52:10 2016
 * 
 * @brief  Tracklet AOD train
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
struct TrackletAODTrain : public TrainSetup
{
  /** 
   * Constructor.  This sets up the available options 
   * 
   * @param name Free form name 
   */
  TrackletAODTrain(const char* name)
    : TrainSetup(name)
  {
    // Define all our options here
    fOptions.Add("max-delta",        "X","Cut on weighted distance",25.);
    fOptions.Add("scale-dtheta",         "Scale dTheta" ,           true);
    fOptions.Add("dphi-window",      "X","dPhi window",             0.06);
    fOptions.Add("dtheta-window",    "X","dTheta window",           0.025);
    fOptions.Add("dphi-shift",       "X","Bending shift",           0.0045);
    fOptions.Add("phi-overlap-cut",  "X","Phi overlap cut",         0.005);
    fOptions.Add("z-eta-overlap-cut","X","Z-Eta overlap cut",       0.05);
    fOptions.Add("copy",         "LIST","',' separated list to copy","cent,v0");
    fOptions.Add("filter-mode",  "MODE","Filter strange clusters",  0);
    fOptions.Add("filter-weight","FILE","File with filter weights", "");
    fOptions.SetDescription("Create branch in AOD with tracklet info");
    
  }
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
    Bool_t needRP = true;
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
   * Create our tasks.  This uses the interpreter to make the object.
   * 
   * @param mgr 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    Info("CreateTasks", "Loading code");
    fRailway->LoadSource("FixPaths.C");
    fRailway->LoadSource("AliAODSimpleHeader.C");
    fRailway->LoadSource("AliSimpleHeaderTask.C");    
    fRailway->LoadSource("AliAODTracklet.C");
    fRailway->LoadSource("AliTrackletWeights.C");
    fRailway->LoadSource("AliTrackletAODUtils.C");
    fRailway->LoadSource("AliTrackletAODTask.C");

    // --- Task to copy header information ---------------------------
    TString cpy = fOptions.Get("copy");
    Info("", "What to copy: %s", cpy.Data());
    CoupleSECar("AddTaskCopyHeader.C", Form("\"%s\"", cpy.Data()),
		AliVEvent::kAny);

    // --- Task to create simple header ------------------------------
    if (!gROOT->ProcessLine("AliSimpleHeaderTask::Create()"))
      return;
    
    // --- Create the task using interpreter -------------------------
    Long_t             ret  =
      gROOT->ProcessLine(Form("AliTrackletAODTask::Create(\"%s\")",
			      fOptions.AsString("filter-weight")));
    AliAnalysisTaskSE* task = reinterpret_cast<AliAnalysisTaskSE*>(ret);
    if (!task) return;
    
    // --- Set various options on task -------------------------------
    FromOption(task, "MaxDelta",	"max-delta",	     25.);
    FromOption(task, "ScaleDTheta",	"scale-dtheta",	     true);
    FromOption(task, "DPhiWindow",	"dphi-window",	     0.06);
    FromOption(task, "DThetaWindow",	"dtheta-window",     0.025);
    FromOption(task, "DPhiShift",	"dphi-shift",	     0.0045);
    FromOption(task, "PhiOverlapCut",	"phi-overlap-cut"  , 0.005);
    FromOption(task, "ZEtaOverlapCut",	"z-eta-overlap-cut", 0.05);
    FromOption(task, "FilterMode",      "filter-mode",       0);
    
    task->Print("");    
  }
  /** 
   * Get the train setup name 
   * 
   * @return The train setup name 
   */
  const char* ClassName() const { return "TrackletAODTrain"; }
};
//
// EOF
//
