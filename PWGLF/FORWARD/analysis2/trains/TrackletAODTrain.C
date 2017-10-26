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
class TrackletAODTrain : public TrainSetup
{
public:
  /** 
   * Constructor. 
   * 
   * @param name     Name of train (free form)
   */
  TrackletAODTrain(const  TString& name) 
    : TrainSetup(name)
  {
    // Other options 
    fOptions.Add("copy", "LIST",    "',' separated list to copy","cent");
    fOptions.Add("mc-tracks", "Enable MC track filter", false);
    fOptions.Add("tpc-ep",    "Use TPC event plane");
    fOptions.Add("cent",      "Use centrality");
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
  const char* ClassName() const { return "TrackletAODTrain"; }
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
