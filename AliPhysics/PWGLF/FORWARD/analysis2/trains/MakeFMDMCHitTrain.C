/**
 * @file   MakeFMDMCHitTrain.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Jun  1 13:54:47 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains_specific
 */
#include "TrainSetup.C"
// #include "ParUtilities.C"

//====================================================================
/**
 * Analysis train to make Forward and Central MC corrections
 * 
 *
 * @ingroup pwglf_forward_mc
 * @ingroup pwglf_forward_trains_specific
 */
class MakeFMDMCHitTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train (free form)
   */
  MakeFMDMCHitTrain(const  char* name) 
    : TrainSetup(name)
  {
    fOptions.Add("use-tuple", "Whether to make an NTuple of hits");
    fOptions.Set("type", "ESD");
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
    AliAnalysisManager::SetCommonFileName("forward_mchits.root");


    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    fRailway->LoadLibrary("Proof");
    fRailway->LoadLibrary("Gui"); // Sigh! CDB depends on GUI!
    fRailway->LoadLibrary("CDB");
    fRailway->LoadLibrary("RAWDatabase");
    fRailway->LoadLibrary("STEER");
    fRailway->LoadLibrary("FMDbase");
    fRailway->LoadLibrary("FMDsim");
    fRailway->LoadLibrary("PWGLFforwardhit");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    if (!mgr->GetMCtruthEventHandler()) 
      Fatal("CreateTasks", "No MC truth handler");
    
    TString args = TString::Format("%d,%d", 
				   fOptions.AsBool("use-tuple"),
				   fOptions.AsInt("verbose"));
    if (!CoupleCar("AddTaskFMDMCHit.C", args))
      Fatal("CreateTasks", "Couldn't add our task");
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
      Fatal("CreatePhysicsSelection", "Couldn't get PhysicsSelection (%p)", ps);

    // --- Ignore trigger class when selecting events.  This means ---
    // --- that we get offline+(A,C,E) events too --------------------
    // ps->SetSkipTriggerClassSelection(true);
  }
  //__________________________________________________________________
  /** 
   * @return 0 - AOD disabled 
   */
  virtual AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  /** 
   * Do not the centrality selection
   */
  // void CreateCentralitySelection(Bool_t) {}
  //__________________________________________________________________
  const char* ClassName() const { return "MakeFMDMCHitTrain"; }
};

//
// EOF
//
