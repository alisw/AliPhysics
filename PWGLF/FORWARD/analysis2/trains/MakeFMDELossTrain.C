/**
 * @file   MakeFMDELossTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:53:02 2012
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_trains
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to do energy loss fits
 * 
 * @ingroup pwglf_forward_trains
 * @ingroup pwglf_forward_eloss
 */
class MakeFMDELossTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train 
   */
  MakeFMDELossTrain(const char* name  = "FMD Energy Loss")
    : TrainSetup(name),
      fUseCent(false)
  {
    SetType(kESD);
  }
  /** 
   * Whether to process per centrality bin
   * 
   * @param use 
   */
  void SetUseCentrality(Bool_t use) { fUseCent = use; }
protected:
  //__________________________________________________________________
  /** 
   * Create the tasks 
   * 
   * @param par  Whether to use par files 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(EMode /*mode*/, Bool_t par, AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_eloss.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2", par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskFMDELoss.C(%d,%d)", mc, fUseCent));
  }
  /** 
   * Create entrality selection if enabled 
   * 
   * @param mc   Whether this is MC or not
   * @param mgr  Analysis manager 
   */
  virtual void CreateCentralitySelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    if (!fUseCent) return;

    const char* name = "CentralitySelection";
    gROOT->Macro("AddTaskCentrality.C");
    AliCentralitySelectionTask* ctask = 
      dynamic_cast<AliCentralitySelectionTask*>(mgr->GetTask(name));
    if (!ctask) return;
    // ctask->SetPass(fESDPass);
    if (mc) ctask->SetMCInput();
  }
  /** 
   * Crete output handler - we don't want one here. 
   * 
   * @return 0
   */
  AliVEventHandler* CreateOutputHandler(EType) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeFMDElossTrain"; }
  //__________________________________________________________________
  void MakeOptions(Runner& r) 
  {
    TrainSetup::MakeOptions(r);
    r.Add(new Option("cent",   "Use centrality"));
  }
  //__________________________________________________________________
  void SetOptions(Runner& r)
  {
    TrainSetup::SetOptions(r);
    Option*   cent	= r.FindOption("cent");
    if (cent) fUseCent  = cent->AsBool();
  }
  Bool_t fUseCent; // Whether to use centrality or not 
};

//
// EOF
//
