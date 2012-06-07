/**
 * @file   MakeQATrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:55:02 2012
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_trains
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to do Quality assurance
 * 
 * @ingroup pwglf_forward_trains
 * @ingroup pwglf_forward_qa
 */
class MakeQATrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train 
   */
  MakeQATrain(const char* name  = "Forward QA") 
    : TrainSetup(name), 
      fUseCent(false)
  {
    SetType(kESD);
  }
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
    AliAnalysisManager::SetCommonFileName("forward_qa.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2", par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwardQA.C(%d,%d)", mc, fUseCent));
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

    gROOT->Macro("AddTaskCentrality.C");
    const char* cname = "CentralitySelection";
    AliCentralitySelectionTask* ctask = 
      dynamic_cast<AliCentralitySelectionTask*>(mgr->GetTask(cname));
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
  const char* ClassName() const { return "MakeQATrain"; }
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
