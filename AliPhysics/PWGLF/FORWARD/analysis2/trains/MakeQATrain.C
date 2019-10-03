/**
 * @file   MakeQATrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:55:02 2012
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_trains_specific
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to do Quality assurance
 * 
 * @ingroup pwglf_forward_trains_specific
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
    : TrainSetup(name)
  {
    fOptions.Add("cent", "Use centrality");
    fOptions.Set("type", "ESD");
    fOptions.Add("corr", "DIR", "Corrections dir", "");
  }
protected:
  //__________________________________________________________________
  /** 
   * Create the tasks 
   * 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_qa.root");

    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;

    // --- Add the task ----------------------------------------------
    if (!CoupleSECar("AddTaskForwardQA.C", 
		     Form("%d,%d", mc, fOptions.Has("cent")),
		     AliVEvent::kAny))
      Fatal("CreateTasks", "Failed to add ForwardQA task");

    TString  cor = "";
    if (fOptions.Has("corr")) cor = fOptions.Get("corr"); 
    if (!cor.IsNull()) {
      fRailway->LoadAux(Form("%s/fmd_corrections.root",cor.Data()), true);
    }
  }
  /** 
   * Create entrality selection if enabled 
   * 
   * @param mc   Whether this is MC or not
   */
  virtual void CreateCentralitySelection(Bool_t mc)
  {
    if (!fOptions.Has("cent")) return;
    TrainSetup::CreateCentralitySelection(mc);    
  }
  /** 
   * Crete output handler - we don't want one here. 
   * 
   * @return 0
   */
  AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeQATrain"; }
};

//
// EOF
//
