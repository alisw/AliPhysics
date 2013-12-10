/**
 * @file   MakeRespMatrTrain.C
 * @author Valentina Zaccolo
 * @date   Fri Jan  11 14:47:26 2013
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains_specific
 * 
 * @ingroup pwglf_forward_multdist
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ Response Matrices@f$
 * 
 *
 * @ingroup pwglf_forward_multdist
 * @ingroup pwglf_forward_trains_specific
 */
class MakeRespMatrTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  
   * 
   * @param name     Name of train (free form)
   */
  MakeRespMatrTrain(const char* name)
    : TrainSetup(name)
  {
    fOptions.Add("trig",  "TYPE",       "Trigger type", "V0AND");
    fOptions.Add("vzMin", "CENTIMETER", "Min Ip Z",     -4);
    fOptions.Add("vzMax", "CENTIMETER", "Max Ip Z",     +4);
  }
protected:
  /** 
   * Make the tasks 
   * 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_response.root");

    // --- Load libraries/pars ---------------------------------------
    fHelper->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Get parameters --------------------------------------------
    TString  trig    = fOptions.Get("trig");
    Double_t vzMin   = fOptions.AsDouble("vzmin", -4);
    Double_t vzMax   = fOptions.AsDouble("vzmax", +4);

    // --- Form arguments --------------------------------------------
    TString args;
    args.Form("\"%s\",%f,%f",
	      trig.Data(), vzMin, vzMax);
    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskCreateRespMatr.C(%s);", args.Data()));
  }
  //__________________________________________________________________
  /** 
   * Do not the centrality selection
   */
  void CreateCentralitySelection(Bool_t, AliAnalysisManager*) {}
  //__________________________________________________________________
  /** 
   * Crete output handler - we don't want one here. 
   * 
   * @return 0
   */
  AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  /** 
   * Never ever make an MC input handler! (we're working on AODs,
   * right?)
   */
  AliVEventHandler* CreateMCHandler(UShort_t type, bool mc) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeRespMatrTrain"; }
  //__________________________________________________________________
  
};
//
// EOF
//
