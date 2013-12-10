/**
 * @file   MakeMultTrain.C
 * @author Valentina Zaccolo
 * @date   Wed Nov  21 12:47:26 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains_specific
 * @ingroup pwglf_forward_multdist
 * 
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ Multiplicity Distributions@f$
 * 
 *
 * @ingroup pwglf_forward_multdist
 * @ingroup pwglf_forward_trains_specific
 */
class MakeMultTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  
   * 
   * @param name     Name of train (free form)
   */
 MakeMultTrain(const char* name)
  : TrainSetup(name)
  {
    fOptions.Add("trig",    "TYPE",       "Trigger type",     "V0AND");
    fOptions.Add("vzMin",   "CENTIMETER", "Min Ip Z",         -4);
    fOptions.Add("vzMax",   "CENTIMETER", "Max Ip Z",         +4);
    fOptions.Add("lowCent", "%",          "Min Centrality",   0); 
    fOptions.Add("highCent","%",          "Max Centrality",   0);
    fOptions.Add("nBins",   "N",          "Max Multiplicity", 500);
  }
protected:
  /** 
   * Create the tasks 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_multiplicity.root");

    // --- Load libraries/pars ---------------------------------------
    fHelper->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Get parameters --------------------------------------------
    TString  trig     = fOptions.Get("trig");
    Double_t vzMin    = fOptions.AsDouble("vzmin", -4);
    Double_t vzMax    = fOptions.AsDouble("vzmax", +4);
    Int_t    lowCent  = fOptions.AsInt("lowCent",  0);
    Int_t    highCent = fOptions.AsInt("highCent", 0);
    Int_t    nBins    = fOptions.AsInt("nBins",    500);

    // --- Form arguments --------------------------------------------
    TString args;
    args.Form("\"%s\",%f,%f,%d,%d,%d",
	      trig.Data(), vzMin, vzMax, lowCent, highCent, nBins);
    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskMultDists.C(%s);", args.Data()));
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
  //__________________________________________________________________
  const char* ClassName() const { return "MakeMultTrain"; }
};
//
// EOF
//
