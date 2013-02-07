/**
 * @file   MakeFlowTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:51:48 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains
 * 
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ flow@f$
 * 
 *
 * @ingroup pwglf_forward_flow
 * @ingroup pwglf_forward_trains
 */
class MakeFlowTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train (free form)
   */
  MakeFlowTrain(const char* name)
  : TrainSetup(name)
  {
    fOptions.Add("mc", "Do MC analysis");
    fOptions.Add("mom", "Flow moments to analyse", "234", "234");
    fOptions.Add("eta-gap", "Whether to use an eta-gap", "[no,yes,both]", "both");
    fOptions.Add("eg-value", "Set value in |eta| of the eta gap", "2.0");
    fOptions.Add("use-cent", "Whether to use the impact parameter for centrality");
    fOptions.Add("afterburner", "What to afterburn", "[eta,phi,b,pid]", "");
    fOptions.Add("ab-type", "Type of afterburner", "1|2|3|4", "");
    fOptions.Add("ab-order", "Order of afterburner", "2|3|4|5|6", "");
    fOptions.Add("sat-vtx", "Whether to use satellite interactions");
    fOptions.Add("outlier-fmd", "Outlier cut for FMD", "NSIGMA", "4.0");
    fOptions.Add("outlier-spd", "Outlier cut for SPD", "NSIGMA", "0.0");
    fOptions.Set("type", "AOD");
    fOptions.Show(std::cout);
  }
protected:
  //__________________________________________________________________
  /** 
   * Create the tasks 
   * 
   * @param mgr Manager 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("AnalysisResults.root");

    // --- Load libraries/pars ---------------------------------------
    fHelper->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Get the parameters ----------------------------------------
    TString  type    = fOptions.Get("mom");
    TString  etaGap  = fOptions.Get("eta-gap");
    Double_t egValue = fOptions.AsDouble("eg-value");
    Bool_t   useCent = fOptions.AsBool("use-cent");
    TString  addFlow = fOptions.Get("afterburner");
    Int_t    abType  = fOptions.AsInt("ab-type");
    Int_t    abOrder = fOptions.AsInt("ab-order");
    Bool_t   satVtx  = fOptions.AsBool("sat-vtx");
    Double_t fmdCut  = fOptions.AsDouble("outlier-fmd");
    Double_t spdCut  = fOptions.AsDouble("outlier-spd");
    Bool_t   mc      = fOptions.AsBool("mc"); 

    // --- Add the task ----------------------------------------------
    if (etaGap.Contains("no") || etaGap.Contains("false") ||
        etaGap.Contains("both"))
          gROOT->Macro(Form("AddTaskForwardFlow.C(\"%s\",%d,%d,%f,%f,%f,%d,%d,\"%s\",%d,%d)",
		      type.Data(),
		      kFALSE,
		      mc, 
		      fmdCut, 
		      spdCut,
		      egValue,
		      useCent,
		      satVtx, 
		      addFlow.Data(), 
		      abType, 
		      abOrder));
    
    if (etaGap.Contains("yes") || etaGap.Contains("true") ||
        etaGap.Contains("both"))
          gROOT->Macro(Form("AddTaskForwardFlow.C(\"%s\",%d,%d,%f,%f,%f,%d,%d,\"%s\",%d,%d)",
		      type.Data(),
		      kTRUE,
		      mc, 
		      fmdCut, 
		      spdCut, 
		      egValue,
		      useCent,
		      satVtx, 
		      addFlow.Data(), 
		      abType, 
		      abOrder));
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
  /** 
   * Create MC input handler.  Since this train does not use the MC
   * information from @c galice.root, @c Kinematics.root, and @c
   * TrackRefs.root directly, we define this member function to return
   * null even when doing MC analysis.  This train _only_ looks at the
   * AOD object of the _processed_ MC data.
   * 
   * @return Always 0
   */
  AliVEventHandler* CreateMCHandler(UShort_t, bool)
  { 
    return 0;
  }
  //__________________________________________________________________
  /** 
   * Get the class name of the train setup 
   * 
   * @return Class name 
   */
  const char* ClassName() const { return "MakeFlowTrain"; }
  //__________________________________________________________________
};
//
// EOF
//
