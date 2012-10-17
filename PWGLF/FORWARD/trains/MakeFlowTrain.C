/**
 * @file   MakeFlowTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:51:48 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains_specific
 * 
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ flow@f$
 * 
 *
 * @ingroup pwglf_forward_flow
 * @ingroup pwglf_forward_trains_specific
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
  : TrainSetup(name),
    fType("123456"), 
    fOutlierCutFMD(4.1), 
    fOutlierCutSPD(4.1),
    fAddFlow(""),
      fAddFType(0),
    fAddFOrder(0)
  {
    fOptions.Add("mom", "MOMENTS", "Flow moments to determine", "123456");
    fOptions.Add("outlier-fmd", "NSIGMA", "Cut on outliers in FMD", "4.1");
    fOptions.Add("outlier-spd", "NSIGMA", "Cut on outliers in SPD", "4.1");
    fOptions.Add("afterburner", "WHAT", "What to afterburn", "eta phi b pid");
    fOptions.Add("ab-type", "1|2|3|4", "Type of afterburner", "");
    fOptions.Add("ab-order", "1|2|3|4|5|6", "Order used by afterburner", "");
    fOptions.Add("satelitte", "Whether to use satelitte interactions");
  }
protected:
  //__________________________________________________________________
  /** 
   * Create the tasks 
   * 
   * @param par  Whether to use par files 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("AnalysisResults.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Get the parameters ----------------------------------------
    TString  type        = fOptions.Get("mom");
    Bool_t   mc          =  mgr->GetMCtruthEventHandler() != 0;
    Bool_t   satelitte   = fOptions.Has("satelitte");
    Double_t outlierFMD  = fOptions.AsDouble("outlier-fmd");
    Double_t outlierSPD  = fOptions.AsDouble("outlier-spd");
    TString  afterburner = fOptions.Get("afterburner");
    Int_t    abType      = fOptions.Get("ab-type");
    Int_t    abOrder     = fOptions.Get("ab-order");
    
    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwardFlow.C(\"%s\",%d,%d,%f,%f,\"%s\",%d,%d)",
		      type.Data(), mc, satelitte, outlierFMD, outlierSPD,
		      afterburner.Data(), abType, abOrder));
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
  AliVEventHandler* CreateMCHandler(UShort_t, bool) { return 0; }
  //__________________________________________________________________
  /** 
   * Get the class name of the train setup 
   * 
   * @return Class name 
   */
  const char* ClassName() const { return "MakeFlowTrain"; }
};
//
// EOF
//
