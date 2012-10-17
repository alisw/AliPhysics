/**
 * @file   MakedNdetaTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:51:26 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains_specific
 * 
 */

#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ dN/d\eta@f$
 * 
 *
 * @ingroup pwglf_forward_dndeta
 * @ingroup pwglf_forward_trains_specific
 */
class MakedNdetaTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  
   * 
   * @param name     Name of train (free form)
   */
  MakedNdetaTrain(const char* name)
  : TrainSetup(name)
  {
    fOptions.Add("trig",     "TYPE", "Trigger type", "INEL");
    fOptions.Add("vzMin",    "CENTIMETER", "Min Ip Z", "-10");
    fOptions.Add("vzMax",    "CENTIMETER", "Max Ip Z", "+10");
    fOptions.Add("scheme",   "SCHEME", "Normalization scheme", "");
    fOptions.Add("trigEff",  "EFFICENCY", "Trigger effeciency", "1");
    fOptions.Add("trigEff0", "EFFICENCY", "0-bin trigger effeciency", "1");
    fOptions.Add("cent",     "Use centrality");
    fOptions.Add("cut-edges", "Cut acceptance edges");
  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param par  Whether to use par files 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_dndeta.root");

    // --- Load libraries/pars ---------------------------------------
    fHelper->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Get parameters --------------------------------------------
    TString  trig   = fOptions.Get("trig");
    TString  scheme = fOptions.Get("scheme");
    Double_t vzMin  = fOptions.AsDouble("vzmin", -10);
    Double_t vzMax  = fOptions.AsDouble("vzmax", +10);
    Double_t effT   = fOptions.AsDouble("trigEff", 1);
    Double_t effT0  = fOptions.AsDouble("trigEff0", 1);
    Bool_t   cent   = fOptions.Has("cent");
    Bool_t   edges  = fOptions.Has("cut-edges");

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwarddNdeta.C(\"%s\",%f,%f,%d,\"%s\",%d,%g,%g)",
		      trig.Data(), vzMin, vzMax, cent, scheme.Data(),
		      edges, effT, effT0));

    gROOT->Macro(Form("AddTaskCentraldNdeta.C(\"%s\",%f,%f,%d,\"%s\",%d,%g,%g)",
		      trig.Data(), vzMin, vzMax, cent, scheme.Data(),
		      edges, effT, effT0));

    gROOT->Macro(Form("AddTaskMCTruthdNdeta.C(\"%s\",%f,%f,%d,\"%s\",%d,%g,%g)",
		      trig.Data(), vzMin, vzMax, cent, scheme.Data(),
		      edges, effT, effT0));
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
  const char* ClassName() const { return "MakedNdetaTrain"; }
};
//
// EOF
//
