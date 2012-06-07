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
  : TrainSetup(name),
    fType("123456"), 
    fOutlierCutFMD(4.1), 
    fOutlierCutSPD(4.1),
    fAddFlow(""),
      fAddFType(0),
    fAddFOrder(0)
  {
    TrainSetup::SetType(kAOD);
  }
  /** 
   * Set the flow moments to calculate.  A string of or more of the
   * numbers from 1 to 6
   * 
   * @param m 
   */
  void SetMoments(const char* m) { fType = m; }
  /** 
   * Set parameters on the MC flow after burner 
   * 
   * @param f  How the flow should be parameterized 
   * @param t  Type of parameterization 
   * @param n  Order of flow to add 
   */
  void SetFlowAfterburner(const char* f, Int_t t, Int_t n) 
  { 
    fAddFlow   = f;
    fAddFType  = t;
    fAddFOrder = n;
  }
  /** 
   * Set wether to use displaced (or satellite) interaction vertices
   * 
   * @param use If true, use displaced interaction vertices 
   */
  void SetUseDispVtx(Bool_t use = kTRUE) { fDispVtx = use; }
  /** 
   * Set the outlier cuts 
   * 
   * @param fmd Cut for FMD 
   * @param spd Cut for SPD 
   */
  void SetOutlierCuts(Double_t fmd, Double_t spd) 
  { 
    fOutlierCutFMD = fmd; 
    fOutlierCutSPD = spd;
  }
protected:
  //__________________________________________________________________
  /** 
   * Create the tasks 
   * 
   * @param par  Whether to use par files 
   */
  void CreateTasks(EMode /*mode*/, Bool_t par, AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("AnalysisResults.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2", par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwardFlow.C(\"%s\",%d,%d,%f,%f,\"%s\",%d,%d)",
		      fType.Data(), 
		      fMC, 
		      fDispVtx, 
		      fOutlierCutFMD, 
		      fOutlierCutSPD, 
		      fAddFlow.Data(), 
		      fAddFType, 
		      fAddFOrder));
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
  AliVEventHandler* CreateOutputHandler(EType) { return 0; }
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
  AliVEventHandler* CreateMCHandler(EType /*type*/, bool /*mc*/)
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
  /** 
   * Declare options 
   * 
   * @param r Runner object to add options to 
   */
  void MakeOptions(Runner& r)
  {
    TrainSetup::MakeOptions(r);
    r.Add(new Option("mom", "Flow moments to analyse", "123456"));
    r.Add(new Option("afterburner", "What to afterburn", "[eta,phi,b,pid]"));
    r.Add(new Option("ab-type", "Type of afterburner", "1|2|3|4"));
    r.Add(new Option("ab-order", "Order of afterburner", "1|2|3|4|5|6"));
    r.Add(new Option("disp-vtx", "Whether to use satellite interactions"));
    r.Add(new Option("outlier-fmd", "Outlier cut for FMD", "NSIGMA"));
    r.Add(new Option("outlier-spd", "Outlier cut for SPD", "NSIGMA"));
  }
  //__________________________________________________________________
  /** 
   * Set the internal parameters from the defined options in the
   * Runner object
   * 
   * @param r Runner object 
   */
  void SetOptions(Runner& r)
  {
    TrainSetup::SetOptions(r);
    Option* mom		= r.FindOption("mom");
    Option* ab_what	= r.FindOption("afterburner");
    Option* ab_type	= r.FindOption("ab-type");
    Option* ab_order	= r.FindOption("ab-order");
    Option* disp        = r.FindOption("disp-vtx");
    Option* out_fmd     = r.FindOption("outlier-fmd");
    Option* out_spd     = r.FindOption("outlier-spd");
    
    if (mom && mom->IsSet())  SetMoments(mom->AsString());
    if (ab_what)              fAddFlow       = ab_what->AsString();
    if (ab_type)              fAddFType      = ab_type->AsInt();
    if (ab_order)             fAddFOrder     = ab_type->AsInt();
    if (disp)                 fDispVtx       = disp->AsBool();    
    if (out_fmd)              fOutlierCutFMD = out_fmd->AsDouble();
    if (out_spd)              fOutlierCutSPD = out_spd->AsDouble();
  }
  TString  fType;          // Type of moments to determine 
  Double_t fOutlierCutFMD; // Outlier cut for FMD 
  Double_t fOutlierCutSPD; // Outlier cut for SPD 
  TString  fAddFlow;       // After-burn MC input 
  Int_t    fAddFType;      // How the afterburner is parameterized
  Int_t    fAddFOrder;     // Order of afterburning
  Bool_t   fDispVtx;       // Whether to use satellite interactions
};
//
// EOF
//
