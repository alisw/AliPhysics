#include "TrainSetup.C"
/**
 * @file   MakeFullTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:53:43 2012
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_trains
 */

//====================================================================
/**
 * Analysis train to make Forward and Central multiplicity, @f$
 * dN/d\eta@f$, flow and @f$\Psi_R@f$ in one loop over the ESDs 
 *
 * @ingroup pwglf_forward_aod
 * @ingroup pwglf_forward_dndete
 * @ingroup pwglf_forward_flow
 * @ingroup pwglf_forward_trains
 */
class MakeFullTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train (free form)
   */
  MakeFullTrain(const  char* name) 
    : TrainSetup(name),
      fSys(0), 
      fSNN(0), 
      fField(0),
      fUseCent(true), 
      fTrig("INEL"), 
      fVzMin(-10),
      fVzMax(+10), 
      fScheme("FULL"),
      fCutEdges(false),
      fFlowMoments("123456"),
      fFlowAfterburnerWhat("pt eta pid cent"), 
      fFlowAfterburnerType(1), 
      fFlowAfterburnerOrder(2),
      fUseTPCEventPlane(false),
      fDebugLvl(0)
  {
    SetType(kESD);
  }
  /** 
   * Set the collision system 
   * 
   * @param sys 1: pp, 2: PbPb, 3: pPb 
   */
  void SetCollisionSystem(UShort_t sys) { fSys = sys; }
  /** 
   * Set the collision energy in GeV
   * 
   * @param sNN @f$\sqrt{s_{NN}}@f$ for PbPb and pPb, @f$\sqrt{s}@f$
   * for pp
   */
  void SetCollisionEnergy(UShort_t sNN) { fSNN = sNN; }
  /** 
   * Set the L3 magnetic field in kilo Gaus
   * 
   * @param fld Field value 
   */
  void SetL3Field(Short_t fld) { fField = fld; }
  /** 
   * Whether to use the centrality information 
   * 
   * @param use Use it or not
   */
  void SetUseCentrality(Bool_t use) { fUseCent = use; } 
  /** 
   * Set the trigger type to generate the final @f$dN/d\eta@f$ for
   *
   * - INEL
   * - INEL>0
   * - NSD 
   *
   * Only makes sense for pp.  
   * 
   * @param trig Trigger string 
   */
  void SetTriggerType(const char* trig) { fTrig = trig; } 
  /** 
   * Set the range for the interaction point @f$l \le v_z \le h@f$
   * 
   * @param l Lower bound
   * @param h Upper bound 
   */
  void SetVertexRange(Double_t l, Double_t h) { fVzMin = l; fVzMax = h; }
  /** 
   * Set the normalization scheme to use.  Only makes sense for pp. 
   * 
   * @param s Scheme to use 
   */
  void SetNormalizationScheme(const char* s) { fScheme = s; }
  /** 
   * Cut off edges of the acceptance 
   * 
   * @param cut 
   */
  void SetCutEdges(Bool_t cut) { fCutEdges = cut; }
  /** 
   * Set the flow moments to calculate.  A string of or more of the
   * numbers from 1 to 6
   * 
   * @param m 
   */
  void SetFlowMoments(const char* m) { fFlowMoments = m; }
  /** 
   * Set parameters on the MC flow after burner 
   * 
   * @param f  How the flow should be parameterized 
   * @param t  Type of parameterization 
   * @param n  Order of flow to add 
   */
  void SetFlowAfterburner(const char* f, Int_t t, Int_t n) 
  { 
    fFlowAfterburnerWhat = f;
    fFlowAfterburnerType = t;
    fFlowAfterburnerOrder = n;
  }
  /** 
   * If set to true, add TPC event plane task. 
   * 
   * @param use Wheter to include TPC event plane task 
   */
  void SetUseTPCEventPlane(Bool_t use) { fUseTPCEventPlane = use; }
  /** 
   * Set the debug level on Forward objects
   * 
   * @param lvl Debug level
   */
  void SetDebugLevel(Int_t lvl) { fDebugLvl = lvl; }
protected:
  /** 
   * Create the tasks 
   * 
   * @param par  Whether to use par files 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(EMode /*mode*/, Bool_t par, AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2", par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;
    
    // --- Add TPC eventplane task
    if (fUseTPCEventPlane) gROOT->Macro("AddTaskEventplane.C");

    // --- Task to copy header information ---------------------------
    gROOT->Macro("AddTaskCopyHeader.C");

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwardMult.C(%d,%d,%d,%d,%d)", 
		      mc, fSys, fSNN, fField, fDebugLvl));
    AddExtraFile(gSystem->Which(gROOT->GetMacroPath(), "ForwardAODConfig.C"));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskCentralMult.C(%d,%d,%d,%d,%d)", 
		      mc, fSys, fSNN, fField, fDebugLvl));
    AddExtraFile(gSystem->Which(gROOT->GetMacroPath(), "CentralAODConfig.C"));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwarddNdeta.C(\"%s\",%f,%f,%d,\"%s\",%d",
		      fTrig.Data(), fVzMin, fVzMax, fUseCent, 
		      fScheme.Data(), fCutEdges));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskCentraldNdeta.C(\"%s\",%f,%f,%d,\"%s\",%d",
		      fTrig.Data(), fVzMin, fVzMax, fUseCent, 
		      fScheme.Data(), fCutEdges));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskMCTruthdNdeta.C(\"%s\",%f,%f,%d,\"%s\")",
		      fTrig.Data(), fVzMin, fVzMax, fUseCent, fScheme.Data()));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskFMDEventPlane.C(%d)", mc));

    // --- Add MC particle task --------------------------------------
    if (mc) gROOT->Macro("AddTaskMCParticleFilter.C");

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwardFlow.C(\"%s\",%d,\"%s\",%d,%d)", 
		      fFlowMoments.Data(), mc, fFlowAfterburnerWhat.Data(), 
		      fFlowAfterburnerType, fFlowAfterburnerOrder));
    

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
      Fatal("CreatePhysicsSelection", "Couldn't get PhysicsSelection (%p)",ps);

    // --- Ignore trigger class when selecting events.  This means ---
    // --- that we get offline+(A,C,E) events too --------------------
    // ps->SetSkipTriggerClassSelection(true);
  }
  //__________________________________________________________________
  /** 
   * Create the centrality selection only if requested
   * 
   * @param mc  Monte-Carlo truth flag 
   * @param mgr Manager
   */
  void CreateCentralitySelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    if (!fUseCent) return;
    TrainSetup::CreateCentralitySelection(mc, mgr);
  }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeFullTrain"; }
  //__________________________________________________________________
 void MakeOptions(Runner& r) 
  {
    TrainSetup::MakeOptions(r);
    r.Add(new Option("sys",    "Collision system",    "1|2|3"));
    r.Add(new Option("sNN",    "Collision energy",    "GEV"));
    r.Add(new Option("field",  "L3 magnetic field", "kG"));
    r.Add(new Option("cent",   "Use centrality"));
    r.Add(new Option("tep",    "Use TPC event plance"));
    r.Add(new Option("trig",   "Trigger class", "INEL|INEL>0|NSD"));
    r.Add(new Option("vzmin",  "Low cut on IP z", "CENTIMETER"));
    r.Add(new Option("vzmax",  "High cut on IP z", "CENTIMETER"));
    r.Add(new Option("scheme", "Normalization scheme", "NONE|FULL"));
    r.Add(new Option("mom", "Flow moments to analyse", "123456"));
    r.Add(new Option("afterburner", "What to afterburn", "[eta,phi,b,pid]"));
    r.Add(new Option("ab-type", "Type of afterburner", "1|2|3|4"));
    r.Add(new Option("ab-order", "Order of afterburner", "1|2|3|4|5|6"));
  }
  //__________________________________________________________________
  void SetOptions(Runner& r)
  {
    TrainSetup::SetOptions(r);
    Option* sys		= r.FindOption("sys");
    Option* sNN		= r.FindOption("sNN");
    Option* field	= r.FindOption("field");
    Option* cent	= r.FindOption("cent");
    Option* tep		= r.FindOption("tep");
    Option* trig	= r.FindOption("trig");
    Option* vzmin	= r.FindOption("vzmin");
    Option* vzmax	= r.FindOption("vzmax");
    Option* scheme	= r.FindOption("scheme");
    Option* mom		= r.FindOption("mom");
    Option* ab_what	= r.FindOption("afterburner");
    Option* ab_type	= r.FindOption("ab-type");
    Option* ab_order	= r.FindOption("ab-order");

    if (sys)   SetCollisionSystem(sys->AsInt());
    if (sNN)   SetCollisionEnergy(sNN->AsInt());
    if (field) SetL3Field(field->AsInt());
    if (cent)  SetUseCentrality(cent->AsBool());
    if (tep)   SetUseTPCEventPlane(tep->AsBool());
    if (trig   && trig->IsSet())   fTrig     = trig->AsString();
    if (vzmin  && vzmin->IsSet())  fVzMin    = vzmin->AsDouble();
    if (vzmax  && vzmax->IsSet())  fVzMax    = vzmax->AsDouble();
    if (scheme && scheme->IsSet()) fScheme   = scheme->AsString();
    if (mom && mom->IsSet())  	   SetFlowMoments(mom->AsString());
    if (ab_what)              fFlowAfterburnerWhat  = ab_what->AsString();
    if (ab_type)              fFlowAfterburnerType  = ab_type->AsInt();
    if (ab_order)             fFlowAfterburnerOrder = ab_type->AsInt();
  }

  UShort_t  fSys;
  UShort_t  fSNN;
  Short_t   fField;
  Bool_t    fUseCent;
  TString   fTrig;
  Double_t  fVzMin;
  Double_t  fVzMax;
  TString   fScheme;
  Bool_t    fCutEdges;
  TString   fFlowMoments;
  TString   fFlowAfterburnerWhat;
  Int_t     fFlowAfterburnerType;
  Int_t     fFlowAfterburnerOrder;
  Bool_t    fUseTPCEventPlane;
  Int_t     fDebugLvl;
};
//
// EOF
//
