/**
 * @file   MakeAODTrain.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:05:30 2011
 * 
 * @brief  Run first pass analysis - make AOD tree
 * 
 * @ingroup pwglf_forward_trains
 */
#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make Forward and Central multiplicity
 * 
 *
 * @ingroup pwglf_forward_aod
 * @ingroup pwglf_forward_trains
 */
class MakeAODTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train (free form)
   */
  MakeAODTrain(const  char* name) 
    : TrainSetup(name),
      fSys(0), 
      fSNN(0), 
      fField(0),
      fUseCent(true), 
      fUseTPCEventPlane(false),
      fForwardConfig("ForwardAODConfig.C"), 
      fCentralConfig("CentralAODConfig.C")
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
   * If set to true, add TPC event plane task. 
   * 
   * @param use Wheter to include TPC event plane task 
   */
  void SetUseTPCEventPlane(Bool_t use) { fUseTPCEventPlane = use; }
  /** 
   * Whether to process per centrality bin
   * 
   * @param use 
   */
  void SetUseCentrality(Bool_t use) { fUseCent = use; }
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
    gROOT->Macro(Form("AddTaskForwardMult.C(%d,%d,%d,%d,\"%s\")", 
		      mc, fSys, fSNN, fField, fForwardConfig.Data()));
    AddExtraFile(gSystem->Which(gROOT->GetMacroPath(), fForwardConfig));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskCentralMult.C(%d,%d,%d,%d,\"%s\")", 
		      mc, fSys, fSNN, fField, fCentralConfig.Data()));
    AddExtraFile(gSystem->Which(gROOT->GetMacroPath(), fCentralConfig));

    // --- Add MC particle task --------------------------------------
    if (mc) gROOT->Macro("AddTaskMCParticleFilter.C");

  }
  //__________________________________________________________________
  /** 
   * Create physics selection , and add to manager
   * 
   * @param mc Whether this is for MC 
   * @param mgr Manager 
   */
  void CreatePhysicsSelection(Bool_t mc, AliAnalysisManager* mgr)
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
  const char* ClassName() const { return "MakeAODTrain"; }
  //__________________________________________________________________
 void MakeOptions(Runner& r) 
  {
    TrainSetup::MakeOptions(r);
    r.Add(new Option("sys", "Collision system",    "1|2|3"));
    r.Add(new Option("sNN", "Collision energy",    "GEV"));
    r.Add(new Option("field", "L3 magnetic field", "kG"));
    r.Add(new Option("cent",  "Use centrality"));
    r.Add(new Option("tep",   "Use TPC event plance"));
    r.Add(new Option("forward-config","Forward configuration script","SCRIPT"));
    r.Add(new Option("central-config","Central configuration script","SCRIPT"));
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
    Option* fwdConfig   = r.FindOption("forward-config");
    Option* cenConfig   = r.FindOption("central-config");
    
    if (sys)   SetCollisionSystem(sys->AsInt());
    if (sNN)   SetCollisionEnergy(sNN->AsInt());
    if (field) SetL3Field(field->AsInt());
    if (cent)  SetUseCentrality(cent->AsBool());
    if (tep)   SetUseTPCEventPlane(tep->AsBool());
    if (fwdConfig && fwdConfig->IsSet()) 
      fForwardConfig = fwdConfig->AsString(); 
    if (cenConfig && cenConfig->IsSet()) 
      fCentralConfig = cenConfig->AsString(); 
  }
  //__________________________________________________________________
  void SaveOptions(std::ostream& o, const char* str, Runner& r)
  {
    TrainSetup::SaveOptions(o, str, r);
    Option* sys		= r.FindOption("sys");
    Option* sNN		= r.FindOption("sNN");
    Option* field	= r.FindOption("field");
    Option* cent	= r.FindOption("cent");
    Option* tep		= r.FindOption("tep");
    Option* fwdConfig   = r.FindOption("forward-config");
    Option* cenConfig   = r.FindOption("central-config");

    if (sys)   sys->Save(o, str, fSys);
    if (sNN)   sNN->Save(o, str, fSNN);
    if (field) field->Save(o, str, fField);
    if (cent)  cent->Save(o, str, fUseCent);
    if (tep)   tep->Save(o, str, fUseTPCEventPlane);
    if (fwdConfig) fwdConfig->Save(o, str, fForwardConfig);
    if (cenConfig) cenConfig->Save(o, str, fCentralConfig);
    
  }
  //__________________________________________________________________
  void SaveSetupROOT(Runner& r, Int_t nEvents)
  {
    TrainSetup::SaveSetup(r, nEvents);
    
    TString base(Form("dndeta_%s", fEscapedName.Data()));
    
    std::ofstream o(Form("%s.C", base.Data()));
    o << "void " << base << "()\n"
      << "{\n"
      << "  TString opts;\n";
    TrainSetup::SaveOptions(o, "opts", r);
    o << "  // Note, `type' from above overwritten here\n"
      << "  opts.Append(\"type=AOD:\");\n"
      << "  opts.Append(\"trig=INEL:\");\n"
      << "  opts.Append(\"vzMin=-10:\");\n"
      << "  opts.Append(\"vzMax=-10:\");\n"
      << "  opts.Append(\"scheme=full:\");\n"
      << "  opts.Append(\"datadir=" << GetOutputDirectory(fExecMode) 
      << ":\");\n\n"
      << "  TString runs(\"";
    for (Int_t i = 0; i < fRunNumbers.GetSize(); i++) 
      o << (i == 0 ? "" : ", ") << fRunNumbers.At(i);
    o << "\");\n\n"
      << "  Int_t   nEvents = " << nEvents << ";\n\n"
      << "  gROOT->LoadMacro(\"$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains/RunTrain.C\");\n"
      << "  RunTrain(\"MakedNdetaTrain\",\"dndeta_" << fName << "\",opts,runs,nEvents);\n"
      << "}\n"
      << "// EOF" << std::endl;
    o.close();
  }

  UShort_t fSys;                // Pre-set collision system
  UShort_t fSNN;                // Pre-set collision energy 
  Short_t  fField;              // Pre-set magnetic field
  Bool_t   fUseCent;            // Whether to use centrality 
  Bool_t   fUseTPCEventPlane;   // Whether to use TPC event plane 
  TString  fForwardConfig;      // Forward configuration script 
  TString  fCentralConfig;      // Central configuration script 

};
//
// EOF
//
