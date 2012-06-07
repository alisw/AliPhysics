/**
 * @file   MakedNdetaTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:51:26 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains
 * 
 */

#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ dN/d\eta@f$
 * 
 *
 * @ingroup pwglf_forward_dndeta
 * @ingroup pwglf_forward_trains
 */
class MakedNdetaTrain : public TrainSetup
{
public:
  /** 
   * Constructor.  Date and time must be specified when running this
   * in Termiante mode on Grid
   * 
   * @param name     Name of train (free form)
   */
  MakedNdetaTrain(const char* name)
  : TrainSetup(name),
      fTrig("INEL"), 
      fVzMin(-10), 
      fVzMax(+10),
      fScheme(""),
      fUseCent(false)
  {
    SetType(kAOD);
  }
  /** 
   * Set the trigger to use (INEL, INEL>0, NSD)
   * 
   * @param trig Trigger to use 
   */
  void SetTrigger(const char* trig) { fTrig = trig; }
  /** 
   * Set the vertex range to accept 
   * 
   * @param min Minimum 
   * @param max Maximum 
   */
  void SetVertexRange(Double_t min, Double_t max) { fVzMin=min; fVzMax=max; }
  /** 
   * Set the normalisation scheme 
   * 
   * @param scheme Normalisation scheme options 
   */
  void SetNormalizationScheme(const char* scheme) { fScheme = scheme; }
  /** 
   * Whether to use centrality or not 
   * 
   * @param use To use the centrality 
   */
  void SetUseCentrality(Bool_t use) { fUseCent = use; }
protected:
  /** 
   * Create the tasks 
   * 
   * @param par  Whether to use par files 
   */
  void CreateTasks(EMode /*mode*/, Bool_t par, AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_dndeta.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2", par, true);
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwarddNdeta.C(\"%s\",%f,%f,%d,\"%s\")",
		      fTrig.Data(), fVzMin, fVzMax, fUseCent, fScheme.Data()));

    gROOT->Macro(Form("AddTaskCentraldNdeta.C(\"%s\",%f,%f,%d,\"%s\")",
		      fTrig.Data(), fVzMin, fVzMax, fUseCent, fScheme.Data()));

    gROOT->Macro(Form("AddTaskMCTruthdNdeta.C(\"%s\",%f,%f,%d,\"%s\")",
		      fTrig.Data(), fVzMin, fVzMax, fUseCent, fScheme.Data()));
  }
  /** 
   * Do not the centrality selection
   */
  void CreateCentralitySelection(Bool_t, AliAnalysisManager*) {}
  /** 
   * Crete output handler - we don't want one here. 
   * 
   * @return 0
   */
  AliVEventHandler* CreateOutputHandler(EType) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "MakedNdetaTrain"; }
  //__________________________________________________________________
  void MakeOptions(Runner& r) 
  {
    TrainSetup::MakeOptions(r);
    r.Add(new Option("trig", "Trigger class", "INEL|INEL>0|NSD"));
    r.Add(new Option("vzmin", "Low cut on IP z", "CENTIMETER"));
    r.Add(new Option("vzmax", "High cut on IP z", "CENTIMETER"));
    r.Add(new Option("scheme", "Normalization scheme", "NONE|FULL"));
    r.Add(new Option("cent",   "Use centrality"));
  }
  //__________________________________________________________________
  void SetOptions(Runner& r)
  {
    TrainSetup::SetOptions(r);
    Option* trig	= r.FindOption("trig");
    Option* vzmin	= r.FindOption("vzmin");
    Option* vzmax	= r.FindOption("vzmax");
    Option* scheme	= r.FindOption("scheme");
    Option* cent	= r.FindOption("cent");

    if (trig   && trig->IsSet())   fTrig     = trig->AsString();
    if (vzmin  && vzmin->IsSet())  fVzMin    = vzmin->AsDouble();
    if (vzmax  && vzmax->IsSet())  fVzMax    = vzmax->AsDouble();
    if (scheme && scheme->IsSet()) fScheme   = scheme->AsString();
    if (cent)                      fUseCent  = cent->AsBool();
  }
  //__________________________________________________________________
  TString  fTrig;      // Trigger to use 
  Double_t fVzMin;     // Least v_z
  Double_t fVzMax;     // Largest v_z
  TString  fScheme;    // Normalisation scheme 
  Bool_t   fUseCent;   // Use centrality 
};
//
// EOF
//
