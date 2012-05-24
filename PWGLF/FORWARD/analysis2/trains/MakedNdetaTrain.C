#include "TrainSetup.C"

//====================================================================
/**
 * Analysis train to make @f$ dN/d\eta@f$
 * 
 * To run, do 
 * @code 
 * gROOT->LoadMacro("TrainSetup.C");
 * // Make train 
 * MakedNdetaTrain t("My Analysis");
 * // Set variaous parameters on the train 
 * t.SetDataDir("/home/of/data");
 * t.AddRun(118506)
 * // Run it 
 * t.Run("LOCAL", "FULL", -1, false, false);
 * @endcode 
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
   * @param trig     Trigger to use 
   * @param vzMin    Least @f$ v_z@f$
   * @param vzMax    Largest @f$ v_z@f$
   * @param scheme   Normalisation scheme 
   * @param useCent  Whether to use centrality 
   * @param dateTime Append date and time to name 
   * @param year     Year     - if not specified, current year
   * @param month    Month    - if not specified, current month
   * @param day      Day      - if not specified, current day
   * @param hour     Hour     - if not specified, current hour
   * @param min      Minutes  - if not specified, current minutes
   */
  MakedNdetaTrain(const char* name, 
		  const char* trig="INEL", 
		  Double_t    vzMin=-10, 
		  Double_t    vzMax=10, 
		  const char* scheme="FULL", 
		  Bool_t      useCent=false,
		  Bool_t      dateTime=false,
		  UShort_t    year  = 0, 
		  UShort_t    month = 0, 
		  UShort_t    day   = 0, 
		  UShort_t    hour  = 0, 
		  UShort_t    min   = 0) 
    : TrainSetup(name, dateTime, year, month, day, hour, min),
      fTrig(trig), 
      fVzMin(vzMin), 
      fVzMax(vzMax),
      fScheme(scheme),
      fUseCent(useCent)
  {}
  /** 
   * Run this analysis 
   * 
   * @param mode     Mode - see TrainSetup::EMode
   * @param oper     Operation - see TrainSetup::EOperation
   * @param nEvents  Number of events (negative means all)
   * @param usePar   If true, use PARs 
   */
  void Run(const char* mode, const char* oper, 
	   Int_t nEvents=-1, Bool_t usePar=false, Int_t dbg=0)
  {
    Exec("AOD", mode, oper, nEvents, false, usePar, dbg);
  }
  /** 
   * Run this analysis 
   * 
   * @param mode     Mode - see TrainSetup::EMode
   * @param oper     Operation - see TrainSetup::EOperation
   * @param nEvents  Number of events (negative means all)
   * @param usePar   If true, use PARs 
   */
  void Run(EMode mode, EOper oper, Int_t nEvents=-1, 
	   Bool_t usePar=false, Int_t dbg=0)
  {
    Exec(kAOD, mode, oper, nEvents, false, usePar, dbg);
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
  void SetScheme(const char* scheme) { fScheme = scheme; }
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
   * @param mode Processing mode
   * @param par  Whether to use par files 
   */
  void CreateTasks(EMode mode, Bool_t par, AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_dndeta.root");

    // --- Load libraries/pars ---------------------------------------
    LoadLibrary("PWGLFforward2", mode, par, true);
    
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
  void WriteConstruction(std::ostream& o, const char* obj) const
  {
    o << "  const char* trig   = \"" << fTrig << "\";\n"
      << "  Double_t    vzMin  = " << fVzMin << ";\n"
      << "  Double_t    vzMax  = " << fVzMax << ";\n"
      << "  const char* scheme = \"" << fScheme << "\";\n"
      << "  Bool_t      cent   = " << fUseCent << ";\n"
      << "  MakedNdetaTrain " << obj << "(\"" 
      << fName << "\",trig,vzMin,vzMax,scheme,cent);\n";
  }
  void WriteRun(std::ostream& o, 
		const char* obj, 
		const char* /* type */, 
		const char* mode, 
		const char* oper, 
		Bool_t      mc, 
		Bool_t      usePar, 
		Int_t       dbg) const
  {
    o << "  " << obj << ".Run(" << mode << ',' << oper << ",-1," 
      << usePar << ',' << dbg << ");" << std::endl;
  }
  TString  fTrig;      // Trigger to use 
  Double_t fVzMin;     // Least v_z
  Double_t fVzMax;     // Largest v_z
  TString  fScheme;    // Normalisation scheme 
  Bool_t   fUseCent;   // Use centrality 
};
//
// EOF
//
