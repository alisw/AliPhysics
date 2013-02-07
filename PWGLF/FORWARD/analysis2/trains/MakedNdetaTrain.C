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
    fOptions.Add("corr-empty", "Correct for empty bins (deprecated)");
  }
protected:
  /** 
   * Create the tasks 
   * 
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
    Bool_t   corrEm = fOptions.Has("corr-empty");
    Bool_t   mc     = fOptions.Has("mc");

    // --- Form arguments --------------------------------------------
    TString args;
    args.Form("\"%s\",%f,%f,%d,\"%s\",%d,%g,%g,%d",
	      trig.Data(), vzMin, vzMax, cent, scheme.Data(),
	      edges, effT, effT0, corrEm);
    // --- Add the task ----------------------------------------------
    gROOT->Macro(Form("AddTaskForwarddNdeta.C(%s);", args.Data()));
    gROOT->Macro(Form("AddTaskCentraldNdeta.C(%s);", args.Data()));
    gROOT->Macro(Form("AddTaskMCTruthdNdeta.C(%s);", args.Data()));
  }
  //__________________________________________________________________
  /** 
   * Do not the centrality selection
   */
  //__________________________________________________________________
  void CreateCentralitySelection(Bool_t, AliAnalysisManager*) {}
  /** 
   * Do not create MC input handler 
   * 
   * @return Always null
   */
  AliVEventHandler* CreateMCHandler(UShort_t, bool) { return 0; }
  //__________________________________________________________________
  /** 
   * Crete output handler - we don't want one here. 
   * 
   * @return 0
   */
  AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  //__________________________________________________________________
  const char* ClassName() const { return "MakedNdetaTrain"; }
  //__________________________________________________________________
  /** 
   * Overloaded to create new draw.C 
   * 
   * @param asShellScript 
   */
  void SaveSetup(Bool_t asShellScript)
  {
    TrainSetup::SaveSetup(asShellScript);

    std::ofstream o("Draw.C");
    if (!o) { 
      Error("MakedNdetaTrain::SaveSetup", "Failed to open Draw.C");
      return;
    }

    o << "// Created by " << ClassName() << "\n"
      << "// \n"
      << "// Will draw dN/deta results from produced file\n"
      << "// \n"
      << "// Options can be specified as needed. To get help, pass the\n"
      << "// string \"help\" for the title\n"
      << "// \n"
      << "void Draw(const TString& title="",\n"
      << "          UShort_t       rebin=5,\n"
      << "          UShort_t       others=0x7,\n"
      << "          UShort_t       flags=0xD87,\n"
      << "          UShort_t       sNN=0,\n"
      << "          UShort_t       sys=0,\n"
      << "          UShort_t       trg=0,\n"
      << "          Float_t        vzMin=999,"
      << "          Float_t        vzMax=-999)\n"
      << "{\n"
      << "   gROOT->LoadMacro(\"$ALICE_ROOT/PWGLF/FORWARD/analysis2/"
      << "DrawdNdeta.C++\");\n\n"
      << "  if (title.EqualTo(\"help\",TString::kIgnoreCase)) {\n"
      << "    DrawdNdeta(\"help\"); // Get the help\n"
      << "    return;\n"
      << "  }\n\n"
      << "  DrawdNdeta(\"forward_dndeta.root\",\n"
      << "             title,\n"
      << "             rebin,\n"
      << "             others,\n"
      << "             flags,\n"
      << "             sNN,\n"
      << "             sys,\n"
      << "             trg,\n"
      << "             vzMin,\n"
      << "             vzMax);\n"
      << "}\n"
      << "//\n"
      << "// EOF\n"
      << "//" << std::endl;
    o.close();
  }
};
//
// EOF
//
