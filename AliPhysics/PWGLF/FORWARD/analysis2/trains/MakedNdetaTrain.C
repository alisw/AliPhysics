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
#include <TF1.h>

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
    fOptions.Add("trig",          "TYPE",      "Trigger type",       "INEL");
    fOptions.Add("filter",        "TYPE",   "Filter type","OUTLIER|PILEUP-BIN");
    fOptions.Add("cent",          "ESTIMATOR", "Use centrality",     "none");
    fOptions.Add("cent-bins",     "BINS",      "Centrality bins",    "default");
    fOptions.Add("ipz-bins",      "BINS",      "IPz bins",           "u15");
    fOptions.Add("scheme",     "SCHEME","Normalization scheme","EVENT,TRIGGER");
    fOptions.Add("trigEff",       "EFFICIENCY","Trigger efficiency", 1.);
    fOptions.Add("trigEff0",      "EFFICIENCY","0-bin trigger effeciency", 1.);
    fOptions.Add("mc",            "Also make dN/deta for MC truth",  false);
    fOptions.Add("satellite",    "Restrict analysis to satellite events",false);
    fOptions.Add("forward-config","FILE", "Forward configuration",
		 "dNdetaConfig.C");
    fOptions.Add("central-config","FILE", "Central configuration", 
		 "dNdetaConfig.C");
    fOptions.Add("truth-config",  "FILE", "MC-Truth configuration", 
		 "dNdetaConfig.C");
    fOptions.Add("mean-ipz",      "MU", "Mean of IPz dist.",         0);
    fOptions.Add("var-ipz",       "SIGMA", "Variance of IPz dist.",  -1);
    fOptions.Add("no-central",    "Do not at SPD cluster task",      false);
    fOptions.Add("abs-min-cent",  "PERCENT", "Absolute least centrality", -1.);
  }
protected:
  Bool_t CoupledNdetaCar(const char* which,
			 const char* cfg)
  {
    UInt_t mask = AliVEvent::kAny;
    AliAnalysisTaskSE* tsk = CoupleSECar("AddTaskdNdeta.C",
					 Form("\"%s\",\"%s\"", which, cfg),
					 mask);
    if (!tsk) {
      Printf("Failed to add task via AddTaskdNdeta.C(%s,%s)", cfg);
      return false;
    }
    FromOption(tsk, "TriggerMask",         "trig",        "INEL");
    FromOption(tsk, "FilterMask",          "filter",      "OUTLIER|PILEUP-BIN");
    FromOption(tsk, "NormalizationScheme", "scheme",      "EVENT,TRIGGER");
    FromOption(tsk, "CentralityMethod",    "cent",        "default");
    FromOption(tsk, "CentralityAxis",      "cent-bins",   "default");
    FromOption(tsk, "IPzAxis",             "ipz-bins",    "u10");
    FromOption(tsk, "TriggerEff",          "trigEff",     1.);
    FromOption(tsk, "TriggerEff0",         "trigEff0",    1.);
    FromOption(tsk, "SatelliteVertices",   "satellite",   false);
    FromOption(tsk, "AbsMinCent",          "abs-min-cent",-1.);

    if (!TString(which).BeginsWith("forward",TString::kIgnoreCase)) return true;

    Double_t muIpz    = fOptions.AsDouble("mean-ipz",0);
    Double_t sigmaIpz = fOptions.AsDouble("var-ipz",-1);
    if (sigmaIpz <= 0) return true;

    TF1* f=new TF1("ipZw","TMath::Gaus(x,[0],[1],1)/TMath::Gaus(x,[2],[3],1)");
    f->SetParNames("#mu_{emp}","#sigma_{emp}","#mu_{this}","#sigma_{this}");
    f->SetParameters(0.592,6.836,muIpz,sigmaIpz);
    // f->SetParErrors(0.023,0.029,eMuIpz,eSigmaIpz);

    Printf("Created re-weight function");
    f->Print();

    gROOT
      ->ProcessLine(Form("((AliBasedNdetaTask*)%p)->SetIpzReweight((TF1*)%p);",
			 tsk, f));
    return true;

  }
  /** 
   * Create the tasks 
   * 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward_dndeta.root");

    // --- Load libraries/pars ---------------------------------------
    fRailway->LoadLibrary("PWGLFforward2");
    // fRailway->LoadLibrary("AOD");
    // gSystem->ListLibraries();
    // gSystem->Load("libAOD");

    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));

    // --- Get parameters --------------------------------------------
    Bool_t   mc      = fOptions.Has("mc");
    TString  fwdCfg  = fOptions.Get("forward-config");
    TString  cenCfg  = fOptions.Get("central-config");
    TString  mcCfg   = fOptions.Get("truth-config");
    if (!mc) mc      = fRailway->IsMC(); 

    // --- Add the task ----------------------------------------------
    CoupledNdetaCar("Forward", fwdCfg);
    Bool_t noCentral = fOptions.Has("no-central");
    if (!noCentral)
      CoupledNdetaCar("Central", cenCfg);
    if (mc) CoupledNdetaCar("MCTruth", mcCfg);

  }
  //__________________________________________________________________
  /** 
   * Do not the centrality selection
   */
  //__________________________________________________________________
  // void CreateCentralitySelection(Bool_t) {}
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

    SaveDraw();
    SaveSummarize();
  }
  void SaveSummarize()
  {
    std::ofstream f("Summarize.C");
    if (!f) { 
      Error("SaveSummarize", "Failed to open Summarize.C script");
      return;
    }
    f << "// Generated by " << ClassName() << "\n"
      << "// WHAT is a bit mask of\n"
      << "//   0x001     Forward\n"
      << "//   0x002     Central\n"
      << "//   0x004     Sums\n"
      << "//   0x008     Results\n"
      << "//   0x010     Only min-bias (no centrality)\n"
      << "//   0x080     Assume simulation results\n"
      << "//   0x100     Landscape\n"
      << "//   0x200     Pause\n"
      << "//   0x400     Also draw single result canvas\n"
      << "//\n"
      << "void Summarize(const char* filename=\"forward_dndeta.root\",\n"
      << "               UShort_t what=0x10F)\n"
      << "{\n"
      << "  const char* fwd=\"$ALICE_PHYSICS/PWGLF/FORWARD/analysis2\";\n"
      << "  gROOT->LoadMacro(Form(\"%s/DrawdNdetaSummary.C\",fwd));\n"
      << "  DrawdNdetaSummary(filename,what & 0x3FF);\n"
      << "\n"
      << "  if (!(what & 0x400)) return;\n"
      << "  gROOT->SetMacroPath(Form(\"../:%s\",gROOT->GetMacroPath()));\n"
      << "  gROOT->Macro(\"Draw.C\");\n"
      << "}\n"
      << "// EOF" << std::endl;
    f.close();
  }
  /** 
   * Make a ROOT script to draw results 
   * 
   */
  void SaveDraw()
  {
    std::ofstream o("Draw.C");
    if (!o) { 
      Error("MakedNdetaTrain::SaveSetup", "Failed to open Draw.C");
      return;
    }

    o << "// Created by " << ClassName() << "\n"
      << "// Will draw dN/deta results from produced file\n"
      << "// \n"
      << "void Draw(UShort_t rebin=5,\n"
      << "          Double_t eff=1,\n"
      << "          Bool_t   raw=false,\n"
      << "          Bool_t   cutEdges=false,\n"
      << "          const char* which=\"Forward\")\n"
      << "{\n"
      << "  TString fwd=\"$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta\";\n"
      << "  gROOT->LoadMacro(Form(\"%s/ExtractGSEs.C\",fwd.Data()));\n"
      << "  ExtractGSEs(\"forward_dndeta.root\",\n"
      << "              rebin,eff,raw,cutEdges,which);\n"
      << "}\n"
      << "//\n"
      << "// EOF\n"
      << "//" << std::endl;
    o.close();
  }
  void PostShellCode(std::ostream& f)
  {
    f << "  echo \"=== Summarizing results ...\"\n"
      << "  aliroot -l -b -q ${prefix}Summarize.C\n"
      << "  echo \"=== Draw results ...\"\n"
      << "  aliroot -l -b -q ${prefix}Draw.C\n"
      << std::endl;
  }
};
//
// EOF
//
