/**
 * @file   MakeAODTrain.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:05:30 2011
 * 
 * @brief  Run first pass analysis - make AOD tree
 * 
 * @ingroup pwglf_forward_trains_specific
 */
#include "TrainSetup.C"
#include <sstream>

//====================================================================
/**
 * Analysis train to make Forward and Central multiplicity
 * 
 *
 * @ingroup pwglf_forward_aod
 * @ingroup pwglf_forward_trains_specific
 */
class MakeAODTrain : public TrainSetup
{
public:
  /** 
   * Constructor. 
   * 
   * @param name     Name of train (free form)
   */
  MakeAODTrain(const  TString& name) 
    : TrainSetup(name)
  {
    fOptions.Add("run",   "NUMBER",  "Run number for corrs", 0);
    fOptions.Add("sys",   "SYSTEM",  "1:pp, 2:PbPb, 3:pPb", 0);
    fOptions.Add("snn",   "ENERGY",  "Center of mass energy in GeV", 0);
    fOptions.Add("field", "STRENGTH","L3 field strength in kG", 0);
    fOptions.Add("forward-config", "FILE", "Forward configuration", 
		 "ForwardAODConfig.C");
    fOptions.Add("central-config", "FILE", "Central configuration", 
		 "CentralAODConfig.C");
    fOptions.Add("cent",  "Use centrality");
    fOptions.Add("tpc-ep", "Use TPC event plane");
    fOptions.Add("satelitte", "Use satelitte interactions");
    fOptions.Add("corr", "DIR", "Corrections dir", "");
    fOptions.Set("type", "ESD");
  }
protected:
  /** 
   * Create the tasks 
   * 
   * @param mgr  Analysis manager 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    // --- Output file name ------------------------------------------
    AliAnalysisManager::SetCommonFileName("forward.root");

    // --- Load libraries/pars ---------------------------------------
    fHelper->LoadLibrary("PWGLFforward2");
    
    // --- Set load path ---------------------------------------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
			     gROOT->GetMacroPath()));

    // --- Check if this is MC ---------------------------------------
    Bool_t mc = mgr->GetMCtruthEventHandler() != 0;
    
    // --- Add TPC eventplane task
    if (fOptions.Has("tpc-ep")) gROOT->Macro("AddTaskEventplane.C");

    // --- Task to copy header information ---------------------------
    gROOT->Macro("AddTaskCopyHeader.C");

    // --- Get options -----------------------------------------------
    ULong_t  run = fOptions.AsInt("run", 0);
    UShort_t sys = fOptions.AsInt("sys", 0);
    UShort_t sNN = fOptions.AsInt("snn", 0);
    UShort_t fld = fOptions.AsInt("field", 0);
    TString  cor = "";
    if (fOptions.Has("corr")) cor = fOptions.Get("corr"); 
    
    // --- Add the task ----------------------------------------------
    TString fwdConfig = fOptions.Get("forward-config");
    gROOT->Macro(Form("AddTaskForwardMult.C(%d,%ld,%d,%d,%d,\"%s\",\"%s\")", 
		      mc, run, sys, sNN, fld, fwdConfig.Data(), cor.Data()));
    fHelper->LoadAux(gSystem->Which(gROOT->GetMacroPath(), fwdConfig), true);

    // --- Add the task ----------------------------------------------
    TString cenConfig = fOptions.Get("central-config");
    gROOT->Macro(Form("AddTaskCentralMult.C(%d,%ld,%d,%d,%d,\"%s\",\"%s\")", 
		      mc, run, sys, sNN, fld, cenConfig.Data(), cor.Data()));
    fHelper->LoadAux(gSystem->Which(gROOT->GetMacroPath(), cenConfig), true);

    // --- Add MC particle task --------------------------------------
    if (mc) gROOT->Macro("AddTaskMCParticleFilter.C");

    if (!cor.IsNull()) {
      fHelper->LoadAux(Form("%s/fmd_corrections.root",cor.Data()), true);
      fHelper->LoadAux(Form("%s/spd_corrections.root",cor.Data()), true);
    }
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

    // --- Special for pPb pilot run Sep. 2012 -----------------------
    UShort_t sys = fOptions.AsInt("sys", 0);
    if (sys == 3) { 
      Warning("CreatePhysicsSelection", 
	      "Special setup for pPb pilot run September, 2012");
      gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/ANALYSIS/macros",
			       gROOT->GetMacroPath()));
      gROOT->LoadMacro("PhysicsSelectionOADB_CINT5_pA.C");
      gROOT->ProcessLine(Form("((AliPhysicsSelection*)%p)"
			      "->SetCustomOADBObjects("
			      "OADBSelection_CINT5_V0A(),0);", ps));
      ps->SetSkipTriggerClassSelection(true);
    }
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
    if (!fOptions.Has("cent")) return;
    TrainSetup::CreateCentralitySelection(mc, mgr);
  }
  //__________________________________________________________________
  const char* ClassName() const { return "MakeAODTrain"; }
  //__________________________________________________________________
  /** 
   * Overloaded to create new dNdeta.C and dndeta.sh in the output 
   * directory
   * 
   * @param asShellScript 
   */
  void SaveSetup(Bool_t asShellScript)
  {
    TrainSetup::SaveSetup(asShellScript);

    if (!fHelper) { 
      Warning("MakeAODTrain::SaveSetup", 
	      "Cannot make dNdeta.C script without helper");
      return;
    }
    
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    Bool_t              mc  = mgr && (mgr->GetMCtruthEventHandler() != 0);
    OptionList          uopts(fHelper->Options());
    OptionList          opts(fOptions);
    TString             cls("MakedNdetaTrain");
    TString             name(fName);
    Int_t               sys = fOptions.AsInt("sys", 0);
    if (name.Contains("aod")) name.ReplaceAll("aod", "dndeta");
    else                      name.Append("_dndeta");
    opts.Remove("forward-config");
    opts.Remove("central-config");
    opts.Remove("run");
    opts.Remove("sys");
    opts.Remove("snn");
    opts.Remove("field");
    opts.Remove("bare-ps");
    opts.Remove("tpc-ep");
    opts.Remove("corr");
    opts.Add("trig", "TRIGGER", "Trigger type", "INEL");
    opts.Add("vzMin", "CENTIMETER", "Lower bound on Ip Z", -10.);
    opts.Add("vzMax", "CENTIMETER", "Upper bound on Ip Z", +10.);
    opts.Add("scheme", "FLAGS", "Normalization scheme", 
	     "TRIGGER EVENT BACKGROUND");
    opts.Add("cut-edges", "Cut edges of acceptance", true);
    opts.Add("trigEff", "EFFICIENCY", "Trigger efficiency", 1.);
    opts.Add("trigEff0", "EFFICIENCY", "0-bin trigger efficiency", 1.);
    opts.Add("mc", "Also analyse MC truth", false);
    
    // Rewrite our URL 
    TString outString = fHelper->OutputLocation();
    if (outString.IsNull()) outString = fEscapedName;
    TUrl    outUrl(outString);
    
    if (uopts.Find("pattern") && outString.EndsWith("AliAOD.root")) 
      uopts.Set("pattern", "*/AliAOD.root");
    if (uopts.Find("concat")) uopts.Set("concat", true);

    std::stringstream s;
    uopts.Store(s, "", "&", false, true);
    outUrl.SetOptions(s.str().c_str());
      
    opts.Set("url", outUrl.GetUrl());
    opts.Set("type", "AOD");
    if (!fDatimeString.IsNull()) opts.Set("date", fDatimeString);

    if (sys != 1) {
      opts.Set("cent", true);
      opts.Set("trig", "");
      opts.Set("scheme", "");
      SaveSetupROOT("dNdeta", cls, name, opts, &uopts);
      if (asShellScript) 
	SaveSetupShell("dndeta", cls, name, opts, &uopts);
    }
    else {
      name.ReplaceAll("dndeta", "dndeta_inel");
      SaveSetupROOT("dNdetaINEL", cls, name, opts, &uopts);
      if (asShellScript) 
	SaveSetupShell("dndeta_inel", cls, name, opts, &uopts);
      
      name.ReplaceAll("inel", "nsd");
      opts.Set("trig", "V0AND");
      opts.Set("scheme", "EVENT TRIGGER");
      SaveSetupROOT("dNdetaNSD", cls, name, opts, &uopts);
      if (asShellScript) 
	SaveSetupShell("dndeta_nsd", cls, name, opts, &uopts);
      
      name.ReplaceAll("nsd", "inelgt0");
      opts.Set("trig", "INELGT0");
      opts.Set("scheme", "EVENT TRIGGER");
      SaveSetupROOT("dNdetaINELGt0", cls, name, opts, &uopts);
      if (asShellScript) 
	SaveSetupShell("dndeta_inelgt0", cls, name, opts, &uopts);
    }

    SaveSummarize();
    if (!fHelper || fHelper->Mode() != Helper::kGrid) return;

    SaveDownloadAODs();
  }
  /** 
   * Write a ROOT script to draw summary 
   * 
   */
  void SaveSummarize()
  {
    std::ofstream f("Summarize.C");
    if (!f) { 
      Error("SaveSummarize", "Failed to open Summarize.C script");
      return;
    }
    f << "// Generated by " << ClassName() << "\n"
      << "// WHAT is a bit mask of\n"
      << "//   0x001     Event inspector\n"
      << "//   0x002     Sharing filter\n"
      << "//   0x004     Density calculator\n"
      << "//   0x008     Corrector\n"
      << "//   0x010     Histogram collector\n"
      << "//   0x020     Analysis step cartoon\n"
      << "//   0x040     Results\n"
      << "//   0x080     Central\n"
      << "//   0x100     Landscape\n"
      << "//   0x200     Pause\n"
      << "//\n"
      << "void Summarize(const char* filename=\"forward.root\",\n"
      << "               UShort_t what=0x1FF)\n"
      << "{\n"
      << "  const char* fwd=\"$ALICE_ROOT/PWGLF/FORWARD/analysis2\";\n"
      << "  gROOT->LoadMacro(Form(\"%s/DrawAODSummary.C\",fwd));\n"
      << "  DrawAODSummary(filename,what);\n"
      << "}\n"
      << "// EOF" << std::endl;
    f.close();
  }
  /** 
   * Make a ROOT Script to download the generated AODs
   * 
   */
  void SaveDownloadAODs()
  {
    std::ofstream f("DownloadAODs.C");
    if (!f) { 
      Error("SaveDownloadAODs", "Failed to open DownloadAODs.C");
      return;
    }
    f << "// Generated by " << ClassName() << "\n"
      << "void DownloadAODs()\n"
      << "{\n"
      << "  if (!TGrid::Connect(\"alien://\")) {\n"
      << "    Error(\"DownloadAODs\",\"Failed to connect to AliEn\");\n"
      << "    return;\n"
      << "  }\n\n"
      << "  TString dir(\"" << fHelper->OutputPath() << "\");\n"
      << "  TString pat(\"*/AliAOD.root\");\n"
      << "  TGridResult* r = gGrid->Query(dir,pat);\n"
      << "  if (!r) {\n"
      << "    Error(\"DownloadAODs\",\"No result from query\");\n"
      << "    return;\n"
      << "  }\n\n"
      << "  Int_t n = r->GetEntries();\n"
      << "  Printf(\"=== Got a total of %d AOD files\",n);\n"
      << "  for (Int_t i = 0; i < n; i++) {\n"
      << "     TString path(r->GetKey(i, \"turl\"));\n"
      << "     TString dir(gSystem->DirName(path));\n"
      << "     TString sub(gSystem->BaseName(dir));\n"
      << "     TString subsub(gSystem->BaseName(gSystem->DirName(dir)));\n"
      << "     TString out = TString::Format(\"AliAOD_%s_%s.root\",\n"
      << "                                   subsub.Data(),sub.Data());\n"
      << "     Printf(\"=== Getting %s %s (%3d/%3d)\",\n"
      << "            subsub.Data(),sub.Data(),i,n);\n"
      << "     if (!TFile::Cp(path, out)) {\n"
      << "       Warning(\"DownloadAODs\",\"Failed to copy %s -> %s\",\n"
      << "               path.Data(), out.Data());\n"
      << "       continue;\n"
      << "     }\n"
      << "   }\n"
      << "}\n"
      << "// EOF\n"
      << std::endl;
    f.close();
  }   
  void PostShellCode(std::ostream& f)
  {
    f << "  echo \"=== Summarizing results ...\"\n"
      << "  aliroot -l -b -q ${prefix}Summarize.C\n"
      << std::endl;
  }
};
//
// EOF
//
