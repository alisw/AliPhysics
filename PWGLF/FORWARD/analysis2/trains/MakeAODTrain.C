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
    fOptions.Add("sys",   "SYSTEM",  "1:pp, 2:PbPb, 3:pPb", "");
    fOptions.Add("snn",   "ENERGY",  "Center of mass energy in GeV", "");
    fOptions.Add("field", "STRENGTH","L3 field strength in kG", "");
    fOptions.Add("forward-config", "FILE", "Forward configuration", 
		 "ForwardAODConfig.C");
    fOptions.Add("central-config", "FILE", "Forward configuration", 
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
    UShort_t sys = fOptions.AsInt("sys", 0);
    UShort_t sNN = fOptions.AsInt("snn", 0);
    UShort_t fld = fOptions.AsInt("field", 0);
    TString  cor = "";
    if (fOptions.Has("corr")) cor = fOptions.Get("corr"); 
    
    // --- Add the task ----------------------------------------------
    TString fwdConfig = fOptions.Get("forward-config");
    gROOT->Macro(Form("AddTaskForwardMult.C(%d,%d,%d,%d,\"%s\",\"%s\")", 
		      mc, sys, sNN, fld, fwdConfig.Data(), cor.Data()));
    fHelper->LoadAux(gSystem->Which(gROOT->GetMacroPath(), fwdConfig), true);

    // --- Add the task ----------------------------------------------
    TString cenConfig = fOptions.Get("central-config");
    gROOT->Macro(Form("AddTaskCentralMult.C(%d,%d,%d,%d,\"%s\",\"%s\")", 
		      mc, sys, sNN, fld, cenConfig.Data(), cor.Data()));
    fHelper->LoadAux(gSystem->Which(gROOT->GetMacroPath(), cenConfig), true);

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
	      "Cannot make dNdeta.C script with helper");
      return;
    }
    
    OptionList  uopts(fHelper->Options());
    
    TString cls("MakedNdetaTrain");
    TString name(fName); name.Append("_dndeta");
    OptionList opts(fOptions);
    opts.Remove("forward-config");
    opts.Remove("central-config");
    opts.Remove("sys");
    opts.Remove("snn");
    opts.Remove("field");
    opts.Remove("bare-ps");
    opts.Remove("tpc-ep");
    opts.Remove("corr");
    opts.Add("trig", "TRIGGER", "Trigger type");
    opts.Add("vzMin", "CENTIMETER", "Lower bound on Ip Z", -10.);
    opts.Add("vzMax", "CENTIMETER", "Upper bound on Ip Z", +10.);
    opts.Add("scheme", "FLAGS", "Normalization scheme", 
	     "TRIGGER EVENT BACKGROUND");
    opts.Add("cut-edges", "Cut edges of acceptance");
    opts.Add("trigEff", "EFFICIENCY", "Trigger efficiency", 1.);
    opts.Add("trigEff0", "EFFICIENCY", "0-bin trigger efficiency", 1.);

    
    // Rewrite our URL 
    TString outString = fHelper->OutputLocation();
    if (outString.IsNull()) outString = fEscapedName;
    TUrl    outUrl(outString);
    
    if (uopts.Find("pattern")) uopts.Set("pattern", "*/AliAOD.root");
    if (uopts.Find("concat")) uopts.Set("concat", true);

    std::stringstream s;
    uopts.Store(s, "", "&", false, true);
    outUrl.SetOptions(s.str().c_str());
      
    opts.Set("url", outUrl.GetUrl());
    opts.Set("type", "AOD");

    SaveSetupROOT("dNdeta", cls, name, opts, &uopts);
    if (asShellScript) 
      SaveSetupShell("dndeta", cls, name, opts, &uopts);

    if (!fHelper || fHelper->Mode() != Helper::kGrid) return;

    SaveDownloadAODs();
  }
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
      << "  for (Int_t i = 0; i < n; i++) {\n"
      << "     TString path(r->GetKey(i, \"turl\"));\n"
      << "     TString sub(gSystem->BaseName(gSystem->DirName(path)));\n"
      << "     TString out = TString::Format(\"AliAOD_%s.root\",sub.Data());\n"
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
};
//
// EOF
//
