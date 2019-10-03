namespace {
  void
  AddPath(const TString& dir, Bool_t prepend=true)
  {
    TString d(gSystem->ExpandPathName(dir.Data()));
    gSystem->AddIncludePath(Form("-I%s", d.Data()));
    const char* oldPath = gROOT->GetMacroPath();
    gROOT->SetMacroPath(Form(".:%s:%s",
			   prepend ? d.Data() : oldPath,
			     prepend ? oldPath  : d.Data()));
  }
}

void
DoOne(UShort_t flags=0x0, const char* var="none", Bool_t forceK=false)
{
  const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
  if (gSystem->Getenv("ANA_SRC"))
    fwd = gSystem->Getenv("ANA_SRC");
  AddPath("$HOME/GraphSysErr");
  AddPath(TString::Format("%s/gse", fwd), false);
  AddPath(TString::Format("%s/dndeta/tracklets", fwd));
  if (!gROOT->GetClass("GraphSysErr")) 
    gROOT->LoadMacro("$HOME/GraphSysErr/GraphSysErr.C+g");
  if (!gROOT->GetListOfGlobals()->FindObject("kCorrectLoaded"))
    gROOT->LoadMacro("Correct.C");
  Correct(flags, "left",   var, forceK);
  Correct(flags, "middle", var, forceK);
  Correct(flags, "right",  var, forceK);
  if (!gROOT->GetListOfGlobals()->FindObject("kCombineLoaded"))
    gROOT->LoadMacro("Combine.C+g");
  Printf("Now combining");
  Combine(flags, var);
}

  
