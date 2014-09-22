TObject* Browse(Bool_t fmd=true, Bool_t rw=false)
{
  if (gROOT->IsBatch()) { 
    Error("", "Cannot run Forward OADB browser in batch mode");
    return;
  }
  const char* fwd = "${ALICE_ROOT}/PWGLF/FORWARD/analysis2";
  if (!gROOT->GetClass("AliOADBForward"))
    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
  gSystem->AddIncludePath(Form("-I%s -I$ALICE_ROOT/include", fwd));
  TString macro = Form("%s/corrs/ForwardOADBGui.C+", fwd);
  Info("", "Loading macro %s", macro.Data());
  gROOT->LoadMacro(macro);
  
  AliOADBForward* db = new AliOADBForward;
  db->Open(fmd ? "fmd_corrections.root" : "spd_corrections.root", "*", rw);
  
  ForwardOADBGui(db);

  return db;
}
