TObject* Browse(Bool_t fmd=true, Bool_t rw=false)
{
  if (gROOT->IsBatch()) { 
    Error("", "Cannot run Forward OADB browser in batch mode");
    return;
  }
  const char* fwd = "${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2";
  if (!gROOT->GetClass("AliOADBForward"))
    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
  gSystem->AddIncludePath(Form("-I%s -I$ALICE_PHYSICS/include", fwd));
  TString macro = Form("%s/corrs/ForwardOADBGui.C+", fwd);
  Info("", "Loading macro %s", macro.Data());
  gROOT->LoadMacro(macro);

  TString fn(fmd ? "fmd_corrections.root" : "spd_corrections.root");
  AliOADBForward* db = new AliOADBForward;
  db->Open(fn, "*", rw);
  
  ForwardOADBGui(db, fn.Data());

  return db;
}
