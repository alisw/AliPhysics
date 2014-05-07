 	TObject* Browse()
{
  const char* fwd = "/opt/alice/aliroot/trunk-inst/PWGLF/FORWARD/analysis2";
  if (!gROOT->GetClass("AliOADBForward"))
    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
  gSystem->AddIncludePath(Form("-I%s -I$ALICE_ROOT/include", fwd));
  gROOT->LoadMacro(Form("%s/corrs/ForwardOADBGui.C", fwd));
  
  AliOADBForward* db = new AliOADBForward;
  db->Open("fmd_corrections.root", "*");
  
  ForwardOADBGui(db);

  return db;
}
