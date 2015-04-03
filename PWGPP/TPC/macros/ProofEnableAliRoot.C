
void ProofEnableAliRoot(const char* location = "/usr/local/grid/AliRoot/v4-05-Release")
{
  // enables a locally deployed AliRoot in a PROOF cluster
  gProof->Exec(Form("TString str(gSystem->ExpandPathName(\"%s\")); gSystem->Setenv(\"ALICE_PHYSICS\", str);", location), kTRUE);

  gProof->AddIncludePath(Form("%s/include", location));
  gProof->AddIncludePath(Form("%s/TPC", location));
  gProof->AddIncludePath(Form("%s/TRD", location));
  gProof->AddIncludePath(Form("%s/TRD/Cal", location));
  gProof->AddIncludePath(Form("%s/PWGPP", location));
  gProof->AddIncludePath(Form("%s/ANALYSIS", location));

  gProof->AddDynamicPath(Form("%s/lib/tgt_linuxx8664gcc", location));

  // load all libraries
  gProof->Exec("gROOT->Macro(\"$ALICE_PHYSICS/PWGPP/TPC/macros/LoadMyLibs.C\")",kTRUE);
}

 
