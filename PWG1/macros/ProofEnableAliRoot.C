
void ProofEnableAliRoot(const char* location = "/usr/local/grid/AliRoot/v4-05-Release")
{
  // enables a locally deployed AliRoot in a PROOF cluster

  /* executes the following commands on each node:
     gSystem->Setenv("ALICE_ROOT", "/home/alicecaf/ALICE/aliroot-head")
     gSystem->AddIncludePath("/home/alicecaf/ALICE/aliroot-head/include");
     gSystem->SetDynamicPath(Form("%s:%s", gSystem->GetDynamicPath(), "/home/alicecaf/ALICE/aliroot-head/lib/tgt_linux"))
     gROOT->Macro("$ALICE_ROOT/macros/loadlibs.C");
  */

  gProof->Exec(Form("TString str(gSystem->ExpandPathName(\"%s\")); gSystem->Setenv(\"ALICE_ROOT\", str);", location), kTRUE);

  gProof->AddIncludePath(Form("%s/include", location));
  gProof->AddIncludePath(Form("%s/TPC", location));
  gProof->AddIncludePath(Form("%s/PWG1", location));
  gProof->AddIncludePath(Form("%s/ANALYSIS", location));

  gProof->AddDynamicPath(Form("%s/lib/tgt_linuxx8664gcc", location));

  // load all libraries
  gProof->Exec("gROOT->Macro(\"$ALICE_ROOT/PWG1/macros/LoadMyLibs.C\")",kTRUE);
  //gProof->Exec("gROOT->Macro(\"$ALICE_ROOT/macros/loadlibsrec.C\")",kTRUE);
}

 
