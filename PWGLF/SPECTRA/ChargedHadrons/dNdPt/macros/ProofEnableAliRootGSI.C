
void ProofEnableAliRootGSI(const char* location = "/usr/local/grid/AliRoot/v4-05-Release")
{
  // enables a locally deployed AliRoot in a PROOF cluster

  printf("Load libraries: %s \n",location);
  gProof->Exec(Form("TString str(gSystem->ExpandPathName(\"%s\")); gSystem->Setenv(\"ALICE_ROOT\", str);", location), kTRUE);

  gProof->AddIncludePath(Form("%s/include", location));
  gProof->AddIncludePath(Form("%s/TPC", location));
  gProof->AddIncludePath(Form("%s/PWGLF", location));
  gProof->AddIncludePath(Form("%s/PWGLF/SPECTRA/ChargedHadrons/dNdPt", location));
  gProof->AddIncludePath(Form("%s/ANALYSIS", location));

  gProof->AddDynamicPath(Form("%s/lib/tgt_linuxx8664gcc", location));

  // load all libraries
  gProof->Exec("gROOT->Macro(\"$ALICE_ROOT/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/LoadMyLibs.C\")",kTRUE);
}
