//
//
//
Bool_t
LoadDummy()
{
  const char*  alice_root = gSystem->Getenv("ALICE_ROOT");
  const char*  dirs[]     = { "", "include", "ITS", 0 };
  const char** d          = dirs;
  TString newpath("-DCOMPILING=1 ");
  newpath += gSystem->GetIncludePath();
  do {
    TString flag(Form("-I%s/%s", alice_root, *d));
    if (newpath.Index(flag) == TString::kNPOS) {
      std::cerr << "Adding " << flag << std::endl;
      newpath += " ";
      newpath += flag;
    }
  } while (*(++d));

  gSystem->SetIncludePath(newpath.Data());
  std::cout << "Include path is\n\t" << gSystem->GetIncludePath() << std::endl;
  gROOT->LoadMacro("Dummy.C+g");
  if (!gROOT->GetClass("Dummy<AliITSvPPRasymmFMD>")) {
    std::cerr << "Failed to make DummyITS" << std::endl;
    return kFALSE;
  }
  return kTRUE;
}
//
// EOF
//
