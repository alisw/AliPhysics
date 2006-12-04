void SETUP() {
  CheckLoadLibrary("libPWG2ebye");

  // Set the include paths
  gROOT->ProcessLine(".include PWG2ebye");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWG2ebye_INCLUDE", "PWG2ebye");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
