void SETUP() {
  CheckLoadLibrary("libPWG2flow");
  
  // Set the include paths
  gROOT->ProcessLine(".include PWG2flow/FLOW");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWG2flow_INCLUDE", "PWG2flow/FLOW");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
