void SETUP() {
  CheckLoadLibrary("libPWGLFforward");

  // Set the include paths
  gROOT->ProcessLine(".include PWGLFforward/FORWARD");
  gROOT->ProcessLine(".include PWGLFforward/FORWARD/analysis");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGLFforward_INCLUDE", "PWGLFforward/FORWARD");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
