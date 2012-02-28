void SETUP() {
  CheckLoadLibrary("libPWGLFforward2");

  // Set the include paths
  gROOT->ProcessLine(".include PWGLFforward2/FORWARD/analysis2");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGLFforward2_INCLUDE", "PWGLFforward2/FORWARD/analysis2");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
