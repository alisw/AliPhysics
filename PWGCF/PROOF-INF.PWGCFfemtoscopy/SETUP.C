void SETUP() {
  CheckLoadLibrary("libPWGCFfemtoscopy");

  // Set the include paths
  gROOT->ProcessLine(".include PWGCFfemtoscopy/FEMTOSCOPY/AliFemto");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGCFfemtoscopy_INCLUDE", "PWGCFfemtoscopy/FEMTOSCOPY/AliFemto");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
