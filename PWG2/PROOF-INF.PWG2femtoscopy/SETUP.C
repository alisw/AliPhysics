void SETUP() {
  CheckLoadLibrary("libPWG2femtoscopyUser");

  // Set the include paths
  gROOT->ProcessLine(".include PWG2femtoscopyUser");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWG2femtoscopyUser_INCLUDE", "PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
