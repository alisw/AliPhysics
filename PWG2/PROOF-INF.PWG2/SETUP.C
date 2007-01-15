void SETUP() {
  // Load some ROOT libraries
  CheckLoadLibrary("libEG");
  CheckLoadLibrary("libGeom");
  
  // Load the ESD library
  CheckLoadLibrary("libESD");

  // Load the PWG2 library
  CheckLoadLibrary("libPWG2");

  // Set the include paths
  gROOT->ProcessLine(".include PWG2");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWG2_INCLUDE", "PWG2");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library
  
  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
