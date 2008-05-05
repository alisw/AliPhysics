void SETUP() {
  // Load some ROOT libraries
  CheckLoadLibrary("libEG");
  CheckLoadLibrary("libGeom");
  
  // Load the ESD library
  CheckLoadLibrary("libESD");

  // Load the ESD library
  CheckLoadLibrary("libAOD");

  // Load the ESD library
  CheckLoadLibrary("libANALYSIS");

  // Load the PWG2 library
  CheckLoadLibrary("libPWG2AOD");

  // Set the include paths
  gROOT->ProcessLine(".include PWG2AOD/AOD");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWG2AOD_INCLUDE", "PWG2AOD/AOD");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library
  
  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
