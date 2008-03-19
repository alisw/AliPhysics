void LoadMyLibs() {

  // Load some ROOT libraries
  CheckLoadLibrary("libEG");
  CheckLoadLibrary("libGeom");
  CheckLoadLibrary("libVMC");
  CheckLoadLibrary("libTree");
  CheckLoadLibrary("libGui");
  
  // Load AliRoot libraries
  
  CheckLoadLibrary("libSTEERBase");
  CheckLoadLibrary("libCDB");
  CheckLoadLibrary("libRAWDatabase");
  CheckLoadLibrary("libRAWDatarec");
  CheckLoadLibrary("libESD");
  CheckLoadLibrary("libSTEER");
  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libPWG0base");
  CheckLoadLibrary("libPWG0dep");
  CheckLoadLibrary("libPWG1");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library
  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  return gSystem->Load(library);
}
