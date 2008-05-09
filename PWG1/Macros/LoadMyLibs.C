void LoadMyLibs() {

  // Load some ROOT libraries
  CheckLoadLibrary("libEG");
  CheckLoadLibrary("libGeom");
  CheckLoadLibrary("libVMC");
  CheckLoadLibrary("libTree");
  CheckLoadLibrary("libGui");
  CheckLoadLibrary("libMinuit");
  
  // Load AliRoot libraries
  
  CheckLoadLibrary("libSTEERBase");
  CheckLoadLibrary("libCDB");
  CheckLoadLibrary("libRAWDatabase");
  CheckLoadLibrary("libRAWDatarec");
  CheckLoadLibrary("libESD");
  CheckLoadLibrary("libSTEER");
  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libANALYSISalice");
  CheckLoadLibrary("libTPCbase");
  CheckLoadLibrary("libTPCsim");
  CheckLoadLibrary("libTRDbase");
  CheckLoadLibrary("libTPCrec");
  CheckLoadLibrary("libTRDrec");
  CheckLoadLibrary("libPWG1");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library
  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  return gSystem->Load(library);
}
