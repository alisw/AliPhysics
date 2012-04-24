void LoadMyLibs() {

  // Load some ROOT libraries
  CheckLoadLibrary("libEG");
  CheckLoadLibrary("libGeom");
  CheckLoadLibrary("libVMC");
  CheckLoadLibrary("libTree");
  CheckLoadLibrary("libGui");
  CheckLoadLibrary("libMinuit");
  CheckLoadLibrary("libSTAT");
  
  // Load AliRoot libraries
  
  CheckLoadLibrary("libSTEERBase");
  CheckLoadLibrary("libCDB");
  CheckLoadLibrary("libRAWDatabase");
  CheckLoadLibrary("libRAWDatarec");
  CheckLoadLibrary("libESD");
  CheckLoadLibrary("libAOD");
  CheckLoadLibrary("libSTEER");
  
  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libANALYSISalice");
  CheckLoadLibrary("libTENDER");
  CheckLoadLibrary("libTPCbase");
  CheckLoadLibrary("libTPCsim");
  CheckLoadLibrary("libTPCrec");
  //CheckLoadLibrary("libTRDbase");
  //CheckLoadLibrary("libTRDcalib");
  //CheckLoadLibrary("libTRDrec");
  CheckLoadLibrary("libITSbase");
  CheckLoadLibrary("libITSrec");

  CheckLoadLibrary("libPWGPP");
  CheckLoadLibrary("libPWGPPdep");
  CheckLoadLibrary("libPWGLFspectra");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not library is loaded
  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  return gSystem->Load(library);
}
