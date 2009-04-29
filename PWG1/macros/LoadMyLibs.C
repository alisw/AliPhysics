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
  CheckLoadLibrary("libESD");
  CheckLoadLibrary("libCDB");
  CheckLoadLibrary("libRAWDatabase");
  CheckLoadLibrary("libRAWDatarec");
  CheckLoadLibrary("libSTEER");
  CheckLoadLibrary("libRAWDatasim");

  CheckLoadLibrary("libTPCbase");
  CheckLoadLibrary("libTPCrec");
  CheckLoadLibrary("libTPCsim");
  CheckLoadLibrary("libITSbase");
  CheckLoadLibrary("libITSsim");
  CheckLoadLibrary("libITSrec");
  CheckLoadLibrary("libTRDbase");
  CheckLoadLibrary("libTRDsim");
  CheckLoadLibrary("libTRDrec");
  CheckLoadLibrary("libTOFbase");
  CheckLoadLibrary("libTOFrec");
  CheckLoadLibrary("libTOFsim");

  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libAOD");
  CheckLoadLibrary("libANALYSISalice");
  CheckLoadLibrary("libTPCcalib");
  CheckLoadLibrary("libPWG1");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library
  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  return gSystem->Load(library);
}
