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

  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libAOD");
  CheckLoadLibrary("libANALYSISalice");
  CheckLoadLibrary("libTENDER");
  CheckLoadLibrary("libSTAT");

  CheckLoadLibrary("libTPCbase");
  CheckLoadLibrary("libTPCsim");
  CheckLoadLibrary("libTPCrec");
  CheckLoadLibrary("libTPCcalib");
  CheckLoadLibrary("libITSbase");
  CheckLoadLibrary("libITSsim");
  CheckLoadLibrary("libITSrec");
  CheckLoadLibrary("libTRDcalib");
  CheckLoadLibrary("libTRDbase");
  CheckLoadLibrary("libTRDrec");
  CheckLoadLibrary("libTRDsim");
  CheckLoadLibrary("libTOFbase");
  CheckLoadLibrary("libTOFrec");
  CheckLoadLibrary("libTOFsim");

  CheckLoadLibrary("libPWG1");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library
  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  return gSystem->Load(library);
}
