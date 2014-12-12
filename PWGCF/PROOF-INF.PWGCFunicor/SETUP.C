//=============================================================================
void SETUP() {
  CheckLoadLibrary("libPhysics");
  CheckLoadLibrary("libEG");
  CheckLoadLibrary("libTree");
  CheckLoadLibrary("libVMC");
  CheckLoadLibrary("libSTEERBase");
  CheckLoadLibrary("libESD");
  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libANALYSISalice");
  CheckLoadLibrary("libPWGCFunicor");

  gROOT->ProcessLine(".include PWGCFunicor/FEMTOSCOPY/UNICOR");
  gSystem->Setenv("PWGCFunicor_INCLUDE", "PWGCFunicor/FEMTOSCOPY/UNICOR");
}
//=============================================================================
Int_t CheckLoadLibrary(const char* library) {

  // load library if not yet done

  if (strlen(gSystem->GetLibraries(library, "", kFALSE)) > 0) return 1;  
  return gSystem->Load(library);
}
//=============================================================================
