//=============================================================================
void SETUP() {
  CheckLoadLibrary("libPhysics.so");
  CheckLoadLibrary("libEG.so");
  CheckLoadLibrary("libTree.so");
  CheckLoadLibrary("libVMC.so"); 
  CheckLoadLibrary("libSTEERBase.so");
  CheckLoadLibrary("libESD.so");
  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libANALYSISalice");
  CheckLoadLibrary("libPWGCFunicor");

  gROOT->ProcessLine(".include PWGCFunicor/UNICOR");
  gSystem->Setenv("PWGCFunicor_INCLUDE", "PWGCFunicor/UNICOR");
}
//=============================================================================
Int_t CheckLoadLibrary(const char* library) {

  // load library if not yet done

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0) return 1;  
  return gSystem->Load(library);
}
//=============================================================================
