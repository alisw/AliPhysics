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
  CheckLoadLibrary("libPWG2unicor");

  gROOT->ProcessLine(".include PWG2unicor/UNICOR");
  gSystem->Setenv("PWG2unicor_INCLUDE", "PWG2unicor/UNICOR");
}
//=============================================================================
Int_t CheckLoadLibrary(const char* library) {

  // load library if not yet done

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0) return 1;  
  return gSystem->Load(library);
}
//=============================================================================
