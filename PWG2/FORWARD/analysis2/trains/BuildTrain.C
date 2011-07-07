/** 
 * Build (compile) a train script 
 * 
 * @param script Script to compile 
 * @param extra  Extra stuff for AcLic ("", "+", or "+g")
 * @param useTmp Use a temporary file 
 * 
 * @return 
 *
 * @ingroup pwg2_forward_trains
 */
Bool_t
BuildTrain(const char* script, const char* extra="", Bool_t useTmp=false)
{
  // --- Load needed libraries ---------------------------------------
  gROOT->LoadClass("TVirtualMC",           "libVMC");
  gROOT->LoadClass("AliVEvent",            "libSTEERBase");
  gROOT->LoadClass("AliESDEvent",          "libESD");
  gROOT->LoadClass("AliAnalysisManager",   "libANALYSIS");
  gROOT->LoadClass("AliAnalysisTaskSE",    "libANALYSISalice");
  gROOT->LoadClass("AliAODForwardMult",    "libPWG2forward2");

  // --- Setup script path -------------------------------------------
  const char* aliPath   = gSystem->ExpandPathName("$ALICE_ROOT");
  const char* fwd2Path  = 
    gSystem->ExpandPathName("$ALICE_ROOT/PWG2/FORWARD/analysis2");
  const char* macroPath = gROOT->GetMacroPath();
  gROOT->SetMacroPath(".:%s:%s/trains:%s/scripts",macroPath,fwd2Path,fwd2Path);
  
  // --- Setup include path ------------------------------------------
  gSystem->AddIncludePath(Form("-I%s", fwd2Path));
  gSystem->AddIncludePath(Form("-I%s", aliPath));
  gSystem->AddIncludePath(Form("-I%s/include", aliPath));

  // --- Check that we can find the script ---------------------------
  TString path = gSystem->Which(gROOT->GetMacroPath(), script);
  if (path.IsNull()) { 
    path = gSystem->Which(gROOT->GetMacroPath(), Form("%s.C", script));
    if (path.IsNull()) { 
      Error("BuildTrain", "Cannot find script %s or %s.C in %s", 
	    script, script, gROOT->GetMacroPath());
      return false;
    }
  }
  
  // --- Possibly make a temporary copy ------------------------------
  if (useTmp) { 
    TString tmp("Train");
    FILE*   fp = gSystem->TempFileName(tmp, ".");
    fclose(fp);
    gSystem->Unlink(tmp);
    tmp.Append(".C");
    
    gSystem->CopyFile(path, tmp);
    path = tmp;
  }

  // --- Now compile the thing ---------------------------------------
  Int_t ret = gROOT->LoadMacro(Form("%s+%s", path.Data(), extra));

  // --- If we did a temporary file, remove stuff --------------------
  if (useTmp) {
    gSystem->Unlink(path);
    path.ReplaceAll(".C",   "_C.d");  gSystem->Unlink(path);
    path.ReplaceAll("_C.d", "_C.so"); gSystem->Unlink(path);
  }

  // --- Warning in case of problems ---------------------------------
  if (ret != 0) 
    Warning("BuildSetup", "Failed to build %s (%s)", path.Data(), script);

  return ret == 0;
}
//
// EOF
//

