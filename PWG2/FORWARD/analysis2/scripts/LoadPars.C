/** 
 * Set-up for a PROOF analysis job.   Make TProof object and load pars. 
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
Bool_t
LoadPars(Int_t nWorkers=4)
{
  const char* option = nWorkers <= 0 ? "" : Form("workers=%d", nWorkers);
  TProof::Open(option);
  const char* pkgs[] = { "STEERBase", "ESD", "AOD", "ANALYSIS", 
			 "ANALYSISalice", "PWG2forward2", 0};
  const char** pkg = pkgs;
  Bool_t ret = true;
  while (*pkg) { 
    if (gProof->UploadPackage(Form("${ALICE_ROOT}/%s.par",*pkg)) < 0) {
      Error("LoadPars", "Failed to upload package %s", *pkg);
      ret = false;
      continue;
    }
    if (gProof->EnablePackage(*pkg) < 0) { 
      Error("LoadPars", "Failed to enable package %s", *pkg);
      ret = false;
      continue;
    }
    pkg++;
  }
  return ret;
}
//
// EOF
//
