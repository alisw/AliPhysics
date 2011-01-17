/** 
 * Set-up for a PROOF analysis job.   Make TProof object and load pars. 
 * 
 */
void
LoadPars(Int_t nWorkers=4)
{
  const char* option = nWorkers <= 0 ? "" : Form("workers=%d", nWorkers);
  TProof::Open(option);
  const char* pkgs[] = { "STEERBase", "ESD", "AOD", "ANALYSIS", 
			 "ANALYSISalice", "PWG2forward2", 0};
  const char** pkg = pkgs;
  while (*pkg) { 
    gProof->UploadPackage(Form("${ALICE_ROOT}/%s.par",*pkg));
    gProof->EnablePackage(*pkg);    
    pkg++;
  }
}
//
// EOF
//
