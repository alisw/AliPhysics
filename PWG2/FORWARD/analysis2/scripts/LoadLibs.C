/** 
 * Load the libraries of PWG2/FORWARD/analsysis2
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
void
LoadLibs()
{
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG2forward2");
}
//
// EOF
//
