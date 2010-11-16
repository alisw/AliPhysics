/** 
 * Run first pass of the analysis - that is read in ESD and produce AOD
 * 
 * @param file           ESD input file
 * @param nEvents        Number of events to process
 * @param nCutBins       Number of additional bins to cut off
 * @param correctionCut  Threshold for when to use secondary map 
 *
 * @ingroup pwg2_forward_analysis_scripts
 */
void
Pass1(const char* file="AliESDs.root", 
      Int_t       nEvents=1000, 
      Int_t       nCutBins=1, 
      Int_t       correctionCut=0.1)
{
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/RunManager.C"); 

  RunManager(file, kFALSE, nEvents, nCutBins, correctionCut);
}
//
// EOF
//
