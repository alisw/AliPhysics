/**
 * @file   MakeMCCorr.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Jul 12 10:06:07 2011
 * 
 * @brief  Generate MC corrections 
 * 
 * @ingroup pwglf_forward_scripts_makers
 */
//====================================================================
/** 
 * Generate MC corrections 
 * 
 * @param name     Name of train
 * @param options  Comma separated list options, pass "help" for list
 * @param runs     Comma separated list of run numbers
 * @param nEvents  Number of events to process.  If 0 or less, then 
 *                 all events are analysed
 *
 * @ingroup pwglf_forward_corr
 */
void MakeMCCorr(TString name    = "mcCorr", 
		TString options = "help", 
		TString runs    = "", 
		Int_t   nEvents = -1)
{
  if (name.IsNull()) Fatal("MakeMCCorr", "Must specify a name");
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains/RunTrain.C");

  RunTrain("MakeMCCorrTrain", name, options, runs, nEvents);
}
//
// EOF
//
