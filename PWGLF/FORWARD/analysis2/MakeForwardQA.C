/**
 * @file   MakeForwardQA.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:08:14 2011
 * 
 * @brief  Generate energy loss fits 
 * 
 * @ingroup pwglf_forward_scripts_makers
 */
/** 
 * Run a pass on ESD data to produce quality assurance data.
 *
 * @param name    Name of train
 * @param options Comma separated list of options, pass "help" for list
 * @param runs    Comma separated list of run numbers
 * @param nEvents Number of events to process.  If 0 or less, then 
 *                all events are analysed
 *
 * @ingroup pwglf_forward
 */
/** 
 * 
 * 
 */void MakeForwardQA(TString name    = "qa", 
		      TString options = "help", 
		      TString runs    = "", 
		      Int_t   nEvents = -1)
{
  if (name.IsNull()) Fatal("MakeForwardQA", "Must specify a name");
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains/RunTrain.C");

  RunTrain("MakeQATrain", name, options, runs, nEvents);
}
//
// EOF
//
