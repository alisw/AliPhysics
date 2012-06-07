/**
 * @file   MakeFlow.C
 * @author Alexander Hansen 
 * @date   Wed Mar 23 12:11:33 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_scripts_makers
 * 
 */
//====================================================================
/**
 * Script to analyse AOD input for flow
 * 
 * Takes either a single (AOD) .root file as input or a .txt
 * The .txt file is expected to contain the path to the files 
 * from the current directory or the absolute path.
 * 
 * @par Inputs: 
 *  
 * 
 * @par Outputs: 
 * - 
 *
 * @param name     Name of train 
 * @param options  Options @see RunTrain 
 * @param runs     Options @see RunTrain 
 * @param nEvents  Number of events to process, negative for all 
 *
 * @ingroup pwglf_forward_flow
 */
void MakeFlow(TString name = "flow", 
	      TString options = "help", 
	      TString runs = "", 
	      Int_t   nEvents   = -1)
{
  if (name.IsNull()) Fatal("MakeFlow", "Must specify a name");
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains/RunTrain.C");

  RunTrain("MakeFlowTrain", name, options, runs, nEvents);
}
//--------------------------------------------------------------------
//
// EOF
//
