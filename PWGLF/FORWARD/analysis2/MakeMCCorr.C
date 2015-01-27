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
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 * @param options    Options string
 * @param url        Execution and input URL
 *
 * @ingroup pwglf_forward_aod
 */
void MakeMCCorr(TString     name       = "aod", 
		TString     url        = "help",
		TString     options    = "help")
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/trains/MakeTrain.C");

  MakeTrain(name, "MakeMCCorrTrain", url, options);
}
//
// EOF
//
