/**
 * @file   MakeELossFits.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:08:14 2011
 * 
 * @brief  Generate energy loss fits 
 * 
 * @ingroup pwglf_forward_scripts_makers
 */
/** 
 * Run a pass on ESD data to produce the energ loss fits 
 *
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 * @param options    Options string
 * @param url        Execution and input URL
 *
 * @ingroup pwglf_forward_aod
 */
void MakeELossFits(TString     name       = "aod", 
	           TString     url        = "help",
	           TString     options    = "help")
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/trains/MakeTrain.C");

  MakeTrain(name, "MakeFMDELossTrain", url, options);
}
//
// EOF
//
