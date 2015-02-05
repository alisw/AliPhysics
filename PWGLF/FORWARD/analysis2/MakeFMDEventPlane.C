/**
 * @file   MakeFlow.C
 * @author Alexander Hansen 
 * @date   Wed Mar 23 12:11:33 2011
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_scripts_makers
 * 
 */
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
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 * @param options    Options string
 * @param url        Execution and input URL
 *
 * @ingroup pwglf_forward_aod
 */
void MakeFMDEventPlane(TString     name       = "aod", 
		       TString     url        = "help",
		       TString     options    = "help")
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/trains/MakeTrain.C");

  MakeTrain(name, "MakeFMDEventPlaneTrain", url, options);
}
//----------------------------------------------------------------
//
// EOF
//
