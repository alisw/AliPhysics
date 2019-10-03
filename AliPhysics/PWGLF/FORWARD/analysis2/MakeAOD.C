/** 
 * @defgroup pwglf_forward_scripts_makers Maker scripts 
 *
 * Scripts to make parts of the analysis.  Users can execute these directly. 
 *
 * @ingroup pwglf_forward_scripts
 */
//====================================================================
/**
 * @file   MakeAOD.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 09:40:10 2011
 * 
 * @brief  Run first pass of the analysis - AOD generation
 * 
 * @ingroup pwglf_forward_scripts_makers
 */
//====================================================================
/** 
 * Run first pass of the analysis - that is read in ESD and produce AOD
 * 
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 * @param options    Options string
 * @param url        Execution and input URL
 *
 * @ingroup pwglf_forward_aod
 */
void MakeAOD(TString     name       = "aod", 
	     TString     url        = "help",
	     TString     options    = "help")
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/trains/MakeTrain.C");

  MakeTrain(name, "MakeAODTrain", url, options);
}
//
// EOF
//
