/**
 * @file   MakedNdeta.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 09:41:56 2011
 * 
 * @brief  Run second pass analysis - make @f$ dN/d\eta@f$
 * 
 * @ingroup pwglf_forward_scripts_makers
 */
//====================================================================
/** 
 * Run second pass analysis - make @f$ dN/d\eta@f$
 * 
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 * @param options    Options string
 * @param url        Execution and input URL
 *
 * @ingroup pwglf_forward_aod
 */
void MakedNdeta(TString     name       = "aod", 
	        TString     url        = "help",
	        TString     options    = "help")
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/trains/MakeTrain.C");

  MakeTrain(name, "MakedNdetaTrain", url, options);
}
//
// EOF
//
