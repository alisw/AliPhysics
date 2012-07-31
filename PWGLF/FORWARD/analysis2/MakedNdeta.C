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
 * If PROOF mode is selected, then Terminate will be run on the master node 
 * in any case. 
 *
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 * @param options    Option string
 * @param runs       List of runs, or file name of file contain runs
 * @param nEvents    Number of events to process.  If 0 or less, then 
 *                   all events are analysed
 *
 * @ingroup pwglf_forward_aod
 */
void MakedNdeta(TString     name     = "dndeta", 
		TString     options  = "help", 
		TString     runs     = "",
	        Int_t       nEvents  = -1)
{
  if (name.IsNull()) Fatal("MakedNdeta", "Must specify a name");
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains/RunTrain.C");

  RunTrain("MakedNdetaTrain", name, options, runs, nEvents);
}
//
// EOF
//
