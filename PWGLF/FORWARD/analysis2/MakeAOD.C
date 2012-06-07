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
void MakeAOD(TString     name       = "aod", 
	     TString     options    = "help",
	     TString     runs       = "",
	     Int_t       nEvents    = -1)
{
  if (name.IsNull()) Fatal("MakeAOD", "Must specify a name");
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains/RunTrain.C");

  RunTrain("MakeAODTrain", name, options, runs, nEvents);
}
//
// EOF
//
