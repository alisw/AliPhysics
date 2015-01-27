/**
 * @file   MakeTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Nov 23 02:23:03 2012
 * 
 * @brief  Run a train
 * 
 * 
 * @ingroup pwglf_forward_scripts_makers
 */
//====================================================================
/** 
 * Run a train
 * 
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 * @param cls        Class name
 * @param options    Options string
 * @param url        Execution and input URL
 *
 * @ingroup pwglf_forward_aod
 */
void MakeTrain(TString     name       = "aod", 
	       TString     cls        = "",
	       TString     url        = "help",
	       TString     options    = "help")
{
  if (cls.IsNull())  Fatal("MakeTrain", "Must specify a class name");
  if (name.IsNull()) Fatal(cls, "Must specify a name");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/trains/RunTrain.C");
  gROOT->SetMacroPath(Form("%s:$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/trains",
			   gROOT->GetMacroPath()));

  RunTrain(name, cls, uri, options);
}
//
// EOF
//
