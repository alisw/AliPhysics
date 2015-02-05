/**
 * @file   MakeFlow.C
 * @author Alexander Hansen 
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_scripts_makers
 * 
 */
//====================================================================
/**
 * Script to analyse AOD input for flow
 * 
 * @param name     Name of train 
 * @param options  Options @see RunTrain 
 * @param mode     Which execution environment 
 * @param datadir  Data directory/proof path 
 * @param urlOpts  URL options
 *
 * @ingroup pwglf_forward_flow
 */
void MakeFlow(TString name    = "flow", 
	      TString options = "help",
	      TString mode    = "lite",
	      TString datadir = "",
	      TString urlOpts = "workers=10&recursive") 
{
  if (name.IsNull()) Fatal("MakeFlow", "Must specify a name");
  gROOT->SetMacroPath(Form("%s:$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/trains",
			   gROOT->GetMacroPath()));
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/trains/RunTrain.C");

  if (!datadir.EndsWith("/") && !mode.Contains("proof")) datadir.Append("/");
  
  TUrl url(datadir.Data());
  url.SetProtocol(mode.Data());
  url.SetAnchor("aodTree");
  url.SetOptions(urlOpts.Data());

  RunTrain(name, "MakeFlowTrain", url, options);
}
//--------------------------------------------------------------------
//
// EOF
//
