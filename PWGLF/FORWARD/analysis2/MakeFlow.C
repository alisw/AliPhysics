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
 * @param mode     Which execution environment 
 * @param datadir  Data directory 
 * @param urlOpts  URL options 
 *
 * @ingroup pwglf_forward_flow
 */
void MakeFlow(TString name    = "flow", 
	      TString options = "help",
	      TString mode    = "lite",
	      TString datadir = "/mnt/Disk2/LHC10h_pass2_flowNoSecCorr/",
	      TString urlOpts = "workers=10&recursive") 
{
  if (name.IsNull()) Fatal("MakeFlow", "Must specify a name");
  gROOT->SetMacroPath(Form("%s:$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains",
			   gROOT->GetMacroPath()));
  
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/trains/RunTrain.C");

  if (!datadir.EndsWith("/")) datadir.Append("/");
  
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
