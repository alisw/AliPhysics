AliInputEventHandler* AddMCGeneratorHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) ::Fatal("handlers", "No analysis manager");
  AliAnalysisAlien *plugin = (AliAnalysisAlien*)mgr->GetGridHandler();
  if (!plugin) ::Fatal("handlers", "The method should be called via: AliAnalysisAlien::CreateAnalysisManager()");
  mgr->SetInputEventHandler(new AliDummyHandler());
  AliMCGenHandler* mcInputHandler = new AliMCGenHandler();

  TString macroPath;
  macroPath.Form("$ALICE_ROOT/%s", gSystem->Getenv("GEN_MACRO_PATH"));
  macroPath = gSystem->ExpandPathName(macroPath.Data());

  TString macroParameters(gSystem->Getenv("GEN_PARAMETERS"));

  TString newlibs;
  Long64_t retval = AliAnalysisAlien::RunMacroAndExtractLibs(macroPath, macroParameters, newlibs);
  if (retval<0) {
    ::Error(Form("The macro %s did not return a valid generator", macroPath.Data()));
    return;
  }
  AliGenerator *gener = reinterpret_cast<AliGenerator*>(retval);

  // customization from LEGO train
  gROOT->LoadMacro("generator_customization.C");
  generator_customization(gener);

  mcInputHandler->SetGenerator(gener);
  mcInputHandler->SetSeedMode(3);

  newlibs += " ";
  newlibs += gSystem->Getenv("GEN_LIBRARIES");
  plugin->SetGeneratorLibs(newlibs);
  mgr->SetMCtruthEventHandler(mcInputHandler);
  
  return mcInputHandler;
}
