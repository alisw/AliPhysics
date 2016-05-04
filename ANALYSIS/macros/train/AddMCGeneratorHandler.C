AliInputEventHandler* AddMCGeneratorHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) ::Fatal("handlers", "No analysis manager");
  
  AliAnalysisAlien *plugin = (AliAnalysisAlien*)mgr->GetGridHandler();
  if (!plugin) ::Fatal("handlers", "The method should be called via: AliAnalysisAlien::CreateAnalysisManager()");

  mgr->SetInputEventHandler(new AliDummyHandler());
  AliMCGenHandler* mcInputHandler = new AliMCGenHandler();

  plugin->SetGeneratorLibs(gSystem->Getenv("GEN_LIBRARIES"));
  mcInputHandler->SetGeneratorMacroPath(gSystem->Getenv("GEN_MACRO_PATH"));
  mcInputHandler->SetGeneratorMacroParameters(gSystem->Getenv("GEN_PARAMETERS"));
  
  TMacro* macro = new TMacro("generator_customization.C");
  mcInputHandler->SetGeneratorCustomization(macro);
  
  mcInputHandler->SetSeedMode(3);

  const char *streamGen = gSystem->Getenv("STREAM_GENERATOR");
  if (streamGen && (atoi(streamGen) == 1)) {
    ::Info("handlers", "creating generator for streaming");
    mcInputHandler->CreateGenerator();
  }

  mgr->SetMCtruthEventHandler(mcInputHandler);
  
  return mcInputHandler;
}
