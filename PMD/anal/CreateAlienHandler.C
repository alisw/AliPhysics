//=========================================================================//
//                                                                         //
//             A template of CreateAlienHandler for PMD analysis           //
//            You can copy it and add features to ir as your need          //
//                      There are many way to do it!!!                     //
//                               Satyajit Jena                             //
//                               sjena@cern.ch                             //
//                                13/04/2012                               //
//                                                                         //
//=========================================================================//


AliAnalysisGrid* CreateAlienHandler(const char *gridmode) {
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(gridmode);
  plugin->SetNtestFiles(1);
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-08-6");
  plugin->SetAliROOTVersion("vAN-20141207");
  
  // for simulation 
  //plugin->SetGridDataDir(" /alice/sim/LHC11a10a_bis");
  // plugin->SetDataPattern("*ESDs.root");
  
  plugin->SetGridDataDir("/alice/data/2010/LHC10h");
  plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");
  plugin->SetRunPrefix("000"); //  for data only IMPORTANT!
    
  //------ Add Run Numbers
  plugin->AddRunNumber(137161);  
  plugin->AddRunNumber(137231);  

  // Where to Store Output 
  //---------------------------------------
  plugin->SetGridWorkingDir("TestRun/pbpb");
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  plugin->SetOutputToRunNo();
  plugin->SetAnalysisSource("AliPMDAnalysisTaskPbPb.cxx");
  plugin->SetAdditionalLibs("AliPMDAnalysisTaskPbPb.h AliPMDAnalysisTaskPbPb.cxx");
  plugin->SetDefaultOutputs(kTRUE);
  plugin->SetAnalysisMacro("pmdPbPbTask.C");
  plugin->SetExecutable("pmdPbPbTask.sh");
  plugin->SetJDLName("pmdPbPbTask.jdl"); 

  //  plugin->SetOverwriteMode();
  
  // No need to change it 
  //___________________________________________
  plugin->SetSplitMaxInputFileNumber(100);
  plugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");
  plugin->SetMergeViaJDL(kTRUE);
  
  return plugin;
}
