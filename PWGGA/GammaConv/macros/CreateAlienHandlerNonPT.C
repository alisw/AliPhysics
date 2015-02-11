AliAnalysisGrid* CreateAlienHandler( TString mode )
{

   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode			(mode.Data());

   plugin->SetAPIVersion	("V1.1x");
   plugin->SetROOTVersion	("v5-34-08");
   plugin->SetAliROOTVersion	("vAN-20150120");
   //plugin->SetProductionMode	(1);
   plugin->SetAdditionalLibs("PWGGAGammaConv.par");
   plugin->SetAdditionalLibs("libTree.so libVMC.so libGui.so libXMLParser.so libMinuit.so libMinuit2.so libProof.so libGeom.so libPhysics.so libSTEERBase.so libESD.so libAOD.so libOADB.so libANALYSIS.so libANALYSISalice.so libCDB.so libRAWDatabase.so libSTEER.so libCORRFW.so libPWGflowBase.so libPWGflowTasks.so libPWGGAGammaConv.so");
   plugin->SetDefaultOutputs	(kTRUE);
   //plugin->SetOutputFiles("GammaConvV1_4.root GammaConvV1_5.root GammaConvV1_14.root GammaConvV1_15.root");
   plugin->SetOutputToRunNo();
   plugin->SetAnalysisMacro	(Form("GammaConv.C"));
   plugin->SetSplitMaxInputFileNumber(20);
   plugin->SetExecutable	(Form("GammaConv.sh"));
   plugin->SetTTL		(70000);
   plugin->SetInputFormat	("xml-single");
   plugin->SetJDLName		(Form("GammaConv.jdl"));
//   plugin->SetFriendChainName("AliAODGammaConversion.root","libPWGGAGammaConv.so");
   plugin->SetKeepLogs();
   plugin->SetExecutableCommand("root -b -q -x");
   plugin->SetPrice		(1);
   plugin->SetNtestFiles	(5);
   plugin->SetNrunsPerMaster(1);
   plugin->SetMergeDirName("");
   plugin->SetMaxMergeFiles(20);
   plugin->SetMergeViaJDL();
   plugin->SetMaxMergeStages(3);
   plugin->SetSplitMode("se");
   return plugin;
}
