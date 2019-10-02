// Macro to generate macros for the ALICE TRAIN
// Authors: Andrea Gheata, Jan Fiete Grosse-Oetringhaus, Costin Grigoras, Markus Zimmermann
//
//
// requires env variables 

void generate(const char *module = "__ALL__")
{
  Int_t nFiles = TString(gSystem->Getenv("TEST_FILES_NO")).Atoi();//TEST_FILES_NO
  TString dataBasePath(gSystem->Getenv("TEST_DIR"));//TEST_DIR
  TString dataAnchor(gSystem->Getenv("FILE_PATTERN"));//FILE_PATTERN
  Int_t splitMaxInputFileNumber = TString(gSystem->Getenv("SPLIT_MAX_INPUT_FILE_NUMBER")).Atoi();//SPLIT_MAX_INPUT_FILE_NUMBER
  Int_t maxMergeFiles = TString(gSystem->Getenv("MAX_MERGE_FILES")).Atoi();//MAX_MERGE_FILES
  Int_t debugLevel = TString(gSystem->Getenv("DEBUG_LEVEL")).Atoi();//DEBUG_LEVEL
  Int_t ttl = TString(gSystem->Getenv("TTL")).Atoi();//TTL
  TString excludeFiles(gSystem->Getenv("EXCLUDE_FILES"));//EXCLUDE_FILES
  TString friendChainNames(gSystem->Getenv("FRIEND_CHAIN_NAMES"));//FRIEND_CHAIN_NAMES
  TString friendChainLibraries(gSystem->Getenv("FRIEND_CHAIN_LIBRARIES"));//FRIEND_CHAIN_LIBRARIES
  TString additionalpackages(gSystem->Getenv("ADDITIONAL_PACKAGES"));//ADDITIONAL_PACKAGES
  Int_t needsAliEn = TString(gSystem->Getenv("ADDTASK_NEEDS_ALIEN")).Atoi();//ADDTASK_NEEDS_ALIEN
  //TString dataFolder("/home/mazimmer/Test2/dataFolder");//put here a local directory
  TString dataFolder("./");
  Int_t AOD = TString(gSystem->Getenv("AOD")).Atoi();//AOD: 0 ESD, 1 AOD, 2 AOD produced together with ESD, 3 Kinematics only, 4 ESD (special WSDD production) 5 ESD cpass1 (Barrel), 6 ESD cpass1 (Outer)
  Int_t isPP = strcmp(gSystem->Getenv("PP"), "false"); //0 false, 1 true //PP
  TString validOutputFiles = gSystem->Getenv("OUTPUT_FILES");//OUTPUT_FILES
  TString periodName(gSystem->Getenv("PERIOD_NAME"));//PERIOD_NAME

   const char *train_name="lego_train";
   
   Bool_t generateProduction = kFALSE;
   
   if (strcmp(module, "__ALL__") == 0)
     module = "";
   
   if (strcmp(module, "__TRAIN__") == 0)
   {
     generateProduction = kTRUE;
     module = "";
   }

   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   
   if (needsAliEn == 1)
   {
     Printf("Connecting to AliEn...");
     //     TGrid::Connect("alien:");
     TGrid::Connect("alien://",0,0,"t");
   }

   TObjArray *arr = AliAnalysisTaskCfg::ExtractModulesFrom("MLTrainDefinition.cfg");
   
   Printf(">>>>>>> Read train configuration");
   arr->Print();
   
   
   AliAnalysisAlien *plugin = new AliAnalysisAlien(train_name);
   // General plugin settings here
   plugin->SetProductionMode();

   plugin->SetAPIVersion("V1.1x");

   // libraries because we start with root!
   plugin->SetAdditionalRootLibs("libVMC.so libPhysics.so libTree.so libMinuit.so libProof.so libSTEERBase.so libESD.so libAOD.so");

   plugin->SetJobTag("test/test");

   plugin->SetMaxMergeFiles(maxMergeFiles);
   plugin->SetTTL(ttl);
   plugin->SetAnalysisMacro(Form("%s.C", train_name));
   plugin->SetInputFormat("xml-single");
   plugin->SetAnalysisMacro(Form("%s.C", train_name));
   plugin->SetValidationScript("validation.sh");
   
   plugin->SetSplitMaxInputFileNumber(splitMaxInputFileNumber);
   
   plugin->SetRegisterExcludes(excludeFiles + " AliAOD.root");
   if (friendChainNames.Length() > 0)
   {
     if (friendChainLibraries.Length () > 0)
       plugin->SetFriendChainName(friendChainNames, friendChainLibraries);
     else
       plugin->SetFriendChainName(friendChainNames);
   }
   plugin->AddExternalPackage(additionalpackages);
   
//    plugin->SetExecutableCommand("root -b -q"); 
   plugin->SetJDLName(Form("%s.jdl", train_name));
   plugin->SetExecutable(Form("%s.sh", train_name));
   plugin->SetSplitMode("se");
   
   plugin->SetGridOutputDir("./");
   plugin->SetGridWorkingDir("./");

   if (!generateProduction)
    plugin->SetKeepLogs(kTRUE);
   
   plugin->SetMergeViaJDL();
   
   if (dataFolder.Length() == 0)
   {
      Printf("ERROR: TRAIN_TESTDATA not set. Exiting...");
      return;
   }
  
   TString archiveName = "root_archive.zip";
   if (AOD == 2){
     Printf(">>>> USING AOD  file format");
     archiveName = "aod_archive.zip";
   }
   if (friendChainNames.Length() > 0)
     archiveName += ";" + friendChainNames;
    
   TString dataTag;
   dataTag.Form("%s_%s_%s_%d", dataBasePath.Data(), archiveName.Data(), dataAnchor.Data(), nFiles);

   dataTag.ReplaceAll(";", "__");
   dataTag.ReplaceAll("/", "__");
   dataTag.ReplaceAll(".root", "");
   dataTag.ReplaceAll(".zip", "");
    
   TString dataFileName(dataTag + ".txt");
   
   Printf("\n>>>> Test files are taken from: %s", dataFileName.Data());
    
   plugin->SetFileForTestMode(Form("%s/%s", dataFolder.Data(), dataFileName.Data()));
   plugin->SetNtestFiles(nFiles);

   // Is MC only?
   if (AOD == 3)
   {
     Printf(">>>> Expecting MC only production");
     plugin->SetUseMCchain();
   }
   
   // Copy dataset locally
   TString cdir = gSystem->WorkingDirectory();
   gSystem->cd(dataFolder);
      
   // check if files are already copied
   if (gSystem->AccessPathName(dataFileName))
     {
       if (!gGrid)
	 {
	   Printf("Connecting to AliEn...");
	  TGrid::Connect("alien:");
	 }
       
       //	gSystem->Exec(Form("touch %s/../downloading_input_files", cdir.Data()));
       
       
       // special treatment for non-officially produced productions	
       Bool_t specialSet = kFALSE;
       if (periodName == "AMPT_LHC12g6" || periodName == "AMPT_LHC12c3")
	 specialSet = kTRUE;
       
       Bool_t result = kFALSE;
       
       if (specialSet)
	 {
	   result = plugin->CopyLocalDataset(dataBasePath, "Kinematics.root", nFiles, dataFileName, "", dataTag);
	   if (!plugin->CopyLocalDataset(dataBasePath, dataAnchor, nFiles, dataFileName, "", dataTag))
	     result = kFALSE;
	 }
       else
	 result = plugin->CopyLocalDataset(dataBasePath, dataAnchor, nFiles, dataFileName, archiveName, dataTag);
       
       if (!result)
	 {
	   gSystem->Unlink(dataFileName);
	   gSystem->Exec(Form("rm -rf %s", dataTag.Data()));
	   //	  gSystem->Exec(Form("rm %s/../downloading_input_files", cdir.Data()));
	   Printf("ERROR: Could not copy test files. Exiting...");
	   return;
	 }
     }
   else
     {
       // mark files as used (for cleanup)
       gSystem->Exec(Form("touch %s", dataFileName.Data()));
     }
   
   gSystem->cd(cdir);
   
   // Load modules here
   plugin->AddModules(arr);
   plugin->CreateAnalysisManager("train","handlers.C");
   
   AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
   mgr->SetDebugLevel(debugLevel);
   mgr->SetNSysInfo((isPP == 1) ? 2000 : 100);
   
   mgr->SetFileInfoLog("fileinfo.log");
   
   // execute custom configuration
   Int_t error = 0;
   gROOT->Macro("globalvariables.C", &error);
   if (error != 0)
   {
      Printf("ERROR: globalvariables.C was not executed successfully...");
      return;
   }
   
   if (generateProduction)
      plugin->GenerateTrain(train_name);
   else
      plugin->GenerateTest(train_name, module);
   
   // check for illegally defined output files
   validOutputFiles += "," + excludeFiles;
   TString outputFiles = plugin->GetListOfFiles("out");
   tokens = outputFiles.Tokenize(",");
   
   Bool_t valid = kTRUE;
   for (Int_t i=0; i<tokens->GetEntries(); i++)
   {
     if (!validOutputFiles.Contains(tokens->At(i)->GetName()))
     {
       Printf("ERROR: Output file %s requested which is not among the defined ones for this train (%s)", tokens->At(i)->GetName(), validOutputFiles.Data());
       valid = kFALSE;
     }
   }
   delete tokens;
   
   if (!valid)
   {
     Printf(">>>>>>>>> Invalid output files requested. <<<<<<<<<<<<");
     gSystem->Unlink("lego_train.C");
   }
}
