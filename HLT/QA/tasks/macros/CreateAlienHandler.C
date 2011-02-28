// $Id$
/*
 * Alien handler macro to run an HLT analysis task for comparing the offline
 * with the HLT esd tree on the GRID.
 *
 * It is loaded in compare_HLT_offline_grid.C.
 *
 * @ingroup alihlt_qa
 * @author Kalliopi.Kanaki@ift.uib.no
 */

AliAnalysisGrid* CreateAlienHandler(TString runNumber, TString dataDir, TString gridWorkingDir, TString gridOutputDir, const char* mode = "full", const char* detectorTask="global"){
  
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  // source /tmp/gclient_env_$UID in the current shell.
  
  //if(!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(mode);
  
  // check the versions available on alien with the command 'packages'
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-27-06d");
  plugin->SetAliROOTVersion("v4-21-16-AN");
  
  cout << "===========================================================================================" << endl;
  cout << "  " << endl;
  cout << " The user is responsible for updating the AliRoot and ROOT versions for running on the GRID."<< endl;  
  cout << "  " << endl; 
  cout << "===========================================================================================" << endl;

  //Allow non-default outputs
  plugin->SetDefaultOutputs(kFALSE);

  // data alien directory
  plugin->SetGridDataDir(dataDir.Data());
  
  // Set data search pattern
  plugin->SetDataPattern("*pass1/*ESDs.root");
  //plugin->SetDataPattern("*/*ESDs.root");
    
  plugin->AddRunNumber(runNumber); 
  //plugin->SetRunRange(xxx,yyy);
 
  // define working and output directories
  plugin->SetGridWorkingDir(gridWorkingDir); // relative to $HOME
  plugin->SetGridOutputDir(gridOutputDir);   // relative to working dir
  plugin->SetOverwriteMode();                // overwrites the contents of the working and output directory
  
  Bool_t bTPC=kFALSE, bPHOS=kFALSE, bEMCAL=kFALSE, bITS=kFALSE, bGLOBAL=kFALSE, bD0=kFALSE, bCB=kFALSE;
 
  TString allArgs = detectorTask;
  TString argument;
 
  TObjArray *pTokens = allArgs.Tokenize(" ");
  if(pTokens){
     for(int i=0; i<pTokens->GetEntries(); i++){
         argument=((TObjString*)pTokens->At(i))->GetString();
         if(argument.IsNull()) continue;

         if(argument.CompareTo("tpc", TString::kIgnoreCase)==0){
	    bTPC = kTRUE;
	    continue;
         }        
         if(argument.CompareTo("phos", TString::kIgnoreCase)==0){
  	    bPHOS = kTRUE;
	    continue;
         }         
         if(argument.CompareTo("emcal", TString::kIgnoreCase)==0){
  	    bEMCAL = kTRUE;
	    continue;
         }         
	 if(argument.CompareTo("its", TString::kIgnoreCase)==0){
  	    bITS = kTRUE;
	    continue;
         }	
	 if(argument.CompareTo("global", TString::kIgnoreCase)==0){
  	    bGLOBAL = kTRUE;
	    continue;
         }        
	 if(argument.CompareTo("D0", TString::kIgnoreCase)==0){
  	    bD0 = kTRUE;
	    continue;
         }  
	 if(argument.CompareTo("cb", TString::kIgnoreCase)==0){
  	    bCB = kTRUE;
	    continue;
         }  
	 if(argument.CompareTo("all",TString::kIgnoreCase)==0){
	    bTPC    = kTRUE;
	    bPHOS   = kTRUE;
	    bEMCAL  = kTRUE;
	    bITS    = kTRUE;
	    bGLOBAL = kTRUE;
	    bD0     = kTRUE;
	    bCB     = kTRUE;
	    continue;
         }
         else break;
    }
  }
    
  // libs, seperate by blanks
  if(bTPC){  
    plugin->SetAnalysisSource("AliAnalysisTaskHLTTPC.cxx");  
    plugin->SetAdditionalLibs("AliAnalysisTaskHLTTPC.h AliAnalysisTaskHLTTPC.cxx");
    plugin->SetOutputFiles("HLT-OFFLINE-TPC-comparison.root");    
  }
  if(bITS){  
    plugin->SetAnalysisSource("AliAnalysisTaskHLTITS.cxx");  
    plugin->SetAdditionalLibs("AliAnalysisTaskHLTITS.h AliAnalysisTaskHLTITS.cxx");
    plugin->SetOutputFiles("HLT-OFFLINE-ITS-comparison.root");    
  }
  if(bPHOS && bEMCAL) {
    plugin->AddIncludePath("-I$ROOTSYS -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/RAW -I$ALICE_ROOT/STEER -I$ALICE_ROOT/PHOS -I$ALICE_ROOT/HLT/BASE -I$ALICE_ROOT/HLT/BASE/util -I$ALICE_ROOT/HLT/global/physics");
    plugin->SetAnalysisSource("AliAnalysisTaskHLTCalo.cxx AliAnalysisTaskHLTPHOS.cxx AliAnalysisTaskHLTEMCAL.cxx");  
    plugin->SetAdditionalLibs("libRAWDatabase.so libProof.so libGui.so libCDB.so libSTEER.so libHLTbase.so libAliHLTUtil.so libAliHLTGlobal.so AliAnalysisTaskHLTCalo.cxx AliAnalysisTaskHLTCalo.h AliAnalysisTaskHLTPHOS.cxx AliAnalysisTaskHLTPHOS.h AliAnalysisTaskHLTEMCAL.cxx AliAnalysisTaskHLTEMCAL.h");  
    plugin->SetOutputFiles("HLT-OFFLINE-PHOS-comparison.root HLT-OFFLINE-EMCAL-comparison.root");    

  } else if(bPHOS){  
    plugin->AddIncludePath("-I$ROOTSYS -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/RAW -I$ALICE_ROOT/STEER -I$ALICE_ROOT/PHOS -I$ALICE_ROOT/HLT/BASE -I$ALICE_ROOT/HLT/BASE/util -I$ALICE_ROOT/HLT/global/physics");
    plugin->SetAnalysisSource("AliAnalysisTaskHLTCalo.cxx AliAnalysisTaskHLTPHOS.cxx");  
    plugin->SetAdditionalLibs("libRAWDatabase.so libProof.so libGui.so libCDB.so libSTEER.so libHLTbase.so libAliHLTUtil.so libAliHLTGlobal.so AliAnalysisTaskHLTCalo.cxx AliAnalysisTaskHLTCalo.h AliAnalysisTaskHLTPHOS.cxx AliAnalysisTaskHLTPHOS.h");  
    plugin->SetOutputFiles("HLT-OFFLINE-PHOS-comparison.root");    
  } else if(bEMCAL){  
    plugin->AddIncludePath("-I$ROOTSYS -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/RAW -I$ALICE_ROOT/STEER -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/HLT/BASE -I$ALICE_ROOT/HLT/BASE/util -I$ALICE_ROOT/HLT/global/physics");
    plugin->SetAnalysisSource("AliAnalysisTaskHLTCalo.cxx AliAnalysisTaskHLTEMCAL.cxx");  
    plugin->SetAdditionalLibs("libRAWDatabase.so libProof.so libGui.so libCDB.so libSTEER.so libHLTbase.so libAliHLTUtil.so libAliHLTGlobal.so AliAnalysisTaskHLTCalo.cxx AliAnalysisTaskHLTCalo.h AliAnalysisTaskHLTEMCAL.cxx AliAnalysisTaskHLTEMCAL.h");  
    plugin->SetOutputFiles("HLT-OFFLINE-EMCAL-comparison.root");    
  }
  if(bGLOBAL){  
    plugin->AddIncludePath("-I$ALICE_ROOT/HLT/BASE");
    plugin->SetAnalysisSource("AliAnalysisTaskHLT.cxx");  
    plugin->SetAdditionalLibs("libHLTbase.so AliAnalysisTaskHLT.h AliAnalysisTaskHLT.cxx"); 
    plugin->SetOutputFiles("HLT-OFFLINE-GLOBAL-comparison.root");
  }
  if(bD0){
    //plugin->AddIncludePath("-I$ROOTSYS -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/RAW -I$ALICE_ROOT/STEER -I$ALICE_ROOT/HLT/BASE -I$ALICE_ROOT/HLT/BASE/util -I$ALICE_ROOT/HLT/global/physics -I$ALICE_ROOT/HLT/trigger");
    //plugin->SetAdditionalLibs("libRAWDatabase.so libProof.so libGui.so libCDB.so libSTEER.so libHLTbase.so libAliHLTUtil.so libAliHLTGlobal.so AliAnalysisTaskD0Trigger.cxx AliAnalysisTaskD0Trigger.h");  
    plugin->SetAnalysisSource("AliAnalysisTaskD0Trigger.cxx");  
    plugin->SetAdditionalLibs("AliAnalysisTaskD0Trigger.h AliAnalysisTaskD0Trigger.cxx"); 
    plugin->SetOutputFiles("HLT-OFFLINE-D0-comparison.root");    
  }
  if(bCB){
    plugin->AddIncludePath("-I$ALICE_ROOT/HLT/BASE");
    plugin->SetAnalysisSource("AliAnalysisTaskHLTCentralBarrel.cxx");  
    plugin->SetAdditionalLibs("AliAnalysisTaskHLTCentralBarrel.h AliAnalysisTaskHLTCentralBarrel.cxx"); 
    plugin->SetOutputFiles("HLT-OFFLINE-CentralBarrel-comparison.root");    
  }
  
  // Optionally define the files to be archived.
  plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("runComparison.C");
  plugin->SetExecutable("comparison.sh");

  plugin->SetSplitMaxInputFileNumber(100);
  
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMaxInitFailed(10);
  
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90); // in %

  plugin->SetTTL(30000);// in sec
  
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
 
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("analysis.jdl");
 
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  
  return plugin;
}
