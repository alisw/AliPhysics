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

AliAnalysisGrid *CreateAlienHandler(TString runNumber, TString dataDir, TString gridWorkingDir, TString gridOutputDir, const char* mode = "full", const char* detectorTask="global"){
  
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  // source /tmp/gclient_env_$UID in the current shell.
  
  //if(!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(mode);
  
  // check the versions available on alien with the command 'packages'
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-26-00b-4");
  plugin->SetAliROOTVersion("v4-19-11-AN");

  // data alien directory
  //plugin->SetGridDataDir("/alice/data/2010/LHC10b");
  plugin->SetGridDataDir(dataDir.Data());
  
  // Set data search pattern
  plugin->SetDataPattern("*ESDs.root");
  //plugin->SetDataPattern("pass1/*tags.root");
  
  //plugin->AddRunNumber(115322); 
  plugin->AddRunNumber(runNumber); 
  //plugin->SetRunRange(xxx,yyy);
 
  // define working and output directories
  //plugin->SetGridWorkingDir("ESDcomparison"); // relative to $HOME
  //plugin->SetGridOutputDir("output");         // relative to working dir
  plugin->SetGridWorkingDir(gridWorkingDir); // relative to $HOME
  plugin->SetGridOutputDir(gridOutputDir);   // relative to working dir
  
  
  Bool_t bTPC=kFALSE, bPHOS=kFALSE, bITS=kFALSE, bGLOBAL=kFALSE;
 
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
	 if(argument.CompareTo("its", TString::kIgnoreCase)==0){
  	    bITS = kTRUE;
	    continue;
         }	
	 if(argument.CompareTo("global", TString::kIgnoreCase)==0){
  	    bGLOBAL = kTRUE;
	    continue;
         }        
	 if(argument.CompareTo("all",TString::kIgnoreCase)==0){
	    bTPC    = kTRUE;
	    bPHOS   = kTRUE;
	    bITS    = kTRUE;
	    bGLOBAL = kTRUE;
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
  if(bPHOS){  
    plugin->SetAnalysisSource("AliAnalysisTaskHLTPHOS.cxx");  
    plugin->SetOutputFiles("HLT-OFFLINE-GLOBAL-comparison.root");
    plugin->SetOutputFiles("HLT-OFFLINE-PHOS-comparison.root");    
  }
  if(bGLOBAL){  
    plugin->SetAnalysisSource("AliAnalysisTaskHLT.cxx");  
    plugin->SetAdditionalLibs("AliAnalysisTaskHLT.h AliAnalysisTaskHLT.cxx"); 
    plugin->SetOutputFiles("HLT-OFFLINE-GLOBAL-comparison.root");
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
