void SetDataFromConfigFile(AliAnalysisAlien *plugin, const char* filename, Int_t AnalysisMC, Bool_t esdAna, Int_t nMaxRuns);


//AliAnalysisGrid* CreateAlienHandler(Int_t nFiles, Bool_t AnalysisMC, Int_t runtype, const char* taskname, const char* mode)
AliAnalysisGrid* CreateAlienHandler(Int_t nMaxRuns, Int_t AnalysisMC, Bool_t esdAna, const char* taskname,  const char* nameoutputs, const char* mode, const char* label, const char* alirootver, Int_t task_num)
{
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  plugin->SetRunMode(mode);
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-30-06");
  plugin->SetAliROOTVersion(alirootver);

  Char_t configfile[64];
  sprintf(configfile,"%s.conf",label);

  cout << "Configuration file: " << configfile << endl;

  if(AnalysisMC){
    if(esdAna)
      sprintf(label, "%s_MC_ESDs_%d", label, AnalysisMC);
    else
      sprintf(label, "%s_MC_AODs_%d", label, AnalysisMC);
  }
  else{
    if(esdAna)
      sprintf(label, "%s_Data_ESDs", label);
    else
      sprintf(label, "%s_Data_AODs", label);
  }
  cout << "Label: " << label << endl;

  // output to run numbers
  plugin->SetOutputToRunNo();

  SetDataFromConfigFile(plugin, configfile, AnalysisMC, esdAna, nMaxRuns);  

  //  sprintf(outputfiles,"%s_%s.root %sDebug.root",taskname,label,taskname);
  //Char_t outputfiles[256];
  //sprintf(outputfiles,"%s_Tree.root",taskname);

// Method 2: Declare existing data files (raw collections, xml
// collections, root file) If no path mentioned data is supposed to be
// in the work directory (see SetGridWorkingDir()) XML collections
// added via this method can be combined with the first method if the
// content is compatible (using or not tags)
// plugin->AddDataFile("tag.xml");
// plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  

// Define alien work directory where all files will be copied. Relative to alien $HOME.
  Char_t tmpname[128];
  sprintf(tmpname,"work_%s_%s",taskname,label);
  plugin->SetGridWorkingDir(tmpname);
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  plugin->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include");
  Char_t anasource[256];
  if(task_num!=4){
    sprintf(anasource, "DebugClasses.C AliAnalysisTask%s.cxx", taskname);
  }
  else{
    sprintf(anasource, "DebugClasses.C AliAnalysisTask%s.cxx AliAnalysisTask%sV0.cxx", taskname, taskname);
  }
  plugin->SetAnalysisSource(anasource);
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  Char_t addlib[256];
  if(task_num!=4){
    sprintf(addlib, "DebugClasses.C AliAnalysisTask%s.h AliAnalysisTask%s.cxx", taskname, taskname);
  }else{
    sprintf(addlib, "DebugClasses.C AliAnalysisTask%s.h AliAnalysisTask%s.cxx AliAnalysisTask%sV0.h AliAnalysisTask%sV0.cxx", taskname, taskname, taskname, taskname);
  }
  plugin->SetAdditionalLibs(addlib);
   
   // Declare the output file names separated by blancs.
   // (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetDefaultOutputs(kFALSE);
   plugin->SetOutputFiles(nameoutputs);

   //plugin->SetMergeViaJDL(kTRUE);


   
   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   sprintf(tmpname,"macro_%s_%s.C",taskname,label);
   plugin->SetAnalysisMacro(tmpname);
   
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(100);
   
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
   plugin->SetMaxInitFailed(500);
   
   // Optionally resubmit threshold.
   plugin->SetMasterResubmitThreshold(90);
   
   // Optionally set time to live (default 30000 sec)
   // plugin->SetTTL(7200);
   plugin->SetTTL(30000);
   
   // Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
   
   // Optionally modify the name of the generated JDL (default analysis.jdl)
   sprintf(tmpname,"%s_%s.jdl",taskname,label);
   plugin->SetJDLName(tmpname);
   
   // Optionally modify the executable name (default analysis.sh)
   sprintf(tmpname,"%s_%s.sh",taskname,label);
   plugin->SetExecutable(tmpname);
   
   // Optionally modify job price (default 1)
   plugin->SetPrice(1); 
   // Merge via JDL
   plugin->SetMergeViaJDL(kTRUE);
   // Use fastread option
   plugin->SetFastReadOption(kTRUE);
   // Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   plugin->SetExecutableCommand("aliroot -b -q");

  

   return plugin;
} 





void SetDataFromConfigFile(AliAnalysisAlien *plugin, const char* fileName, Int_t AnalysisMC, Bool_t esdAna, Int_t nMaxRuns)
{


  FILE* file = fopen(fileName,"r");
  if(!file) {
    cout << "File " << fileName << " not found!" << endl;
    return;
  }

  Char_t dummy[128];

  Char_t runperiodpattern[128];
  if(AnalysisMC)
    sprintf(runperiodpattern,"Run period MC%d: %s", AnalysisMC, "%s %s %s");
  else
    sprintf(runperiodpattern,"Run period: %s", "%s %s %s");

  cout << "PATTERN: " << runperiodpattern << endl;

  Int_t nRuns = 0;
  while (fgets(dummy,128,file) != NULL && (nRuns<nMaxRuns||nMaxRuns<=0)) {     
    char runperiod[128], pass[64], aodDir[64];
    Int_t run, a, b;
    if(sscanf(dummy, runperiodpattern, &runperiod, &pass, &aodDir)){
      //if(sscanf(dummy, "Run period: %d %d",&a,&b)){
      Char_t griddatadir[256];
      Char_t datapattern[256];
      
      if(AnalysisMC){
	sprintf(griddatadir, "/alice/sim/%s", runperiod);
	if(esdAna) {
	  sprintf(datapattern,"/*/AliESDs.root");
	} else {
	  sprintf(datapattern,"/%s/*/AliAOD.root", aodDir);
	}
      }
      else{
	plugin->SetRunPrefix("000");
	Int_t year = 0;
	sscanf(runperiod, "LHC%d", &year);
	sprintf(griddatadir, "/alice/data/20%d/%s", year, runperiod);
	if(esdAna) {
	  sprintf(datapattern,"*ESDs/%s/*/AliESDs.root",pass);
	  //sprintf(datapattern,"*ESDs/%s/10000121040024.40/AliESDs.root",pass);
	} else {
	  sprintf(datapattern,"*ESDs/%s/%s/*/AliAOD.root",pass, aodDir);
	}
      }
      cout << "GridDataDir: " << griddatadir << endl;
      cout << "DataPatter: " << datapattern << endl;
      plugin->SetGridDataDir(griddatadir);
      plugin->SetDataPattern(datapattern);
      continue;
    }
    if(sscanf(dummy,"Run: %d %s", &run)){
      plugin->AddRunNumber(run);
      nRuns++;
      continue;
    }
  }



}
