
AliAnalysisGrid* CreateAlienHandler(TString dataset, const char * runList, const char * suffix, TList * listToLoad, const char * mode = "full", Bool_t isMC = 0)
{
  
  TGrid::Connect("alien:"); // Connecto to alien

  AliAnalysisAlien * handler = new AliAnalysisAlien("test");
  // input
  handler->SetGridDataDir(dataset);
  handler->AddRunNumber(runList);
  if(dataset.Contains("sim")) handler->SetDataPattern("*AliESDs.root");
  else handler->SetDataPattern("*pass2/*AliESDs.root");
  handler->SetAnalysisMacro("MultPb.C");
  handler->SetJDLName("MultPb.jdl");
  handler->SetAdditionalRootLibs("libCore libTree libGeom libVMC libPhysics libSTEERBase libESD libAOD libANALYSIS libOADB libANALYSISalice");   
  handler->SetRunMode(mode);
  TIterator * iter = listToLoad->MakeIterator();
  TObjString * name = 0;
  TString sources = "";
  TString sourcescxxOnly = "";
  while (name = (TObjString *)iter->Next()) {
     gSystem->ExpandPathName(name->String());
     name->String().ReplaceAll("+","");     
     sources = sources + gSystem->BaseName(name->String().Data()) + " ";
     TString header = gSystem->BaseName(name->String().Data());
     header.ReplaceAll("cxx","h");
     sources = sources + header.Data() + " ";
     if(name->String().Contains("cxx")) sourcescxxOnly = sourcescxxOnly + gSystem->BaseName(name->String().Data()) + " ";
   }
  handler->SetAnalysisSource(sourcescxxOnly.Data());
  handler->SetAdditionalLibs(sources.Data());   

  // Alirootversion
  //VO_ALICE@GEANT3::v1-12  
  handler->SetAPIVersion("V1.1x");
  handler->SetROOTVersion   ("v5-30-03-1")    ;
  handler->SetAliROOTVersion("v5-02-12-AN") ;

  // output
  TString runs = runList;
  runs.ReplaceAll(" ","_");
  handler->SetGridWorkingDir(TString("MultPbPb_") + gSystem->BaseName(dataset)+suffix+"_"+runs);
  handler->SetGridOutputDir("out");
  //  handler->SetOutputFiles("sec.root");
  handler->SetMergeExcludes("sec.root EventStat_temp.root event_stat.root multPbPbtracks.root"); // don't merge files; FIXME: if processing many bin, the name of files not to be merged should be changed

  return handler;

}
