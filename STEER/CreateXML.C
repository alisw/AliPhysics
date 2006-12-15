//_________________________________________________________________________
// Macro that creates esd xml collections by querying the tags.
// It addresses the following use cases:
// o) The tag files are stored locally.
//   - One queries the tags by using simple string statements.
//   - One queries the tags by using the corresponding aliroot classes.
// o) The tag files are stored in the file catalog. 
//    In this case the first thing we do is to query the f.c.
//    and extract a collection (xml) of tag files. 
//   - One queries the tags by using simple string statements.
//   - One queries the tags by using the corresponding aliroot classes.
//                                             
// In all cases you create the xml file by using the CreateXMLCollection
// of the AliTagAnalysisClass. The first argument of this method is the
// name of the output xml collection which is stored locally.
//_________________________________________________________________________
Bool_t CreateXML() {
  TStopwatch timer;
  timer.Start();
  
  //needed in the case of the string statements
  gSystem->Load("libTreePlayer.so");

  // Create A tag analysis object and impose some selection criteria
  AliTagAnalysis *TagAna = new AliTagAnalysis(); 

  //Case where the tag files are stored locally
  //TagAna->ChainLocalTags(".");

  //Case where the tag files are stored in the file catalog
  //pp.xml is the xml collection of tag files that was produced 
  //by querying the file catalog.
  TGrid::Connect("alien://pcapiserv01.cern.ch:10000","pchrist");  
  //TGrid::Connect("alien://"); 
  TAlienCollection* coll = TAlienCollection::Open("pp.xml");
  TGridResult* TagResult = coll->GetGridResult("");
  TagAna->ChainGridTags(TagResult);

  //__________________________//
  //Usage of string statements//
  //__________________________//
  /*const char* fRunCuts = "fAliceRunId == 340";
  const char* fEventCuts = "(fEventTag.fTopPtMin >= 1.0)&&(fEventTag.fNumberOfTracks >= 11)&&(fEventTag.fNumberOfTracks <= 12)";
  TagAna->CreateXMLCollection("global",fRunCuts,fEventCuts);*/

  //________________________________________________//
  //Usage of AliRunTagCuts & AliEventTagCuts classes//
  //________________________________________________//
  // create a RunTagCut object
  AliRunTagCuts *RunCuts = new AliRunTagCuts();
  RunCuts->SetRunId(340);
  // create an EventTagCut object
  AliEventTagCuts *EvCuts = new AliEventTagCuts();
  EvCuts->SetMultiplicityRange(11,12);
  EvCuts->SetTopPtMin(1.0);
  TagAna->CreateXMLCollection("global2",RunCuts,EvCuts);

  timer.Stop();
  timer.Print();

  return kTRUE;
}
