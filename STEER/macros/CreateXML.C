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
//
//  coll_in:    xml collection of tag files 
//           or path to the tag files
//
//  coll_out:   name of the output xml collection
//_________________________________________________________________________
Bool_t CreateXML(char * coll_in="pp.xml", char * coll_out="global2")
{
  TStopwatch timer;
  timer.Start();
  
  
  gSystem->Load("libTreePlayer");    //needed in the case of the string statements
  gSystem->Load("libANALYSIS");      //needed by libANALYSISalice
  gSystem->Load("libANALYSISalice"); //needed AliTagAnalysis

  // Create A tag analysis object and impose some selection criteria
  AliTagAnalysis *TagAna = new AliTagAnalysis(); 

  //Case where the tag files are stored locally
  //TagAna->ChainLocalTags(coll_in);

  //Case where the tag files are stored in the file catalog
  //coll_in (pp.xml) is the xml collection of tag files that was produced 
  //by querying the file catalog:
  // find -x pp /alice/sim/PDC_08/LHC08r/270001/* *tag.root > pp.xml

  TGrid::Connect("alien://pcapiserv01.cern.ch:10000","elopez");
  //TGrid::Connect("alien://"); 
  TAlienCollection* coll = TAlienCollection::Open(coll_in);
  TGridResult* TagResult = coll->GetGridResult("",0,0);
  cout << endl << "Chain Grid Tags..."  << endl;
  TagAna->ChainGridTags(TagResult);

  //__________________________//
  //Usage of string statements//
  //__________________________//
/*
  const char* runTagCuts = "fAliceRunId==270001";
  const char* lhcTagCuts = "fLHCTag.fLHCState==LHC08r";
  const char* detTagCuts = "fDetectorTag.fTPC==1";
  const char* evTagCuts  = "(fEventTag.fTopPtMin >= 1.0)&&(fEventTag.fNumberOfTracks >= 11)&&(fEventTag.fNumberOfTracks <= 12)";
*/

  //________________________________________________//
  //Usage of AliRunTagCuts & AliEventTagCuts classes//
  //________________________________________________//
  // create a RunTagCut object
  AliRunTagCuts *runTagCuts = new AliRunTagCuts();
  // runTagCuts->SetRunId(270001);

  // create an LHCTagCuts object
  AliLHCTagCuts *lhcTagCuts = new AliLHCTagCuts();

  // create an DetectorTagCuts object
  AliDetectorTagCuts *detTagCuts  = new AliDetectorTagCuts();

  // create an EventTagCut object
  AliEventTagCuts *evTagCuts = new AliEventTagCuts();
  //evTagCuts->SetMultiplicityRange(11,12);
  //evTagCuts->SetTopPtMin(1.0);

  // Query the Tags and create the xml collection 
  cout << endl << "Create XML Collection..."  << endl;
  TagAna->CreateXMLCollection(coll_out, runTagCuts, lhcTagCuts, detTagCuts, evTagCuts);

  timer.Stop();
  timer.Print();

  return kTRUE;
}
