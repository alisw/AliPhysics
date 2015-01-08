void  runFlowReaderExample( Int_t maxNumberOfEvents = 1000, const char* inputDataFiles="/Users/snelling/alice_data/jthomas/testData/")
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libPWGflowBase");

  TStopwatch timer;
  timer.Start();

  //define reference particles
  AliStarTrackCuts* rpCuts = AliStarTrackCuts::StandardCuts();
  rpCuts->SetPtMin(0.05);
  rpCuts->SetPtMax(10.);

  //define particles of interest
  AliStarTrackCuts* poiCuts = AliStarTrackCuts::StandardCuts();
  poiCuts->SetPtMin(0.05);
  poiCuts->SetPtMax(10.);

  //define event cuts
  AliStarEventCuts* starEventCuts = AliStarEventCuts::StandardCuts();
  starEventCuts-> SetCentralityIDMax(3);
  starEventCuts-> SetCentralityIDMin(3);

  AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
  AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();   //mcep->SetHarmonic(2); // default is v2

  mcep->Init(); qc->Init();  

  Int_t i=0;
  AliStarEventReader starReader(inputDataFiles) ;
  while ( starReader.GetNextEvent() )                                // Get next event
  {
    AliStarEvent* starEvent = starReader.GetEvent();
    if ( !starEventCuts->PassesCuts(starEvent) ) continue;              // Test if the event is good
    AliFlowEventSimple* flowEvent = new AliFlowEventStar(starEvent,rpCuts,poiCuts);  // make a flow event from a star event 
    // do flow analysis for various methods

    mcep->Make(flowEvent);
    qc->Make(flowEvent);

    delete flowEvent;

    i++;
    cout <<"Event: " << i << "\r"; cout.flush();
    if (i>maxNumberOfEvents) break;
  }

  //---------------------------------------------------------------------------------------  
  // create a new file which will hold the final results of all methods:
  TString outputFileName = "AnalysisResults.root";  
  TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
  // create a new file for each method wich will hold list with final results:
  const Int_t nMethods = 2;
  TString method[nMethods] = {"MCEP","QC"};
  TDirectoryFile *dirFileFinal[nMethods] = {NULL};
  TString fileName[nMethods]; 
  for(Int_t i=0;i<nMethods;i++)
    {
      // form a file name for each method:
      fileName[i]+="output";
      fileName[i]+=method[i].Data();
      fileName[i]+="analysis";
      dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
    } 
  
  // calculating and storing the final results of default flow analysis:
  mcep->Finish();    mcep->WriteHistograms(dirFileFinal[0]);
  qc->Finish();      qc->WriteHistograms(dirFileFinal[1]);
  //---------------------------------------------------------------------------------------  
  
  outputFile->Close();
  delete outputFile;
  
  timer.Stop();
  cout << endl;
  timer.Print();
  
  delete rpCuts;
  delete poiCuts;
  delete starEventCuts;
}
