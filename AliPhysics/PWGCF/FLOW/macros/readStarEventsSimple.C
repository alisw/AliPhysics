//example script on what to do with the star events
//run e.g. like this:
//                    root readStarEventSimple.C

void  readStarEventsSimple()
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libPWGflowBase");

  Int_t maxNumberOfEvents = 1000;

  //define reference particles
  AliStarTrackCuts* rpCuts = AliStarTrackCuts::StandardCuts();

  //define particles of interest
  AliStarTrackCuts* poiCuts = AliStarTrackCuts::StandardCuts();
  poiCuts->SetPtMin(1.0);

  //define event cuts
  AliStarEventCuts* starEventCuts = AliStarEventCuts::StandardCuts();

  Int_t i=0;
  AliStarEventReader starReader("/data/alice3/jthomas/testData/") ;
  while ( starReader.GetNextEvent() )                                // Get next event
  {
    AliStarEvent* starEvent = starReader.GetEvent();
    if ( !starEventCuts->PassesCuts(starEvent) ) continue;              // Test if the event is good

    AliFlowEventSimple* flowEvent = new AliFlowEventStar(starEvent,rpCuts,poiCuts);  // make a flow event from a star event (aka "the magic")

    /////analysis here////////////////

    

    //////////////////////////////////

    //starEvent->Print("all");
    flowEvent->Print();

    delete flowEvent;

    i++;
    if (i>maxNumberOfEvents) break;
  }
  delete rpCuts;
  delete poiCuts;
  delete starEventCuts;
}
