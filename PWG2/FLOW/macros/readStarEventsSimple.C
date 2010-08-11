//example script on what to do with the star events
//run e.g. like this:
//                    root readStarEventSimple.C

void  readStarEventsSimple()
{
  gSystem->Load("libTree.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libPWG2flowCommon");

  AliStarEventReader*  starReader = new AliStarEventReader( "/data/alice3/jthomas/testData/") ;

  while ( starReader->GetNextEvent() )                                // Get next event
  {
    AliStarEvent* starEvent = starReader->GetEvent();
    if ( !starReader->AcceptEvent(starEvent) ) continue;              // Test if the event is good

    AliFlowEventSimple* flowEvent = new AliFlowEventStar(starEvent);  // make a flow event from a star event (aka "the magic")

    /////analysis here////////////////

    

    //////////////////////////////////

    starEvent->Print("all");
    flowEvent->Print("all");

    delete flowEvent;
    break;
  }

  delete starReader ;
}
