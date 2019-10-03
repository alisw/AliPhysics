// Title:   Sample Analysis Macro for looking at data on STAR NTuples
// Author:  Jim Thomas     jhthomas@lbl.gov
// Modified: Mikolaj Krzewicki  mikolaj.krzewicki@cern.ch
// Date:    04-Aug-2010
//

#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLeaf.h"

void  readStarEvents()
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libPWGflowBase");

  Long64_t EventCounter = 0 ;
  Int_t    TrackCounter = 0 ;
  Int_t    PrintInfo    = 1 ;                               // Print event and track data (1/0) (on/off)
  Int_t    PDG_ID           ;                               // for PID test using Particle Data Group PID names

  TH1F* myHistogram = new TH1F("Pt","Transverse Momentum", 100, 0.0, 3.0) ;
  TCanvas* myCanvas = new TCanvas("c1","c1",150,50,500,500) ;
  myCanvas -> cd() ;
  myHistogram -> Draw() ;                // Prepare the histogram and canvas for updates

  AliStarEventReader*  starReader = new AliStarEventReader( "/data/alice3/jthomas/testData/") ;
  AliStarEvent* starEvent = starReader->GetEvent();

  while ( starReader->GetNextEvent() )                                // Get next event
  {
    if ( !starReader->AcceptEvent(starEvent) ) continue;                           // Test if the event is good

    if ( PrintInfo == 1 ) starEvent->Print() ;  // Print basic information for this event

    if ((EventCounter%100) == 0)
    {
      TrackCounter = 0 ;

      //loop over track
      for ( Int_t j = 0 ; j < starEvent->GetNumberOfTracks() ; j++ )
      {
        AliStarTrack* starTrack = starEvent->GetTrack(j);             // Get next track
        if ( !starReader->AcceptTrack(starTrack) ) continue;               // Test if the track is good

        // Do something useful with the track
        if ( TrackCounter<5 && PrintInfo==1 ) starTrack->Print(); ; // Print the track data

        PDG_ID = starReader->ParticleID(starTrack);         // Assign Particle Data Book ID number to track (very simple PID algorithm)
        // Note this is an overly simplified ID algorithm and should really be extended for serious work
        if ( PDG_ID == 211 ) myHistogram->Fill( starTrack->GetPt() ) ;
        TrackCounter ++       ;               // Count the number of accepted and analyzed tracks
      }
      
      //starEvent->Print("all"); //if you want to print info for all tracks
      
      
      if ( PrintInfo == 1 )                         // Talk to the user and decide whether to continue, or stop.
      {
        cout << endl <<  "Enter 0 to quit, 1 to continue, 2 to continue without printing " << endl ;
        Int_t k = 0 ;
        cin >> k ;
        if ( k == 2 ) PrintInfo = 0 ;
        if (k <= 0) return ;
      }
    }
    EventCounter ++ ;
    if ( (EventCounter%10000) == 0 )   cout << EventCounter << endl ;
    if ( (EventCounter%10000) == 0 )   myCanvas->Update() ;
  }

  myCanvas->Update() ;                                      // Final update of histograms
  delete starReader ;
  starReader = NULL ;                              // Prepare to exit
}
