#ifndef ALIHBTREADERPPPROD_H
#define ALIHBTREADERPPPROD_H

#include "AliHBTReader.h"
#include "AliHBTReaderTPC.h"

//This reader reads tracks AliTPCtracks.root
//                  particles form tpc_good_tracks 
//I am aware that this file is temporary however we do not have any other PID yet
//Piotr.Skowronski@cern.ch

#include <TString.h>
class TFile;

class AliHBTReaderPPprod: public AliHBTReader
{
  public:
    AliHBTReaderPPprod(const Char_t* trackfilename = "AliTPCtracks.root",
                       const Char_t* clusterfilename = "AliTPCclusters.root", 
                       const Char_t* goodtracksfilename = "good_tracks_tpc",
                       const Char_t* galicefilename = "");

    virtual ~AliHBTReaderPPprod();
    
    Int_t Read(AliHBTRun* particles, AliHBTRun *tracks); //reads tracks and particles and puts them in runs
    
    AliHBTEvent* GetParticleEvent(Int_t); //returns pointer to event with particles 
    AliHBTEvent* GetTrackEvent(Int_t);//returns pointer to event with particles 
    Int_t GetNumberOfPartEvents(); //returns number of particle events
    Int_t GetNumberOfTrackEvents();//returns number of track events
    
  protected:
    //in the future this class is will read global tracking

    
    Int_t OpenFiles(); //opens files to be read
    void CloseFiles(); //close files
    
    AliHBTRun* fParticles; //!simulated particles
    AliHBTRun* fTracks; //!reconstructed tracks (particles)

    TString fTrackFileName; //name of the file with tracks
    TString fClusterFileName;//name of the file with clusters
    TString fGAliceFileName;//name of the file with galice.root
    TString fGoodTPCTracksFileName; //name of text file with good tracks
    
    TFile *fTracksFile; //file with tracks
    TFile *fClustersFile;//file with clusters
    
    Bool_t fIsRead; //flag indicating if the data are already read
    
    
  private:
  public:
    ClassDef(AliHBTReaderPPprod,1)
};

struct GoodTrack //data of good tracks produced by AliTPCComparison.C
 {
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
 };


class AliGoodTracksPP
 { 
   //container for good tracks
   //this class is for internal use only
   
   friend class AliHBTReaderPPprod;
   
   private:
     AliGoodTracksPP(const TString& infilename = TString("good_tracks_tpc")); //constructor
     virtual ~AliGoodTracksPP(); //dctor
   
     const GoodTrack& GetTrack(Int_t event, Int_t n) const; //returns reference to the nth good track in event "event"

     Int_t  fNevents;  //Number of events
     Int_t* fGoodInEvent; //Numbers of good track in event
     struct GoodTrack **fData;
 };


#endif
