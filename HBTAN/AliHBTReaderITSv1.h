#ifndef ALIHBTREADERITSV1_H
#define ALIHBTREADERITSV1_H

#include "AliHBTReader.h"

#include <TString.h>


class AliHBTReaderITSv1: public AliHBTReader
{
  public:    
    AliHBTReaderITS(const Char_t* goodtracksfilename = "itsgood_tracks");
    virtual ~AliHBTReaderITS();
    
    Int_t Read(AliHBTRun* particles, AliHBTRun *tracks);//reads tracks and particles and puts them in runs
    
    AliHBTEvent* GetParticleEvent(Int_t);//returns pointer to event with particles
    AliHBTEvent* GetTrackEvent(Int_t);//returns pointer to event with particles
    Int_t GetNumberOfPartEvents();//returns number of particle events
    Int_t GetNumberOfTrackEvents();//returns number of track events
    
  protected:
    TString fGoodITSTracksFileName;
    
    AliHBTRun* fParticles; //!simulated particles
    AliHBTRun* fTracks; //!reconstructed tracks (particles)
    
    Bool_t fIsRead;//flag indicating if the data are already read
    
  private:
  public:
    ClassDef(AliHBTReaderITS,1)
};

struct GoodTrackITSv1 //good tracks produced by ITSComparison V1
{
  Int_t fEventN; //event number
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
  Float_t pxg,pyg,pzg,ptg;
  Bool_t flag;
};

class AliGoodTracksITSv1
 { 
   //container for good tracks ITS tracking V1
   //this class is for internal use only
   friend class AliHBTReaderITSv1;
   
   private:
     AliGoodTracksITSv1(const TString& infilename = TString("itsgood_tracks"));
     ~AliGoodTracksITSv1();
   
     const GoodTrackITSv1& GetTrack(Int_t event, Int_t n) const;

     Int_t  fNevents;  //Number of events
     Int_t* fGoodInEvent; //Numbers of good track in event
     struct GoodTrack **fData;
 };

