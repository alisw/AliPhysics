#ifndef AliHBTReaderTPC_H
#define AliHBTReaderTPC_H

#include "AliHBTReader.h"

//Multi file reader for TPC
//
//This reader reads tracks AliTPCtracks.root
//                  particles form gAlice
//Piotr.Skowronski@cern.ch
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html

#include <TString.h>
class TFile;

class AliHBTReaderTPC: public AliHBTReader
{
  public:
    AliHBTReaderTPC();
    AliHBTReaderTPC(const Char_t* galicefilename);
    AliHBTReaderTPC(TObjArray* dirs, const Char_t* galicefilename = "galice.root");

    virtual ~AliHBTReaderTPC();
    
    Int_t Read(AliHBTRun* particles, AliHBTRun *tracks);//reads tracks and particles and puts them in runs
    
    AliHBTEvent* GetParticleEvent(Int_t);//returns pointer to event with particles
    AliHBTEvent* GetTrackEvent(Int_t);//returns pointer to event with particles 
    Int_t GetNumberOfPartEvents();//returns number of particle events
    Int_t GetNumberOfTrackEvents();//returns number of track events
    
  protected:
    //in the future this class is will read global tracking
    
    AliHBTRun* fParticles; //!simulated particles
    AliHBTRun* fTracks; //!reconstructed tracks (particles)

    TString fFileName;//name of the file with galice.root

    Bool_t fIsRead;//!flag indicating if the data are already read
  private:
  public:
    ClassDef(AliHBTReaderTPC,2)
};


#endif
