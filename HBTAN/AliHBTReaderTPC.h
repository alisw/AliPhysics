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
    AliHBTReaderTPC(const Char_t* trackfilename = "AliTPCtracks.root",
                      const Char_t* clusterfilename = "AliTPCclusters.root",
	  const Char_t* galicefilename = "galice.root");

    AliHBTReaderTPC(TObjArray* dirs,
                      const Char_t* trackfilename = "AliTPCtracks.root",
                      const Char_t* clusterfilename = "AliTPCclusters.root",
	  const Char_t* galicefilename = "galice.root");

    virtual ~AliHBTReaderTPC();
    
    Int_t Read(AliHBTRun* particles, AliHBTRun *tracks);//reads tracks and particles and puts them in runs
    
    AliHBTEvent* GetParticleEvent(Int_t);//returns pointer to event with particles
    AliHBTEvent* GetTrackEvent(Int_t);//returns pointer to event with particles 
    Int_t GetNumberOfPartEvents();//returns number of particle events
    Int_t GetNumberOfTrackEvents();//returns number of track events
    
    void SetMagneticField(Float_t mf){fMagneticField=mf;}
    void UseMagneticFieldFromRun(Bool_t flag = kTRUE){fUseMagFFromRun=flag;}
  protected:
    //in the future this class is will read global tracking

    
    Int_t OpenFiles(TFile*&,TFile*&,TFile*&,Int_t);//opens files to be read for given event
    void CloseFiles(TFile*&,TFile*&,TFile*&);//close files

    
    
    AliHBTRun* fParticles; //!simulated particles
    AliHBTRun* fTracks; //!reconstructed tracks (particles)
    

    TString    fTrackFileName;//name of the file with tracks
    TString    fClusterFileName;//name of the file with clusters
    TString    fGAliceFileName;//name of the file with galice.root

    Bool_t     fIsRead;//!flag indicating if the data are already read
    
    Float_t    fMagneticField;//magnetic field value that was enforced while reading
    Bool_t     fUseMagFFromRun;//flag indicating if using field specified in gAlice (kTRUE)
                               // or enforece other defined by fMagneticField
  private:
    ClassDef(AliHBTReaderTPC,2)
};


#endif
