#ifndef ALIHBTREADERINTERNAL_H
#define ALIHBTREADERINTERNAL_H

#include "AliHBTReader.h"

//Multi file reader for Internal Data Format
//
//This reader reads data created by
//                  
//Piotr.Skowronski@cern.ch
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html


#include <TString.h>
class TFile;
class TClonesArray;

class AliHBTReaderInternal: public AliHBTReader
{
  public:
    AliHBTReaderInternal();
    AliHBTReaderInternal(const char *filename);
    AliHBTReaderInternal(TObjArray* dirs, const char *filename);
    virtual ~AliHBTReaderInternal();
    
    Int_t Read(AliHBTRun* particles, AliHBTRun *tracks);//reads tracks and particles and puts them in runs
    static Int_t Write(AliHBTReader* reader,const char* outfile = "data.root", Bool_t multcheck = kFALSE);//reads tracks from runs and writes them to file
    
    AliHBTEvent* GetParticleEvent(Int_t);//returns pointer to event with particles
    AliHBTEvent* GetTrackEvent(Int_t);//returns pointer to event with particles 
    Int_t GetNumberOfPartEvents();//returns number of particle events
    Int_t GetNumberOfTrackEvents();//returns number of track events
    
  protected:
    AliHBTRun* fParticles; //!simulated particles
    AliHBTRun* fTracks; //!reconstructed tracks (particles)
    Bool_t     fIsRead;//!flag indicating if the data are already read    
    TString    fFileName;//name of the file with tracks

    Int_t      OpenFile(TFile*& aFile,Int_t event);//opens file to be read for given event 
    static Bool_t FindIndex(TClonesArray* arr,Int_t idx);
    
    ClassDef(AliHBTReaderInternal,1)
};
#endif 
