#ifndef ALIHBTREADERITSV2_H
#define ALIHBTREADERITSV2_H

#include "AliHBTReader.h"

#include <TString.h>

class TFile;

class AliHBTReaderITSv2: public AliHBTReader
{
  public:
    
    AliHBTReaderITSv2();
    AliHBTReaderITSv2(const Char_t* galicefilename);
    AliHBTReaderITSv2(TObjArray* dirs, const Char_t* galicefilename = "galice.root");

    virtual ~AliHBTReaderITSv2();
    
    Int_t Read(AliHBTRun*, AliHBTRun*);//reads tracks and particles and puts them in runs
    
    AliHBTEvent* GetParticleEvent(Int_t);//returns pointer to event with particles
    AliHBTEvent* GetTrackEvent(Int_t);//returns pointer to event with particles 
    Int_t GetNumberOfPartEvents();//returns number of particle events
    Int_t GetNumberOfTrackEvents();//returns number of track events

    void SetMagneticField(Float_t mf){fMagneticField=mf;}
    void UseMagneticFieldFromRun(Bool_t flag = kTRUE){fUseMagFFromRun=flag;}
    
  protected:
    
    AliHBTRun* fParticles; //!simulated particles
    AliHBTRun* fTracks; //!reconstructed tracks (particles)
    
    TString fFileName;//name of the file with galice.root
  
    Bool_t fIsRead;//!flag indicating if the data are already read

    Float_t    fMagneticField;//magnetic field value that was enforced while reading
    Bool_t     fUseMagFFromRun;//flag indicating if using field specified in gAlice (kTRUE)
                               // or enforece other defined by fMagneticField
    ClassDef(AliHBTReaderITSv2,1)
};

#endif
