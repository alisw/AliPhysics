#ifndef ALIHBTREADERITSV1_H
#define ALIHBTREADERITSV1_H
//____________________________________________________________________
//////////////////////////////////////////////////////////////////////
//                                                                  //
//  class AliHBTReaderITSv1                                         //
//                                                                  //
//  Reader for ITSv1 tracks. Not maintained since v1 is not         //
//  supposed to be used                                             //
//                                                                  //
//////////////////////////////////////////////////////////////////////


#include "AliHBTReader.h"

#include <TString.h>

class TObjArray;
class TFile;
class AliHBTReaderITSv1: public AliHBTReader
{
  public:    
    AliHBTReaderITSv1(const Char_t* tracksfilename="itstracks.root",
                      const Char_t* galicefilename="galice.root");
    AliHBTReaderITSv1(TObjArray* dirs,
                      const Char_t* tracksfilename="itstracks.root",
                      const Char_t* galicefilename="galice.root");    
    
    
    virtual ~AliHBTReaderITSv1();
    
    Int_t Read(AliHBTRun* particles, AliHBTRun *tracks);//reads tracks and particles and puts them in runs
    
    AliHBTEvent* GetParticleEvent(Int_t);//returns pointer to event with particles
    AliHBTEvent* GetTrackEvent(Int_t);//returns pointer to event with particles
    Int_t GetNumberOfPartEvents();//returns number of particle events
    Int_t GetNumberOfTrackEvents();//returns number of track events
    
  protected:
    TString fITSTracksFileName; //name of the file with tracks
    TString fGAliceFileName;//name of the file containing Run Loader
    
    AliHBTRun* fParticles; //!simulated particles
    AliHBTRun* fTracks; //!reconstructed tracks (particles)
    
    Bool_t fIsRead;//flag indicating if the data are already read
   
    TFile* OpenTrackFile(Int_t);//opens files to be read for given directoru nomber in fDirs Array
    TFile* OpenGAliceFile(Int_t);
    
  private:
    ClassDef(AliHBTReaderITSv1,1)
};

#endif
