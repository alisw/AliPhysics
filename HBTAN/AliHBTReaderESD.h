#ifndef ALIHBTREADERESD_H
#define ALIHBTREADERESD_H

#include "AliHBTReader.h"
//___________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Multi file reader for ESD                                               //
//                                                                         //
// This reader reads tracks from Event Summary Data                        //
// do not read particles                                                   //
// Piotr.Skowronski@cern.ch                                                //
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include <TString.h>
class TFile;

class AliHBTReaderESD: public AliHBTReader
{
  public:
    AliHBTReaderESD(const Char_t* esdfilename = "AliESDs.root");

    AliHBTReaderESD(TObjArray* dirs,const Char_t* esdfilename = "AliESDs.root");

    virtual ~AliHBTReaderESD();
    
    Int_t Read(AliHBTRun* particles, AliHBTRun *tracks);//reads tracks and particles and puts them in runs
    
    AliHBTEvent* GetParticleEvent(Int_t);//returns pointer to event with particles
    AliHBTEvent* GetTrackEvent(Int_t);//returns pointer to event with particles 
    Int_t        GetNumberOfPartEvents();//returns number of particle events
    Int_t        GetNumberOfTrackEvents();//returns number of track events
    
    enum ESpecies {kESDElectron = 0, kESDMuon, kESDPion, kESDKaon, kESDProton, kNSpecies};
    static Int_t GetSpeciesPdgCode(ESpecies spec);//skowron
  protected:
      
    TFile*       OpenFile(Int_t evno);//opens files to be read for given event
    void         CloseFiles(TFile*);//close files
    
    AliHBTRun*   fParticles; //!simulated particles
    AliHBTRun*   fTracks; //!reconstructed tracks (particles)

    TString      fESDFileName;//name of the file with tracks
    Bool_t       fIsRead;//!flag indicating if the data are already read
    
  private:
    ClassDef(AliHBTReaderESD,1)
};


#endif
