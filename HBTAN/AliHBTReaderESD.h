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
class AliRunLoader;

class AliHBTReaderESD: public AliHBTReader
{
  public:
    AliHBTReaderESD(const Char_t* esdfilename = "AliESDs.root", const Char_t* galfilename = "galice.root");

    AliHBTReaderESD(TObjArray* dirs,const Char_t* esdfilename = "AliESDs.root", const Char_t* galfilename = "galice.root");

    virtual ~AliHBTReaderESD();
    
    void          Rewind();
    
    void          ReadParticles(Bool_t flag){fReadParticles = flag;}
    Bool_t        ReadsTracks() const {return kTRUE;}
    Bool_t        ReadsParticles() const {return fReadParticles;}
    
    enum ESpecies {kESDElectron = 0, kESDMuon, kESDPion, kESDKaon, kESDProton, kNSpecies};
    static Int_t  GetSpeciesPdgCode(ESpecies spec);//skowron
    
  protected:
    Int_t         ReadNext();
    TFile*        OpenFile(Int_t evno);//opens files to be read for given event
    void          CloseFiles(TFile*);//close files

    TString       fESDFileName;//name of the file with tracks
    TString       fGAlFileName;//name of the file with tracks
    TFile*        fFile;//! pointer to current ESD file
    AliRunLoader* fRunLoader;//!Run Loader
    Bool_t        fReadParticles;//flag indicating wether to read particles from kinematics
    
  private:
    ClassDef(AliHBTReaderESD,2)
};


#endif
