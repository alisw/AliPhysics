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
class AliRunLoader;
class AliTPCLoader;

class AliHBTReaderTPC: public AliHBTReader
{
  public:
    AliHBTReaderTPC();
    AliHBTReaderTPC(const Char_t* galicefilename);
    AliHBTReaderTPC(TObjArray* dirs, const Char_t* galicefilename = "galice.root");

    virtual ~AliHBTReaderTPC();
    
    void          Rewind();
    
    Bool_t        ReadsTracks() const {return kTRUE;}
    Bool_t        ReadsParticles() const {return kTRUE;}
    
    void          SetMagneticField(Float_t mf){fMagneticField=mf;}
    void          UseMagneticFieldFromRun(Bool_t flag = kTRUE){fUseMagFFromRun=flag;}
    
  protected:
    //in the future this class is will read global tracking
    Int_t         ReadNext();
    Int_t         OpenNextSession();
    void          DoOpenError(const char* msgfmt, ...);
    
    TString       fFileName;//name of the file with galice.root
    AliRunLoader* fRunLoader;//!RL
    AliTPCLoader* fTPCLoader;//!TPCLoader
    Float_t       fMagneticField;//magnetic field value that was enforced while reading
    Bool_t        fUseMagFFromRun;//flag indicating if using field specified in gAlice (kTRUE)
                               // or enforece other defined by fMagneticField

  private:
  public:
    ClassDef(AliHBTReaderTPC,3)
};


#endif
