#ifndef ALIHBTREADERITSV2_H
#define ALIHBTREADERITSV2_H

#include "AliHBTReader.h"


class AliLoader;
class AliRunLoader;
class TString;

class AliHBTReaderITSv2: public AliHBTReader
{
  public:
    
    AliHBTReaderITSv2();
    AliHBTReaderITSv2(const Char_t* galicefilename);
    AliHBTReaderITSv2(TObjArray* dirs, const Char_t* galicefilename = "galice.root");

    virtual ~AliHBTReaderITSv2();
    
    void          Rewind();
    
    Bool_t        ReadsTracks() const {return kTRUE;}
    Bool_t        ReadsParticles() const {return kTRUE;}
    
    void          SetMagneticField(Float_t mf){fMagneticField=mf;}
    void          UseMagneticFieldFromRun(Bool_t flag = kTRUE){fUseMagFFromRun=flag;}
    
  protected:
    
    Int_t         ReadNext();//reads tracks and particles and puts them in runs
    Int_t         OpenNextFile();
    void          DoOpenError( const char *va_(fmt), ...);
        
    TString       fFileName;//name of the file with galice.root
    AliRunLoader* fRunLoader;//!Run Loader
    AliLoader*    fITSLoader;//! ITS Loader
        
    Float_t       fMagneticField;//magnetic field value that was enforced while reading
    Bool_t        fUseMagFFromRun;//flag indicating if using field specified in gAlice (kTRUE)
                               // or enforece other defined by fMagneticField
    
    ClassDef(AliHBTReaderITSv2,2)
};

#endif
