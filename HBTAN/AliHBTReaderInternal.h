#ifndef ALIHBTREADERINTERNAL_H
#define ALIHBTREADERINTERNAL_H
//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//                                                                       //
// class AliHBTReaderInternal                                            //
//                                                                       //
// Multi file reader for Internal Data Format                            //
//                                                                       //
// This reader reads data created by itself                              //
//   (method AliHBTReaderInternal::Write)                                //
//                                                                       //
// Piotr.Skowronski@cern.ch                                              //
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliHBTReader.h"
#include <TString.h>

class TFile;
class TTree;
class TBranch;
class TClonesArray;

class AliHBTReaderInternal: public AliHBTReader
{
  public:
    AliHBTReaderInternal();
    AliHBTReaderInternal(const char *filename);
    AliHBTReaderInternal(TObjArray* dirs, const char *filename);
    AliHBTReaderInternal(const AliHBTReaderInternal& in);
    
    virtual ~AliHBTReaderInternal();
    
    AliHBTReaderInternal& operator=(const AliHBTReaderInternal& in);
    
    void          Rewind();
    
    Bool_t        ReadsTracks() const {return kTRUE;}
    Bool_t        ReadsParticles() const {return kTRUE;}
    
    static Int_t Write(AliHBTReader* reader,const char* outfile = "data.root", Bool_t multcheck = kFALSE);//reads tracks from runs and writes them to file
    
  protected:
    
    TString       fFileName;//name of the file with tracks
    TBranch*      fPartBranch;//!branch with particles
    TBranch*      fTrackBranch;//!branch with tracks
    TTree*        fTree;//!tree
    TFile*        fFile;//!file
    TClonesArray* fPartBuffer;//!buffer array that tree is read to
    TClonesArray* fTrackBuffer;//!
        
    Int_t         ReadNext();//reads next event
    Int_t         OpenNextFile();//opens file to be read for given event 
    static Bool_t FindIndex(TClonesArray* arr,Int_t idx);
    
    ClassDef(AliHBTReaderInternal,2)
};
#endif 
