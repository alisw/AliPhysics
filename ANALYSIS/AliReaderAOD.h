#ifndef ALIREADERAOD_H
#define ALIREADERAOD_H

#include "AliReader.h"

class AliReaderAOD: public AliReader
{
  public:
    AliReaderAOD(const Char_t* aodfilename = "AliAOD.root"){}
    virtual ~AliReaderAOD(){}

    void          ReadSimulatedData(Bool_t flag){fReadSim = flag;}//switches reading MC data
    Bool_t        ReadsRec() const {return kTRUE;}
    Bool_t        ReadsSim() const {return fReadSim;}


    static Int_t WriteAOD(AliReader* reader, const char* outfilename = "AliAOD.root", Bool_t multcheck = kFALSE);//reads tracks from runs and writes them to file
    
  protected:
  private:
    TString fFileName;//File name
    
    Bool_t  fReadSim;//indicates if to read simulated data
    ClassDef(AliReaderAOD,1)
};

#endif
