#ifndef ALIREADERAOD_H
#define ALIREADERAOD_H

#include "AliReader.h"

class TTree;
class TFile;

class AliReaderAOD: public AliReader
{
  public:
    AliReaderAOD(const Char_t* aodfilename = "AOD.root");
    virtual ~AliReaderAOD();

    void          ReadSimulatedData(Bool_t flag){fReadSim = flag;}//switches reading MC data
    Bool_t        ReadsRec() const {return kTRUE;}
    Bool_t        ReadsSim() const {return fReadSim;}

    void          Rewind();


    static Int_t WriteAOD(AliReader* reader, const char* outfilename = "AliAOD.root", //reads tracks from runs and writes them to file
                          const char* pclassname = "AliAODParticle", Bool_t multcheck = kFALSE);
    
  protected:
    virtual Int_t         ReadNext();
    virtual Int_t        OpenFile(Int_t evno);//opens files to be read for given event

    static const TString  fgkTreeName;//name of branch holding simulated data
    static const TString  fgkRecosntructedDataBranchName;//name of branch holding reconstructed data
    static const TString  fgkSimulatedDataBranchName;//name of branch holding simulated data
    
  
  private:
    TString fFileName;//File name
    
    Bool_t  fReadSim;//indicates if to read simulated data

    TTree*        fTree;//!tree
    TFile*        fFile;//!file
    AliAOD*       fSimBuffer;//!buffer array that tree is read to
    AliAOD*       fRecBuffer;//!
    
    
    ClassDef(AliReaderAOD,1)
};

#endif
