#ifndef AliReaderESDTree_H
#define AliReaderESDTree_H
//_______________________________________________________________________
/////////////////////////////////////////////////////////////////////////
//
// class AliReaderESDTree
//
// Reader for ESD Tree 
//
// Ch. Finck
//
/////////////////////////////////////////////////////////////////////////
#include "AliReaderESD.h"
#include <TString.h>

class TFile;
class TTree;

class AliReaderESDTree: public AliReaderESD
 {
   public:

    AliReaderESDTree(const Char_t* esdfilename = "AliESDs.root", 
                     const Char_t* galfilename = "galice.root");

    virtual ~AliReaderESDTree();


   protected:
    Int_t         ReadNext();//reads tracks and particles and puts them in runs
    TFile*        OpenFile(Int_t evno);//opens files to be read for given event
   
    TTree*        fTree;// tree pointer
    
   private:
    ClassDef(AliReaderESDTree,1)
 };

#endif
