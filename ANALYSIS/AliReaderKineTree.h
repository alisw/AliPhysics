#ifndef AliReaderKineTree_H
#define AliReaderKineTree_H
//_______________________________________________________________________
/////////////////////////////////////////////////////////////////////////
//
// class AliReaderKineTree
//
// Reader for Kinematics
//
// Piotr.Skowronski@cern.ch
//
/////////////////////////////////////////////////////////////////////////
#include "AliReader.h"
#include <TString.h>

class TFile;
class AliStack;
class AliRunLoader;

class AliReaderKineTree: public AliReader
 {
   public:
    AliReaderKineTree();
    
    AliReaderKineTree(TString&);
    AliReaderKineTree(TObjArray*,const Char_t *filename="galice.root");
    AliReaderKineTree(const AliReaderKineTree& in);

    virtual ~AliReaderKineTree();

    AliReaderKineTree& operator=(const AliReaderKineTree& in);
    
    void          Rewind();
    
    Bool_t        ReadsRec() const {return kFALSE;}
    Bool_t        ReadsSim() const {return kTRUE;}
    
   protected:
    Int_t         ReadNext();//reads tracks and particles and puts them in runs
    Int_t         OpenNextFile();
   
    TString       fFileName;//file name 
    AliRunLoader* fRunLoader;//!Pointer to loader
    
    static const TString fgkEventFolderName; //Event folder name that session are mounter
    
   private:
     ClassDef(AliReaderKineTree,2)
 };

#endif
