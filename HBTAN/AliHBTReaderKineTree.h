#ifndef ALIHBTREADERKINETREE_H
#define ALIHBTREADERKINETREE_H

#include "AliHBTReader.h"
#include <TString.h>

class TFile;
class AliStack;
class AliRunLoader;

class AliHBTReaderKineTree: public AliHBTReader
 {
   public:
    AliHBTReaderKineTree();
    
    AliHBTReaderKineTree(TString&);
    AliHBTReaderKineTree(TObjArray*,const Char_t *filename="galice.root");

    virtual ~AliHBTReaderKineTree();
    
    void          Rewind();
    
    Bool_t        ReadsTracks() const {return kFALSE;}
    Bool_t        ReadsParticles() const {return kTRUE;}
    
   protected:
    Int_t         ReadNext();//reads tracks and particles and puts them in runs
    Int_t         OpenNextFile();
   
    TString       fFileName;//file name 
    AliRunLoader* fRunLoader;
    
    static const TString fgkEventFolderName;
    
   private:
   
   public:
     ClassDef(AliHBTReaderKineTree,2)
 };

#endif
