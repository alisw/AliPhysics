#ifndef ALIHBTREADERKINETREE_H
#define ALIHBTREADERKINETREE_H
//_______________________________________________________________________
/////////////////////////////////////////////////////////////////////////
//
// class AliHBTReaderKineTree
//
// Reader for Kinematics
//
// Piotr.Skowronski@cern.ch
//
/////////////////////////////////////////////////////////////////////////
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
    AliHBTReaderKineTree(const AliHBTReaderKineTree& in);

    virtual ~AliHBTReaderKineTree();

    AliHBTReaderKineTree& operator=(const AliHBTReaderKineTree& in);
    
    void          Rewind();
    
    Bool_t        ReadsTracks() const {return kFALSE;}
    Bool_t        ReadsParticles() const {return kTRUE;}
    
   protected:
    Int_t         ReadNext();//reads tracks and particles and puts them in runs
    Int_t         OpenNextFile();
   
    TString       fFileName;//file name 
    AliRunLoader* fRunLoader;//!Pointer to loader
    
    static const TString fgkEventFolderName; //Event folder name that session are mounter
    
   private:
     ClassDef(AliHBTReaderKineTree,2)
 };

#endif
