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

    virtual ~AliHBTReaderKineTree(){}
    
    Int_t        Read(AliHBTRun* particles, AliHBTRun* /*tracks*/);//reads tracks and particles and puts them in runs

    AliHBTEvent* GetParticleEvent(Int_t);//returns pointer to event with particles
    AliHBTEvent* GetTrackEvent(Int_t){return 0x0;}//returns pointer to event with particles
    Int_t        GetNumberOfPartEvents();//returns number of particle events
    Int_t        GetNumberOfTrackEvents(){return 0;}//returns number of track events

    
   protected:
    TString       fFileName;
    AliHBTRun*    fParticles; //!simulated particles

    AliRunLoader* OpenFile(Int_t);

    Bool_t fIsRead;//!flag indicating if the data are already read
    
    static const TString fgkEventFolderName;
    
   private:
   
   public:
     ClassDef(AliHBTReaderKineTree,1)
 };

#endif
