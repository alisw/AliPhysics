#ifndef ALIHBTREADER_H
#define ALIHBTREADER_H

#include <TObject.h>

//Reader Base class (reads particles and tracks and
//puts it to the AliHBTRun objects
//Piotr.Skowronski@cern.ch

class AliHBTRun;
class AliHBTEvent;
class AliHBTParticleCut;
class TObjArray;
class AliHBTParticle;
class TString;

class AliHBTReader: public TObject

{
  public:
    AliHBTReader();
    AliHBTReader(TObjArray*);
    virtual ~AliHBTReader();
    //in the future this class is will read global tracking
    virtual Int_t Read(AliHBTRun* particles, AliHBTRun *tracks) = 0;
    
    virtual AliHBTEvent* GetParticleEvent(Int_t) = 0;
    virtual AliHBTEvent* GetTrackEvent(Int_t) = 0;
    virtual Int_t GetNumberOfPartEvents() = 0;
    virtual Int_t GetNumberOfTrackEvents() = 0;
    
    void AddParticleCut(AliHBTParticleCut* cut);
    
    void SetDirs(TObjArray* dirs){fDirs = dirs;} //sets array directories names

  protected:
    
    TObjArray *fCuts;//array with particle cuts
    TObjArray *fDirs;//arry with directories to read data from

    Bool_t Pass(AliHBTParticle*);
    Bool_t Pass(Int_t pid);

    TString& GetDirName(Int_t);
    
  private:
  
  public:
    ClassDef(AliHBTReader,1)
};

#endif
