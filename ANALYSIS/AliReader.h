#ifndef ALIREADER_H
#define ALIREADER_H
//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//
// class AliReader
//
// Reader Base class 
// Reads particles and tracks and
// puts it to the AliAOD objects and eventuall buffers in AliAODRuns
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TObjArray.h>

class AliAODRun;
class AliAOD;
class AliAODParticleCut;
class AliAODParticle;
class TString;
class TH1I;
 
class AliReader: public TNamed
{
  public:
    AliReader();
    AliReader(TObjArray*);
    AliReader(const AliReader& in);
    virtual ~AliReader();
    
    AliReader& operator=(const AliReader& in);
        
    virtual Int_t        Next();
    virtual void         Rewind() = 0; //
    
    virtual Bool_t       ReadsSim() const = 0; //specifies if reader is able to read simulated particles
    virtual Bool_t       ReadsRec() const = 0;//specifies if reader is able to read reconstructed tracks
    
    void                 AddParticleCut(AliAODParticleCut* cut);
    
    virtual Int_t        Read(AliAODRun* particles, AliAODRun *tracks);
    
    virtual AliAOD*      GetEventRec() const {return fEventRec;}//
    virtual AliAOD*      GetEventSim() const {return fEventSim;}//can not be const because position randomizer overloads it
    
    virtual AliAOD*      GetEventRec(Int_t n);
    virtual AliAOD*      GetEventSim(Int_t n);
    
    virtual Int_t        GetNumberOfRecEvents();
    virtual Int_t        GetNumberOfSimEvents();
    
    void                 SetDirs(TObjArray* dirs){fDirs = dirs;} //sets array directories names
    void                 SetEventBuffering(Bool_t flag){fBufferEvents = flag;}//switches on/off buffering - read data are kept in local buffer
    void                 SetBlend(Bool_t flag = kTRUE){fBlend=flag;} //set blending - randomizing particle order
    virtual Int_t        GetNumberOfDirs() const {return (fDirs)?fDirs->GetEntries():0;}
    void                 ReadEventsFromTo(Int_t first,Int_t last){fFirst = first; fLast = last;}
    virtual TH1I*        GetTrackCounter() const {return fTrackCounter;}
    virtual void         WriteTrackCounter() const;
    
  protected:
    
    TObjArray*           fCuts;//array with particle cuts
    TObjArray*           fDirs;//arry with directories to read data from
    
    Int_t                fCurrentEvent;//!  number of current event in current directory
    Int_t                fCurrentDir;//! number of current directory
    
    Int_t                fNEventsRead;//!total 
        
    AliAOD*              fEventRec;    //! tracks read from current event
    AliAOD*              fEventSim;    //! particles read from current event
    
    AliAODRun*           fRunSim; //!simulated particles
    AliAODRun*           fRunRec; //!reconstructed tracks
    
    Bool_t               fIsRead;//!flag indicating if the data are already read
    Bool_t               fBufferEvents;//flag indicating if the data should be bufferred
    
    Bool_t               fBlend;// flag indicating if randomly change positions of the particles after reading
    
    Int_t                fFirst;//first event to return (all are before are skipped)
    Int_t                fLast;//last

    TH1I*                fTrackCounter; //histogram with number of tracks read
    
    virtual Int_t        ReadNext() = 0; //this methods reads next event and put result in fTracksEvent and/or fParticlesEvent
    Bool_t               Pass(AliAODParticle* p);
    Bool_t               Pass(Int_t pid);
    void                 Blend();
    
    TString&             GetDirName(Int_t entry);
    
  private:
  
    ClassDef(AliReader,1)//
};

#endif
