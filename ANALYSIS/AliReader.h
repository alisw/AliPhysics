#ifndef ALIREADER_H
#define ALIREADER_H
//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//
// class AliReader
//
// Reader Base class 
// Reads particles and tracks and
// puts them to the AliAOD objects and eventually, if needed, buffers AliAODs in AliAODRun(s)
//
// User loops over events calling method Next. In case of success this method returns 0.
// In case of error or if there is no more events to read, non-0 value is returned
//
// Reading can be rewound to the beginning using method Rewind.
//
// Tracks are read to the fEventRec (contains reconstructed tracks) 
// and fEventSim (corresponding MC simulated data) data members,
// that are of the type AliAOD. 
//
// If a given reader has ability of reading both, reconstructed and simulated data, 
// these are structured in AODs so a "n'th" simulated particle 
// (the one stored in the fEventSim at slot n) 
// corresponds to the n'th reconstructed track (the one stored in the fEventRec at slot n).
//
// The same reconstructed track can be present more than ones in the AOD,
// but with a different PID. In this case
// pointer to the corresponding MC simulated particles is also present more than ones.
// This situation happens if you want to read all particles 
// with PID probability of being , e.g.,  pion higher than 60%
// and being kaon higher than 40%. Than, if a given track has probability Ppid(pi)=52% and Ppid(K)=48% 
// than it is read twise.
//
// Provides functionality for both buffering and non-buffering reading
// This can be switched on/off via method SetEventBuffering(bool)
// The main method that inheriting classes need to implement is ReadNext()
// that read next event in queue.
//
// The others are:
// Bool_t  ReadsRec() const; specifies if reader is able to read simulated particles
// Bool_t  ReadsSim() const; specifies if reader is able to read reconstructed tracks
// void    Rewind(); rewind reading to the beginning
//
// This class provides full functionality for reading from many sources
// User can provide TObjArray of TObjStrings (SetDirs method or via parameter 
// in the constructor) which desribes paths of directories to search data in.
// If none specified current directory is searched.
// 
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TObjArray.h>

class TGliteXmlEventlist;
    
class AliAODRun;
class AliAOD;
class AliAODParticleCut;
class AliVAODParticle;
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

    void                 AddParticleCut(AliAODParticleCut* cut);//adds a particle cut to the list of cuts
    
    virtual AliAOD*      GetEventRec() const {return fEventRec;}//returns current event with reconstructed tracks
    virtual AliAOD*      GetEventSim() const {return fEventSim;}//returns current event with simulated particles
    
    virtual AliAOD*      GetEventRec(Int_t n);//returns event number n
    virtual AliAOD*      GetEventSim(Int_t n);
    
    virtual Int_t        Read(const char * name) {return TObject::Read(name);}
    virtual Int_t        Read(AliAODRun* particles, AliAODRun *tracks);//Reads all available evenets and stores them in 'particles' and 'tracks'

    virtual Int_t        GetNumberOfRecEvents();//Returns number of available events -> usually conncected with reading all events
                                                //may be time consuming
    virtual Int_t        GetNumberOfSimEvents();// 
     
    void                 SetEventList(TGliteXmlEventlist* evl){fEventList = evl;}
         
    void                 SetDirs(TObjArray* dirs){fDirs = dirs;} //sets array directories names
    void                 SetEventBuffering(Bool_t flag){fBufferEvents = flag;}//switches on/off buffering - read data are kept in local buffer
    void                 SetBlend(Bool_t flag = kTRUE){fBlend=flag;} //set blending - randomizing particle order
    virtual Int_t        GetNumberOfDirs() const {return (fDirs)?fDirs->GetEntries():0;}
    void                 ReadEventsFromTo(Int_t first,Int_t last){fFirst = first; fLast = last;}
    virtual TH1I*        GetTrackCounter() const {return fTrackCounter;}
    virtual void         WriteTrackCounter() const;//Writes the track counting histigram 
    
  protected:
    
    TGliteXmlEventlist*  fEventList;//Event list delivered by GLite/AliEn
    
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
    
    Int_t                fFirst;//first event to return (all before are skipped)
    Int_t                fLast;//the last one

    TH1I*                fTrackCounter; //histogram with number of tracks read
    
    virtual Int_t        ReadNext() = 0; //this methods reads next event and put result in fTracksEvent and/or fParticlesEvent
    Bool_t               Rejected(AliVAODParticle* p);//Checks if a given particle agains cuts
    Bool_t               Rejected(Int_t pid);//Checks if a given pid passes cuts
    void                 Blend();//Mixes current events in a symmetric way so after mixing thy are consistent
    
    TString              GetDirName(Int_t entry);
    
  private:
  
    ClassDef(AliReader,1)//
};

#endif
