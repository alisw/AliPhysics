#ifndef ALIHBTEVENT_H
#define ALIHBTEVENT_H
//__________________________________________________________
///////////////////////////////////////////////////////////////////
//
// class AliHBTEvent
//
// This class is container for paticles coming from one event
// -----------------------------------------------------------------
// -----------------------------------------------------------------
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
//
// Piotr.Skowronski@cern.ch 
//
///////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliHBTParticle;
class TParticle;
class AliHBTEvent: public TObject
 {
  public:
    AliHBTEvent();
    AliHBTEvent(const AliHBTEvent& source);
    virtual ~AliHBTEvent();
    
    AliHBTEvent& operator=(const AliHBTEvent& source);
        
    static const UInt_t fgkInitEventSize; //initial number of the array
                                          //if expanded, this size is used also
    AliHBTParticle* GetParticle(Int_t n);  //gets particle 
    AliHBTParticle* GetParticleSafely(Int_t n); //gets particle with index check
    
    void    AddParticle(AliHBTParticle* hbtpart); //adds particle to the event
    void    AddParticle(TParticle* part, Int_t idx); //adds particle to the event
    void    AddParticle(Int_t pdg, Int_t idx, Double_t px, Double_t py, Double_t pz, Double_t etot,
                        Double_t vx, Double_t vy, Double_t vz, Double_t time);
    
    Int_t   GetNumberOfParticles() const;
    void    Reset(); //deletes all entries
    void    SetOwner(Bool_t owns = kTRUE){ fOwner = owns; }
    Bool_t  IsOwner() const {return fOwner;}
    void    SetRandomized(Bool_t rd = kTRUE){fRandomized = rd;}
    Bool_t  IsRandomized()const {return fRandomized;}
    void    SwapParticles(Int_t i, Int_t j);//swaps particles positions; used by AliHBTEvent::Blend
    
  protected:
    Int_t             fSize;       //!current size of the array
    AliHBTParticle ** fParticles; //!array of pointers to the particles
    Int_t             fNParticles; //!number of particles in Event
    Bool_t            fOwner;      //flag if that event owns the 
    Bool_t            fRandomized; //!flag indicating if particles positions has been already randomizd
    void              Expand();    //expands the array if necessary

  private:
    ClassDef(AliHBTEvent,1)
 };
/**************************************************************************/ 
 
inline 
AliHBTParticle* AliHBTEvent::GetParticle(Int_t n)
 {
 //Gets particle without boundary check
   return fParticles[n];
 }

/**************************************************************************/ 
 
inline 
Int_t  AliHBTEvent::GetNumberOfParticles() const
 {
 //reurns number of particles in this event
   return fNParticles;
 }
 
#endif
