#ifndef ALIHBTEvent_H
#define ALIHBTEvent_H
//This class sters HBT perticles for one event
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html

#include <TObject.h>

class AliHBTParticle;
class TParticle;
class AliHBTEvent: public TObject
 {
  public:
    AliHBTEvent();
    virtual ~AliHBTEvent();
    const static UInt_t fgkInitEventSize; //initial number of the array
                                         //if expanded, this size is used also
    AliHBTParticle* GetParticle(Int_t n);  //gets particle 
    AliHBTParticle* GetParticleSafely(Int_t n); //gets particle with index check
    
    void    AddParticle(AliHBTParticle*); //adds particle to the event
    void    AddParticle(TParticle*); //adds particle to the event
    void    AddParticle(Int_t pdg, Double_t px, Double_t py, Double_t pz, Double_t etot,
                        Double_t vx, Double_t vy, Double_t vz, Double_t time);
    
    Int_t   GetNumberOfParticles() const;
    void    Reset(); //deletes all entries
    void    SetOwner(Bool_t owns = kTRUE){ fOwner = owns; }
    Bool_t  IsOwner() {return fOwner;}
  protected:
    AliHBTParticle ** fParticles; //!array of pointers to the particles
    Int_t  fNParticles; //!number of particles in Event
    Int_t  fSize;       //!current size of the array
    Bool_t fOwner;      //flag if that event owns the 
    void   Expand();    //expands the array if necessary
  private:
    
  public:
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
