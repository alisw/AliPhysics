#ifndef ALIHBTRUN_H
#define ALIHBTRUN_H

#include "AliHBTEvent.h"
#include <TObjArray.h>

//class describing set of events (the run)
//designed for fast acces
//Piotr.Skowronski@cern.ch

class AliHBTParticle;

class AliHBTEvent;

class AliHBTRun: public TObject
 {
  public:
    AliHBTRun();
    virtual ~AliHBTRun();

    void            AddParticle(Int_t event, AliHBTParticle* part); //inerface to AliHBTEvent::AddParticle(AliHBTParticle*) 
    void            AddParticle(Int_t event, TParticle* part);//inerface to AliHBTEvent::AddParticle(TParticle*) 
    
    //inerface to AliHBTEvent::AddParticle(Int_t.Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t)
    void            AddParticle(Int_t event, Int_t pdg, 
                                Double_t px, Double_t py, Double_t pz, Double_t etot,
                                Double_t vx, Double_t vy, Double_t vz, Double_t time); 
    
    AliHBTParticle* GetParticle(Int_t event, Int_t n); //returns nth particle from event
    AliHBTEvent*    GetEvent(Int_t event) const; //returns AliHBTEvent number "event"
    
    Int_t           GetNumberOfEvents() const; //returns number of events
    Int_t         	GetNumberOfParticlesInEvent(Int_t event) const; //returns number of particles in event number "event"
    void            Reset();//clears all events in the array (deletes)
  protected:
    TObjArray* fEvents;//!Array containig AliHBTEvents
  private:
    
  public:
    ClassDef(AliHBTRun,1)
 };
 
 
/**************************************************************************/

inline
AliHBTEvent* AliHBTRun::GetEvent(Int_t event) const
{
//returns pointer to AliHBTEvent number "event"
  
  return (AliHBTEvent*)fEvents->At(event);
}
/**************************************************************************/
inline
AliHBTParticle* AliHBTRun::GetParticle(Int_t event, Int_t n) 
{
 //returns nth particle from event number event
  return GetEvent(event)->GetParticle(n);
}

/**************************************************************************/

inline
Int_t AliHBTRun::GetNumberOfEvents() const
 {
//returns number of events in collection

   return fEvents->GetEntriesFast(); //there may be empty slots but we do not care
                                     //Analysis checks it if return is not NULL
 }
/**************************************************************************/

inline
Int_t AliHBTRun::GetNumberOfParticlesInEvent(Int_t event) const
{
//returns number of Particles in event 
  return GetEvent(event)->GetNumberOfParticles();
}
 
#endif
