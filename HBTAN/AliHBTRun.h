#ifndef ALIHBTRUN_H
#define ALIHBTRUN_H
//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTRun                                         //
//                                                   //
// Class storing and managing set events             //
// designed for fast acces                           //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////


#include "AliHBTEvent.h"
#include <TObjArray.h>

class AliHBTParticle;

class AliHBTEvent;

class AliHBTRun: public TObject
 {
  public:
    AliHBTRun();
    virtual ~AliHBTRun();

    void            AddParticle(Int_t event, AliHBTParticle* part); //inerface to AliHBTEvent::AddParticle(AliHBTParticle*) 
    void            AddParticle(Int_t event, TParticle* part, Int_t idx);//inerface to AliHBTEvent::AddParticle(TParticle*) 
    
    //inerface to AliHBTEvent::AddParticle(Int_t.Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t)
    void            AddParticle(Int_t event, Int_t pdg, Int_t idx,
                                Double_t px, Double_t py, Double_t pz, Double_t etot,
                                Double_t vx, Double_t vy, Double_t vz, Double_t time); 
    
    void            SetEvent(Int_t number, AliHBTEvent* event);
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
  //check if array is enough big - protect from error massage from array "Out Of Bounds"
  if (event>=fEvents->GetSize()) return 0x0;//WARNING, that line relies 
                                            //on index of first object in TObjArray is 0
		    //== LowerBound = 0
  return (AliHBTEvent*)fEvents->At(event);
}
/**************************************************************************/
inline
AliHBTParticle* AliHBTRun::GetParticle(Int_t event, Int_t n) 
{
 //returns nth particle from event number event
  AliHBTEvent* e = GetEvent(event);
  return (e)?e->GetParticle(n):0x0;
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
  AliHBTEvent* e = GetEvent(event);
  return (e)?e->GetNumberOfParticles():0x0;
}
 
#endif
