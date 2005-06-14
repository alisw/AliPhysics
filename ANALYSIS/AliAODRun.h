#ifndef ALIAODRUN_H
#define ALIAODRUN_H
//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliAODRun                                         //
//                                                   //
// Class storing and managing set events             //
// designed for fast acces                           //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://aliweb.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////


#include "AliAOD.h"
#include <TObjArray.h>

class AliVAODParticle;
class TParticle;

class AliAODRun: public TObject
 {
  public:
    AliAODRun();
    virtual ~AliAODRun();

    void            AddParticle(Int_t event, AliVAODParticle* part); //inerface to AliAOD::AddParticle(AliVAODParticle*) 
    void            AddParticle(Int_t event, TParticle* part, Int_t idx);//inerface to AliAOD::AddParticle(TParticle*) 
    
    //inerface to AliAOD::AddParticle(Int_t.Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t)
    void            AddParticle(Int_t event, Int_t pdg, Int_t idx,
                                Double_t px, Double_t py, Double_t pz, Double_t etot,
                                Double_t vx, Double_t vy, Double_t vz, Double_t time); 
    
    void            SetEvent(Int_t number, AliAOD* event);
    AliVAODParticle* GetParticle(Int_t event, Int_t n); //returns nth particle from event
    AliAOD*    GetEvent(Int_t event) const; //returns AliAOD number "event"
    
    Int_t           GetNumberOfEvents() const; //returns number of events
    Int_t         	GetNumberOfParticlesInEvent(Int_t event) const; //returns number of particles in event number "event"
    void            Reset();//clears all events in the array (deletes)
  protected:
    TObjArray* fEvents;//!Array containig AliAODs
  private:
    
  public:
    ClassDef(AliAODRun,1)
 };
 
 
/**************************************************************************/

inline
AliAOD* AliAODRun::GetEvent(Int_t event) const
{
//returns pointer to AliAOD number "event"
  //check if array is enough big - protect from error massage from array "Out Of Bounds"
  if (event>=fEvents->GetSize()) return 0x0;//WARNING, that line relies 
                                            //on index of first object in TObjArray is 0
		    //== LowerBound = 0
  return (AliAOD*)fEvents->At(event);
}
/**************************************************************************/
inline
AliVAODParticle* AliAODRun::GetParticle(Int_t event, Int_t n) 
{
 //returns nth particle from event number event
  AliAOD* e = GetEvent(event);
  return (e)?e->GetParticle(n):0x0;
}

/**************************************************************************/

inline
Int_t AliAODRun::GetNumberOfEvents() const
 {
//returns number of events in collection

   return fEvents->GetEntriesFast(); //there may be empty slots but we do not care
                                     //Analysis checks it if return is not NULL
 }
/**************************************************************************/

inline
Int_t AliAODRun::GetNumberOfParticlesInEvent(Int_t event) const
{
//returns number of Particles in event 
  AliAOD* e = GetEvent(event);
  return (e)?e->GetNumberOfParticles():0x0;
}
 
#endif
