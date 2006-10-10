#include "AliAODRun.h"
//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliAODRun                                         //
//                                                   //
// Class storing and managing events                 //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://aliweb.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

#include <TObjArray.h>

ClassImp(AliAODRun)
/**************************************************************************/ 

AliAODRun::AliAODRun():
  fEvents(new TObjArray())
{ 
 //contructor
  if(!fEvents) Fatal("AliAODRun::AliAODRun","Can not allocate memory");
  fEvents->SetOwner(); //array is an owner: when is deleted or cleared it deletes objects that it contains
}
/**************************************************************************/

AliAODRun::~AliAODRun()
{
  //destructor
  delete fEvents;//delete array with events
}
/**************************************************************************/

void AliAODRun::Reset()
 { 
   fEvents->Clear();//clear an array with events. 
                    //All events are deleted because array is an owner (set in constructor)
 }
/**************************************************************************/

void AliAODRun::AddParticle(Int_t event, AliVAODParticle* part)
{
 //Adds particle to event
 //if there is no event of this number, crate it and add to the collection
 if(!GetEvent(event))  fEvents->AddAtAndExpand(new AliAOD, event);
 
 GetEvent(event)->AddParticle(part);
}
/**************************************************************************/

void AliAODRun::AddParticle(Int_t event, TParticle* part, Int_t idx)
{
 //if there is no event of this number, crate it and add to the collection
 if(!GetEvent(event))  fEvents->AddAtAndExpand(new AliAOD, event);
 GetEvent(event)->AddParticle(part,idx);
}
/**************************************************************************/

void AliAODRun::AddParticle(Int_t event, Int_t pdg, Int_t idx,
                            Double_t px, Double_t py, Double_t pz, Double_t etot,
                            Double_t vx, Double_t vy, Double_t vz, Double_t time)
{
 //if there is no event of this number, crate it and add to the collection
 if(!GetEvent(event))  fEvents->AddAtAndExpand(new AliAOD, event);
 GetEvent(event)->AddParticle(pdg,idx,px,py,pz,etot,vx,vy,vz,time);
}
/**************************************************************************/ 

void AliAODRun::SetEvent(Int_t number, AliAOD* event)
{
  //adds an event to the run
  if (event == 0x0)
   {
     delete fEvents->RemoveAt(number);
     return;
   }
  AliAOD* ev = GetEvent(number);
  if (ev == event) return;
  
  delete fEvents->RemoveAt(number);
  fEvents->AddAtAndExpand(event, number);
  
}
/**************************************************************************/ 

