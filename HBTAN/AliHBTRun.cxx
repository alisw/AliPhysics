#include "AliHBTRun.h"

#include <TObjArray.h>

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTRun                                         //
//                                                   //
// Class storing and managing events                 //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

ClassImp(AliHBTRun)
/**************************************************************************/ 

AliHBTRun::AliHBTRun()
{ 
 //contructor
  fEvents = new TObjArray();//create array for AliHBTEvents
  if(!fEvents) Fatal("AliHBTRun::AliHBTRun","Can not allocate memory");
  fEvents->SetOwner(); //array is an owner: when is deleted or cleared it deletes objects that it contains
}
/**************************************************************************/

AliHBTRun::~AliHBTRun()
{
  //destructor
  delete fEvents;//delete array with events
}
/**************************************************************************/

void AliHBTRun::Reset()
 { 
   fEvents->Clear();//clear an array with events. 
                    //All events are deleted because array is an owner (set in constructor)
 }
/**************************************************************************/

void AliHBTRun::AddParticle(Int_t event, AliHBTParticle* part)
{
 //Adds particle to event
 //if there is no event of this number, crate it and add to the collection
 if(!GetEvent(event))  fEvents->AddAtAndExpand(new AliHBTEvent, event);
 
 GetEvent(event)->AddParticle(part);
}
/**************************************************************************/

void AliHBTRun::AddParticle(Int_t event, TParticle* part, Int_t idx)
{
 //if there is no event of this number, crate it and add to the collection
 if(!GetEvent(event))  fEvents->AddAtAndExpand(new AliHBTEvent, event);
 GetEvent(event)->AddParticle(part,idx);
}
/**************************************************************************/

void AliHBTRun::AddParticle(Int_t event, Int_t pdg, Int_t idx,
                            Double_t px, Double_t py, Double_t pz, Double_t etot,
                            Double_t vx, Double_t vy, Double_t vz, Double_t time)
{
 //if there is no event of this number, crate it and add to the collection
 if(!GetEvent(event))  fEvents->AddAtAndExpand(new AliHBTEvent, event);
 GetEvent(event)->AddParticle(pdg,idx,px,py,pz,etot,vx,vy,vz,time);
}
/**************************************************************************/ 

void AliHBTRun::SetEvent(Int_t number, AliHBTEvent* event)
{
  //adds an event to the run
  // makes an own copy of the event!
  if (event == 0x0)
   {
     delete fEvents->RemoveAt(number);
     return;
   }
  AliHBTEvent* ev = GetEvent(number);
  if (ev) *ev = *event;
  else fEvents->AddAtAndExpand(new AliHBTEvent(*event), number);
}
/**************************************************************************/ 

