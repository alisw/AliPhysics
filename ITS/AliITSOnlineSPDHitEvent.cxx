/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// One object for each half stave and step in a scan. It keeps //
// the nr of events with at least one pixel hit in each chip.  //
// It also keeps the value for the all the 10 chips together.  //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscan class.                                  //
/////////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDHitEvent.h"

ClassImp(AliITSOnlineSPDHitEvent)

AliITSOnlineSPDHitEvent::AliITSOnlineSPDHitEvent(){
  // constructor, sets the nr of events hit for each chip to 0
  for (Int_t i=0; i<11; i++) {
    fHitEvent[i]=0;
  }
}

AliITSOnlineSPDHitEvent* AliITSOnlineSPDHitEvent::CloneThis() const {
  // makes a copy of this object and returns it
  AliITSOnlineSPDHitEvent* returnpointer = new AliITSOnlineSPDHitEvent();
  for (Int_t i=0; i<11; i++) {
    returnpointer->SetHitEvent(i,fHitEvent[i]);
  }
  return returnpointer;
}

void   AliITSOnlineSPDHitEvent::IncrementHitEvent(UInt_t chip) {
  // increment the nr of hit events for chip 'chip'
  if (chip<=10) {
    fHitEvent[chip]++;
  }
}
void   AliITSOnlineSPDHitEvent::SetHitEvent(UInt_t chip, UInt_t events) {
  // set the nr of hit events for chip 'chip'
  if (chip<=10) {
    fHitEvent[chip] = events;
  }
}
UInt_t AliITSOnlineSPDHitEvent::GetHitEvent(UInt_t chip) const {
  // get the nr of hit events for chip 'chip'
  if (chip<=10) {
    return fHitEvent[chip];
  }
  else {
    return 0;
  }
}
