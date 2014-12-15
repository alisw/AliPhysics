#ifndef ALI_ITS_ONLINESPDHITEVENT_H
#define ALI_ITS_ONLINESPDHITEVENT_H

/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// One object for each half stave and step in a scan. It keeps //
// the nr of events with at least one pixel hit in each chip.  //
// It also keeps the value for the all the 10 chips together.  //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscan class.                                  //
/////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliITSOnlineSPDHitEvent : public TObject {

 public:
  AliITSOnlineSPDHitEvent();
  virtual ~AliITSOnlineSPDHitEvent(){}
  void   IncrementHitEvent(UInt_t chip);
  void   SetHitEvent(UInt_t chip, UInt_t events);
  UInt_t GetHitEvent(UInt_t chip) const;
  AliITSOnlineSPDHitEvent* CloneThis() const;

 private:
  UInt_t fHitEvent[11];  // nr of events with at least one hit in a chip
                         // index 10 is for all 10 chips together

  ClassDef(AliITSOnlineSPDHitEvent,1)
   };

#endif
