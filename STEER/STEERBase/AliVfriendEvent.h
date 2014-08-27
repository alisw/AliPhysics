#ifndef ALIVFRIENDEVENT_H
#define ALIVFRIENDEVENT_H

#include "Rtypes.h"
class AliVfriendTrack;

//_____________________________________________________________________________
class AliVfriendEvent {
public:
  AliVfriendEvent() {}
  virtual ~AliVfriendEvent() {}

  virtual Int_t GetNclustersTPC(UInt_t /*sector*/) const = 0;
  virtual Int_t GetNclustersTPCused(UInt_t /*sector*/) const = 0;

  //used in calibration
  virtual Bool_t TestSkipBit() const = 0;
  virtual Int_t GetNumberOfTracks() const = 0;
  virtual const AliVfriendTrack *GetTrack(Int_t /*i*/) const = 0;

private: 
  AliVfriendEvent(const AliVfriendEvent &);
  AliVfriendEvent& operator=(const AliVfriendEvent& esd);  

//  ClassDef(AliVfriendEvent,0);
};

#endif


