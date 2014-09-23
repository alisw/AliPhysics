#ifndef ALIVFRIENDEVENT_H
#define ALIVFRIENDEVENT_H

#include "Rtypes.h"
#include "AliVMisc.h"
class AliVfriendTrack;

//_____________________________________________________________________________
class AliVfriendEvent {
public:
  AliVfriendEvent() {}
  virtual ~AliVfriendEvent() {}

  // constructor and method for reinitialisation of virtual table
  AliVfriendEvent( AliVConstructorReinitialisationFlag );
  void Reinitialize(){} // do nothing

  virtual Int_t GetNumberOfTracks() const = 0;
  virtual const AliVfriendTrack *GetTrack(Int_t /*i*/) const = 0;
  virtual Int_t GetEntriesInTracks() const = 0;

  // AliESDVZEROfriend *GetVZEROfriend();
  // AliESDTZEROfriend *GetTZEROfriend();

  virtual void Ls() const = 0;
  virtual void Reset() = 0;

  // bit manipulation for filtering
  virtual void SetSkipBit(Bool_t skip) = 0;
  virtual Bool_t TestSkipBit() const = 0;

 //TPC cluster occupancy
  virtual Int_t GetNclustersTPC(UInt_t /*sector*/) const = 0;
  virtual Int_t GetNclustersTPCused(UInt_t /*sector*/) const = 0;

private: 

  AliVfriendEvent(const AliVfriendEvent &);
  AliVfriendEvent& operator=(const AliVfriendEvent& esd);
};

#pragma GCC diagnostic ignored "-Weffc++" 
inline AliVfriendEvent::AliVfriendEvent(AliVConstructorReinitialisationFlag ) {} // do nothing
#pragma GCC diagnostic warning "-Weffc++" 

#endif

