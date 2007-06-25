#ifndef ALIMUONMCDATAINTERFACE_H
#define ALIMUONMCDATAINTERFACE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONMCDataInterface
/// \brief Easy to use data access to MC information
/// 
// Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliLoader;
class AliMUONVDigitStore;
class AliMUONVHitStore;
class AliMUONVStore;
class AliMUONVTriggerStore;
class AliStack;
class TClonesArray;

class AliMUONMCDataInterface : public TObject
{
public:
  AliMUONMCDataInterface(const char* filename="galice.root");
  virtual ~AliMUONMCDataInterface();

  void Open(const char* filename);
  
  Bool_t IsValid() const;

  Int_t NumberOfEvents() const;
  Int_t NumberOfTracks(Int_t event);
  Int_t NumberOfTrackRefs(Int_t event);
  
  AliMUONVHitStore* HitStore(Int_t event, Int_t track);
  AliMUONVDigitStore* SDigitStore(Int_t event);
  AliMUONVDigitStore* DigitStore(Int_t event);
  AliStack* Stack(Int_t event);  
  TClonesArray* TrackRefs(Int_t event, Int_t track);
  AliMUONVTriggerStore* TriggerStore(Int_t event);
  
  void DumpDigits(Int_t event, Bool_t sorted=kTRUE);
  void DumpSDigits(Int_t event, Bool_t sorted=kTRUE);
  void DumpHits(Int_t event);
  void DumpKine(Int_t event);
  void DumpTrackRefs(Int_t event);
  void DumpTrigger(Int_t event);
  
private:
  /// Not implemented
  AliMUONMCDataInterface(const AliMUONMCDataInterface&);
  /// Not implemented
  AliMUONMCDataInterface& operator=(const AliMUONMCDataInterface&);

  void DumpSorted(const AliMUONVStore& store) const;
  Int_t LoadEvent(Int_t event);
  
private:
  
  AliLoader* fLoader; //!< Tree accessor
  AliMUONVHitStore* fHitStore; //!< current hit store (owner)
  AliMUONVDigitStore* fSDigitStore; //!< current sdigit store (owner)
  AliMUONVDigitStore* fDigitStore; //!< current digit store (owner)
  AliMUONVTriggerStore* fTriggerStore; //!< current trigger store (owner)
  TClonesArray* fTrackRefs; //!< current trackrefs (owner)
  Int_t fCurrentEvent; //!< Current event we've read in
  Bool_t fIsValid; //!< whether we were initialized properly or not
  
  static Int_t fgInstanceCounter; //!< To build unique folder name for each instance

  ClassDef(AliMUONMCDataInterface,0) // Easy to use MC data accessor
};

#endif
