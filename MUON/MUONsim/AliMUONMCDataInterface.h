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
//
// Moved parts of old AliMUONDataInterface interface to AliMUONMCDataInterface
//  Artur Szostak <artursz@iafrica.com> (University of Cape Town)

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVStore;
class AliMUONVHitStore;
class AliMUONVDigitStore;
class AliMUONVTriggerStore;
class AliMUONHit;
class AliMUONVDigit;
class AliMUONLocalTrigger;
class AliMUONRegionalTrigger;
class AliMUONGlobalTrigger;

class AliLoader;
class AliStack;

class TIterator;
class TClonesArray;
class TParticle;

class AliMUONMCDataInterface : public TObject
{
public:
  AliMUONMCDataInterface(const char* filename="galice.root");
  virtual ~AliMUONMCDataInterface();

  void Open(const char* filename);
  
  /// Returns true if the data interface was able to open the root file correctly.
  Bool_t IsValid() const { return fIsValid; };

  Int_t NumberOfEvents() const;

  /// Returns the index number of the current event loaded.
  /// This is the event number as was used in the last calls to any of the methods
  /// in this interface that have 'Int_t event' in the parameter list.
  /// GetEvent(Int_t event) for example.
  Int_t   CurrentEvent() const { return fCurrentEvent; }

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
  
  Bool_t GetEvent(Int_t event = 0);
  
  // Note the following methods can be extremely slow. Remember they are only
  // here for end user convenience for his/her small tests and macros.
  // If you want speed then don't use these methods. If you really want peak
  // performance then you should be talking to the AliRunLoader and Store
  // objects directly.
  Int_t NumberOfParticles();
  TParticle* Particle(Int_t index);
  Int_t NumberOfTracks();
  Int_t NumberOfHits(Int_t track);
  AliMUONHit* Hit(Int_t track, Int_t index);
  Int_t NumberOfSDigits(Int_t detElemId);
  AliMUONVDigit* SDigit(Int_t detElemId, Int_t index);
  Int_t NumberOfSDigits(Int_t chamber, Int_t cathode);
  AliMUONVDigit* SDigit(Int_t chamber, Int_t cathode, Int_t index);
  Int_t NumberOfDigits(Int_t detElemId);
  AliMUONVDigit* Digit(Int_t detElemId, Int_t index);
  Int_t NumberOfDigits(Int_t chamber, Int_t cathode);
  AliMUONVDigit* Digit(Int_t chamber, Int_t cathode, Int_t index);
  Int_t NumberOfLocalTriggers();
  AliMUONLocalTrigger* LocalTrigger(Int_t index);
  Int_t NumberOfRegionalTriggers();
  AliMUONRegionalTrigger* RegionalTrigger(Int_t index);
  AliMUONGlobalTrigger* GlobalTrigger();
  Int_t NumberOfTrackRefs();
  TClonesArray* TrackRefs(Int_t track);
  
private:

  /// The various identifiers for the type of iterator constructed.
  enum IteratorType
  {
    kNoIterator,  ///< No iterator was constructed.
    kHitIterator,  ///< An iterator to iterate over the hits.
    kSDigitIteratorByDetectorElement,  ///< A summable digit iterator to iterate over the detector elements.
    kSDigitIteratorByChamberAndCathode,  ///< A summable digit iterator to iterate over chambers and cathodes.
    kDigitIteratorByDetectorElement,  ///< An iterator for simulated digits to iterate over the detector elements.
    kDigitIteratorByChamberAndCathode,  ///< An iterator for simulated digits to iterate over chambers and cathodes.
    kLocalTriggerIterator,  ///< An iterator for iterating over the simulated local triggers.
    kRegionalTriggerIterator  ///< An iterator for iterating over the simulated regional triggers.
  };
  
  /// Not implemented
  AliMUONMCDataInterface(const AliMUONMCDataInterface&);
  /// Not implemented
  AliMUONMCDataInterface& operator=(const AliMUONMCDataInterface&);

  void DumpSorted(const AliMUONVStore& store) const;
  Bool_t LoadEvent(Int_t event);
  
  void ResetStores();
  
  TIterator* GetIterator(IteratorType type, Int_t x = 0, Int_t y = 0);
  void ResetIterator();
  
  Int_t CountObjects(TIterator* iter);
  TObject* FetchObject(TIterator* iter, Int_t index);
  
  
  AliLoader* fLoader; //!<! Tree accessor
  AliMUONVHitStore* fHitStore; //!<! current hit store (owner)
  AliMUONVDigitStore* fSDigitStore; //!<! current sdigit store (owner)
  AliMUONVDigitStore* fDigitStore; //!<! current digit store (owner)
  AliMUONVTriggerStore* fTriggerStore; //!<! current trigger store (owner)
  TClonesArray* fTrackRefs; //!<! current trackrefs (owner)
  Int_t fCurrentEvent; //!<! Current event we've read in
  Bool_t fIsValid; //!<! whether we were initialized properly or not
  
  IteratorType fCurrentIteratorType;  //!<! The type of iterator that is currently set.
  Int_t fCurrentIndex;  //!<! A current index number maintained for certain iteration operations.
  Int_t fDataX; //!<! Extra data parameter about the iterator, can be the chamber number, detector element or track number.
  Int_t fDataY; //!<! Extra data parameter about the iterator, can be the cathode number.
  TIterator* fIterator; //!<! Iterator for various iteration operations.
  
  static Int_t fgInstanceCounter; //!<! To build unique folder name for each instance

  ClassDef(AliMUONMCDataInterface,0) // Easy to use MC data accessor
};

#endif
