#ifndef ALIMUONDATAINTERFACE_H
#define ALIMUONDATAINTERFACE_H
/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Includes revised 07/05/2004
//
/// \ingroup evaluation
/// \class AliMUONDataInterface
/// \brief An easy to use interface to MUON data

// Author: Artur Szostak (University of Cape Town)
//  email: artursz@iafrica.com
//
// Updated to MUON module w/o MUONData by Laurent Aphecetche, Subatech
//

#include <TObject.h>

class TIterator;
class AliLoader;
class AliMUONVStore;
class AliMUONVDigitStore;
class AliMUONVClusterStore;
class AliMUONVTrackStore;
class AliMUONVTriggerStore;
class AliMUONVTriggerTrackStore;
class AliMUONVDigit;
class AliMUONVCluster;
class AliMUONTrack;
class AliMUONLocalTrigger;
class AliMUONRegionalTrigger;
class AliMUONGlobalTrigger;
class AliMUONTriggerTrack;


class AliMUONDataInterface : public TObject
{
public:
  
  AliMUONDataInterface(const char* filename="galice.root");
  virtual ~AliMUONDataInterface();
  
  /// Returns true if the data interface was able to open the root file correctly.
  Bool_t IsValid() const { return fIsValid; };

  void Open(const char* filename);

  Int_t NumberOfEvents() const;
  
  /// Returns the index number of the current event as used int GetEvent(Int_t).
  Int_t   CurrentEvent() const { return fCurrentEvent; }

  AliMUONVDigitStore* DigitStore(Int_t event);  
  AliMUONVClusterStore* ClusterStore(Int_t event);
  AliMUONVTrackStore* TrackStore(Int_t event);
  AliMUONVTriggerStore* TriggerStore(Int_t event, const char* treeLetter="R");
  AliMUONVTriggerTrackStore* TriggerTrackStore(Int_t event);

  /// Dump the clusters for a given event, sorted if so required
  void DumpClusters(Int_t event, Bool_t sorted=kTRUE)  { return DumpRecPoints(event,sorted); }
  void DumpRecPoints(Int_t event, Bool_t sorted=kTRUE);
  void DumpDigits(Int_t event, Bool_t sorted=kTRUE);
  void DumpTracks(Int_t event, Bool_t sorted=kFALSE);  
  void DumpTrigger(Int_t event, const char* treeLetter="R");  
  void DumpTriggerTracks(Int_t event, Bool_t sorted=kFALSE);
  
  Bool_t GetEvent(Int_t event = 0);
  
  // Note the following methods can be extremely slow. Remember they are only
  // here for end user convenience for his/her small tests and macros.
  // If you want speed then don't use these methods. If you really want peak
  // performance then you should be talking to the AliRunLoader and Store
  // objects directly.
  Int_t NumberOfDigits(Int_t detElemId);
  AliMUONVDigit* Digit(Int_t detElemId, Int_t index);
  Int_t NumberOfDigits(Int_t chamber, Int_t cathode);
  AliMUONVDigit* Digit(Int_t chamber, Int_t cathode, Int_t index);
  Int_t NumberOfRawClusters(Int_t chamber);
  AliMUONVCluster* RawCluster(Int_t chamber, Int_t index);
  Int_t NumberOfTracks();
  AliMUONTrack* Track(Int_t index);
  Int_t NumberOfLocalTriggers();
  AliMUONLocalTrigger* LocalTrigger(Int_t index);
  Int_t NumberOfRegionalTriggers();
  AliMUONRegionalTrigger* RegionalTrigger(Int_t index);
  AliMUONGlobalTrigger* GlobalTrigger();
  Int_t NumberOfTriggerTracks();
  AliMUONTriggerTrack* TriggerTrack(Int_t index);
  
private:

  /// The various identifiers for the type of iterator constructed.
  enum IteratorType
  {
    kNoIterator,  ///< No iterator was constructed.
    kDigitIteratorByDetectorElement,  ///< A digit iterator for iterating over detector elements.
    kDigitIteratorByChamberAndCathode,  ///< A digit iterator for iterating over chambers and cathodes.
    kRawClusterIterator,  ///< A raw cluster iterator.
    kTrackIterator,  ///< An iterator for iterating over reconstructed tracks.
    kLocalTriggerIterator,  ///< An iterator for iterating over reconstructed local triggers.
    kRegionalTriggerIterator,  ///< An iterator for iterating over reconstructed regional triggers.
    kTriggerTrackIterator  ///< An iterator for iterating over reconstructed trigger tracks.
  };
    
  void DumpSorted(const AliMUONVStore& store) const;

  Bool_t LoadEvent(Int_t event);

  void NtupleTrigger(const char* treeLetter);
  
  void ResetStores();
  
  TIterator* GetIterator(IteratorType type, Int_t x = 0, Int_t y = 0);
  void ResetIterator();
  
  Int_t CountObjects(TIterator* iter);
  TObject* FetchObject(TIterator* iter, Int_t index);
  
  /// Not implemented
  AliMUONDataInterface(const AliMUONDataInterface& rhs);
  /// Not implemented
  AliMUONDataInterface& operator=(const AliMUONDataInterface& rhs);
  
  
  AliLoader* fLoader; //!< Tree accessor
  AliMUONVDigitStore* fDigitStore; //!< current digit store (owner)
  AliMUONVTriggerStore* fTriggerStore; //!< current trigger store (owner)
  AliMUONVClusterStore* fClusterStore; //!< current cluster store (owner)
  AliMUONVTrackStore* fTrackStore; //!< current track store (owner)
  AliMUONVTriggerTrackStore* fTriggerTrackStore; //!< current trigger track store (owner)
  Int_t fCurrentEvent; //!< Current event we've read in
  Bool_t fIsValid; //!< whether we were initialized properly or not
  
  IteratorType fCurrentIteratorType;  //!< The type of iterator that is currently set.
  Int_t fCurrentIndex;  //!< A current index number maintained for certain iteration operations.
  Int_t fDataX; //!< Extra data parameter about the iterator, can be the chamber number or detector element.
  Int_t fDataY; //!< Extra data parameter about the iterator, can be the cathode number.
  TIterator* fIterator; //!< Iterator for various iteration operations.
  
  static Int_t fgInstanceCounter; //!< To build unique folder name for each instance
  
  ClassDef(AliMUONDataInterface, 0)  // An easy to use interface to MUON reconstructed data
};
    

#endif // ALIMUONDATAINTERFACE_H
