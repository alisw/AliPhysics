#ifndef ALI_MUON_DATA_INTERFACE_H
#define ALI_MUON_DATA_INTERFACE_H
/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Includes revised 07/05/2004
//
/// \ingroup evaluation
/// \class AliMUONDataInterface
/// \brief An easy to use interface to MUON data

// Author: Artur Szostak
//  email: artur@alice.phy.uct.ac.za
//
// Updated to MUON module w/o MUONData by Laurent Aphecetche, Subatech
//

#include <TObject.h>

class AliMUONDataManager;
class AliMUONVDigitStore;
class AliMUONVTriggerStore;
class AliMUONVTriggerTrackStore;
class AliMUONVClusterStore;
class AliMUONVTrackStore;

// >>LA should not needed (once we remove deprecated methods)
class AliMUONHit; 
class AliMUONDigit;
class AliMUONLocalTrigger;
class AliMUONGlobalTrigger;
class AliMUONRawCluster;
class AliMUONTrack;
class TParticle;
// <<LA
class TList;

class AliMUONDataInterface : public TObject
{
 public:
  
  AliMUONDataInterface(const char* filename="galice.root");
  virtual ~AliMUONDataInterface();
  
  Bool_t IsValid() const;

  AliMUONVClusterStore* ClusterStore(Int_t event) const;
  void DumpRecPoints(Int_t event, Bool_t sorted=kFALSE) const;

  /** Return the digit store for one event. The returned pointer should be
    deleted by the client.
    */
  AliMUONVDigitStore* DigitStore(Int_t event) const;  
  void DumpDigits(Int_t event, Bool_t sorted=kFALSE) const;
  TList* DigitStoreAsList(Int_t event) const;
  
  Int_t NumberOfEvents() const;

  AliMUONVDigitStore* SDigitStore(Int_t event) const;
  void DumpSDigits(Int_t event, Bool_t sorted=kFALSE) const;
  
  AliMUONVTrackStore* TrackStore(Int_t event) const;
  void DumpTracks(Int_t event, Bool_t sorted=kFALSE) const;
  
  /// Get the triggerStore from the given tree (can be "D" or "R").
  AliMUONVTriggerStore* TriggerStore(Int_t event, const char* treeLetter) const;
  void DumpTrigger(Int_t event, const char* treeLetter="R") const;
  
  AliMUONVTriggerTrackStore* TriggerTrackStore(Int_t event) const;
  void DumpTriggerTracks(Int_t event, Bool_t sorted=kFALSE) const;
  
  // all the methods below are deprecated.
  
  // Sets all internal pointers to NULL without releasing the current runloader.
  void Reset();
  
  Bool_t UseCurrentRunLoader();
  
  Int_t NumberOfEvents(TString filename, TString foldername);
  
  Int_t NumberOfParticles(TString filename, TString foldername, Int_t event);
  TParticle* Particle(TString filename, TString foldername, Int_t event, Int_t particle);
  
  Int_t NumberOfTracks(TString filename, TString foldername, Int_t event);
  Int_t NumberOfHits(TString filename, TString foldername, Int_t event, Int_t track);
  AliMUONHit* Hit(TString filename, TString foldername, Int_t event, Int_t track, Int_t hit);
  
  Int_t NumberOfSDigits(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cathode);
  AliMUONDigit* SDigit(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cathode, Int_t sdigit);
  
  Int_t NumberOfDigits(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cathode);
  AliMUONDigit* Digit(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cathode, Int_t digit);
  
  Int_t NumberOfRawClusters(TString filename, TString foldername, Int_t event, Int_t chamber);
  AliMUONRawCluster* RawCluster(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cluster);
  
  Int_t NumberOfLocalTriggers(TString filename, TString foldername, Int_t event);
  AliMUONLocalTrigger* LocalTrigger(TString filename, TString foldername, Int_t event, Int_t trigger);
  
  Bool_t SetFile(TString filename = "galice.root", TString foldername = "MUONFolder");
  Bool_t GetEvent(Int_t event = 0);
  
  Int_t NumberOfParticles();
  TParticle* Particle(Int_t particle);
  
  Int_t NumberOfTracks();
  Int_t NumberOfHits(Int_t track);
  AliMUONHit* Hit(Int_t track, Int_t hit);
  
  Int_t NumberOfSDigits(Int_t chamber, Int_t cathode);
  AliMUONDigit* SDigit(Int_t chamber, Int_t cathode, Int_t sdigit);
  
  Int_t NumberOfDigits(Int_t chamber, Int_t cathode);
  AliMUONDigit* Digit(Int_t chamber, Int_t cathode, Int_t digit);
  
  Int_t NumberOfRawClusters(Int_t chamber);
  AliMUONRawCluster* RawCluster(Int_t chamber, Int_t cluster);
  
  Int_t NumberOfLocalTriggers();
  AliMUONLocalTrigger* LocalTrigger(Int_t trigger);
  
  Int_t NumberOfGlobalTriggers();
  AliMUONGlobalTrigger* GlobalTrigger(Int_t trigger);
  // Returns the name of the currently selected file.
  
  Int_t NumberOfRecTracks();
  AliMUONTrack* RecTrack(Int_t rectrack);
  
  /// Returns the file name from which we are fetching data
  TString CurrentFile() const { return ""; }
  
  /// Returns the name of the currently selected folder.
  TString CurrentFolder() const { return ""; }
  
  /// Returns the number of the currently selected event.
  Int_t   CurrentEvent() const { return 0; }
  
  /// Returns the currently selected track.
  Int_t   CurrentTrack() const { return 0; }
  
  /// Returns the currently selected cathode in TreeS.
  Int_t   CurrentSCathode() const { return 0; }
  
  /// Returns the currently selected cathode in TreeD.
  Int_t   CurrentDCathode() const { return 0; }
  
 private:
    
  void DumpIt(const char* treeLetter, const char* what, Int_t event, Bool_t sorted) const;
  
  /// Not implemented
  AliMUONDataInterface(const AliMUONDataInterface& rhs);
  /// Not implemented
  AliMUONDataInterface& operator=(const AliMUONDataInterface& rhs);

  AliMUONDataManager* fDataManager; //!< Internal data accessor
  
  ClassDef(AliMUONDataInterface, 0)  // An easy to use interface to MUON data
};
    

#endif // ALI_MUON_DATA_INTERFACE_H
