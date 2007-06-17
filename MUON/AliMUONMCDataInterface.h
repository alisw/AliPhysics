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

class AliStack;
class AliMUONDataManager;
class AliMUONVHitStore;
class TClonesArray;

class AliMUONMCDataInterface : public TObject
{
public:
  AliMUONMCDataInterface(const char* filename="galice.root");
  virtual ~AliMUONMCDataInterface();
  
  AliMUONVHitStore* HitStore(Int_t event, Int_t track) const;
  void DumpHits(Int_t event) const;
  
  Bool_t IsValid() const;
  
  Int_t NumberOfEvents() const;
  
  Int_t NumberOfTracks(Int_t event) const;

  Int_t NumberOfTrackRefs(Int_t event) const;

  AliStack* Stack(Int_t event) const;
  void DumpKine(Int_t event) const;
  
  TClonesArray* TrackRefs(Int_t event, Int_t track) const;
  void DumpTrackRefs(Int_t event) const;
  
private:
  /// Not implemented
  AliMUONMCDataInterface(const AliMUONMCDataInterface&);
  /// Not implemented
  AliMUONMCDataInterface& operator=(const AliMUONMCDataInterface&);

private:
  
  AliMUONDataManager* fDataManager; //!< internal data accessor
  
  ClassDef(AliMUONMCDataInterface,0) // Easy to use MC data accessor
};

#endif
