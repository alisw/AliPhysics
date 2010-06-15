#ifndef ALIMUONPAINTERDATAREGISTRY_H
#define ALIMUONPAINTERDATAREGISTRY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id: AliMUONPainterDataRegistry.h 26812 2008-06-20 15:22:59Z laphecet $

/// \ingroup graphics
/// \class AliMUONPainterDataRegistry
/// \brief Registry for painter data sources
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TQObject
#  include <TQObject.h>
#endif

class TObjArray;
class AliMUONVTrackerData;
class AliMUONVTrackerDataMaker;

class AliMUONPainterDataRegistry : public TObject, public TQObject
{
public:
  virtual ~AliMUONPainterDataRegistry();

  AliMUONVTrackerDataMaker* DataMaker(Int_t i) const;

  AliMUONVTrackerData* DataSource(Int_t i) const;
  
  AliMUONVTrackerData* DataSource(const char* name) const;

  AliMUONVTrackerData* InteractiveReadOutConfig() const;
  
  void DataSourceWasRegistered(const AliMUONVTrackerData* data); // *SIGNAL*
  
  void DataSourceWasUnregistered(const AliMUONVTrackerData* data); // *SIGNAL*

  void DataMakerWasRegistered(const AliMUONVTrackerDataMaker* reader); // *SIGNAL*
  
  void DataMakerWasUnregistered(const AliMUONVTrackerDataMaker* reader); // *SIGNAL*
    
  static AliMUONPainterDataRegistry* Instance();
  
  Int_t NumberOfDataMakers() const;

  /// Number of data sources = data makers
  Int_t NumberOfDataSources() const { return NumberOfDataMakers(); }

  void Print(Option_t* opt) const;
  
  void Register(AliMUONVTrackerDataMaker* reader);

  Bool_t Unregister(AliMUONVTrackerDataMaker* reader);

  void DeleteZombies();
  
private:
  /// Not implemented
  AliMUONPainterDataRegistry();
  /// Not implemented
  AliMUONPainterDataRegistry(const AliMUONPainterDataRegistry&);
  /// Not implemented
  AliMUONPainterDataRegistry& operator=(const AliMUONPainterDataRegistry&);
  
  void CreateInteractiveReadOutConfig() const;
  
private:
  static AliMUONPainterDataRegistry* fgInstance; ///< unique instance
  TObjArray* fDataMakers; ///< data makers
  TObjArray* fZombies; ///< data readers to be deleted
  mutable AliMUONVTrackerData* fInteractiveReadOutConfig; ///< clickable readout configuration
  
  ClassDef(AliMUONPainterDataRegistry,1) // Registry for AliMUONVTrackerDataMaker objects
};

#endif
