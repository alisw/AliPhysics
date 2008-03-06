#ifndef ALIMUONPAINTERREGISTRY_H
#define ALIMUONPAINTERREGISTRY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterRegistry
/// \brief Registry for a bunch of AliMUONVPainter related stuff
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TQObject
#  include <TQObject.h>
#endif

class TGPopupMenu;
class TObjArray;
class AliMUONVTrackerData;
class AliMUONVTrackerDataMaker;
class AliMUONPainterMatrix;
class TGMenuBar;

class AliMUONPainterRegistry : public TObject, public TQObject
{
public:
  virtual ~AliMUONPainterRegistry();

  AliMUONVTrackerDataMaker* DataReader(Int_t i) const;

  AliMUONVTrackerData* DataSource(Int_t i) const;
  
  void DataSourceWasRegistered(AliMUONVTrackerData* data); // *SIGNAL*
  
  void DataSourceWasUnregistered(AliMUONVTrackerData* data); // *SIGNAL*

  void DataReaderWasRegistered(AliMUONVTrackerDataMaker* reader); // *SIGNAL*
  
  void DataReaderWasUnregistered(AliMUONVTrackerDataMaker* reader); // *SIGNAL*
  
  AliMUONVTrackerData* FindDataSource(const char* name) const;

  Int_t FindIndexOf(AliMUONPainterMatrix* group) const;
  
  Int_t FindIndexOf(AliMUONVTrackerData* data) const;
  
  AliMUONPainterMatrix* FindPainterMatrix(const char* groupName) const;
  
  void HistoryMenuActivated(Int_t i);
  
  static AliMUONPainterRegistry* Instance();
  
  Int_t NumberOfDataReaders() const;

  Int_t NumberOfDataSources() const;

  Int_t NumberOfPainterMatrices() const;

  AliMUONPainterMatrix* PainterMatrix(Int_t i) const;
  
  void AddToHistory(AliMUONPainterMatrix* group);
    
  void PainterMatrixWasRegistered(AliMUONPainterMatrix* group); // *SIGNAL*
  
  void PainterMatrixWasUnregistered(AliMUONPainterMatrix* group); // *SIGNAL*
  
  void PainterMatrixWantToShow(AliMUONPainterMatrix* group); // *SIGNAL*
  
  void Print(Option_t* opt) const;
  
  Int_t Register(AliMUONPainterMatrix* group);
  
  void Register(AliMUONVTrackerData* data);

  void Register(AliMUONVTrackerDataMaker* reader);

  /// Set the menu bar where to put the history menu
  void SetMenuBar(TGMenuBar* bar) { fMenuBar = bar; }
  
  Bool_t Unregister(AliMUONPainterMatrix* group);
  
  Bool_t Unregister(AliMUONVTrackerData* data);

  Bool_t Unregister(AliMUONVTrackerDataMaker* reader);

private:
  /// Not implemented
  AliMUONPainterRegistry();
  /// Not implemented
  AliMUONPainterRegistry(const AliMUONPainterRegistry&);
  /// Not implemented
  AliMUONPainterRegistry& operator=(const AliMUONPainterRegistry&);
  
private:
  static AliMUONPainterRegistry* fgInstance; ///< unique instance
  TObjArray* fDataSources; ///< data sources
  TObjArray* fPainterMatrices; ///< painter matrices
  TObjArray* fDataReaders; ///< data readers
  TGPopupMenu* fHistoryMenu; ///< history menu
  TGMenuBar* fMenuBar; ///< Menu bar where to put the history menu
  Int_t fHistoryCounter; ///< index to get back history menu items
  
  ClassDef(AliMUONPainterRegistry,1) // Registry for AliMUONVPainter related stuff
};

#endif
