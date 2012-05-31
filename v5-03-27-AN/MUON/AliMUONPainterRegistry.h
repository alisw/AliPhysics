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
class AliMUONPainterMatrix;
class TGMenuBar;

class AliMUONPainterRegistry : public TObject, public TQObject
{
public:
  virtual ~AliMUONPainterRegistry();

  Int_t FindIndexOf(AliMUONPainterMatrix* group) const;
  
  void HistoryMenuActivated(Int_t i);
  
  static AliMUONPainterRegistry* Instance();
  
  Int_t NumberOfPainterMatrices() const;

  AliMUONPainterMatrix* PainterMatrix(Int_t i) const;
  
  AliMUONPainterMatrix* PainterMatrix(const char* groupName) const;
  
  void AddToHistory(AliMUONPainterMatrix* group);
    
  void PainterMatrixWasRegistered(const AliMUONPainterMatrix* group); // *SIGNAL*
  
  void PainterMatrixWasUnregistered(const AliMUONPainterMatrix* group); // *SIGNAL*
  
  void PainterMatrixWantToShow(const AliMUONPainterMatrix* group); // *SIGNAL*
  
  void Print(Option_t* opt) const;
  
  Int_t Register(AliMUONPainterMatrix* group);
  
  /// Set the menu bar where to put the history menu
  void SetMenuBar(TGMenuBar* bar) { fMenuBar = bar; }
  
  Bool_t Unregister(AliMUONPainterMatrix* group);
  
private:
  /// Not implemented
  AliMUONPainterRegistry();
  /// Not implemented
  AliMUONPainterRegistry(const AliMUONPainterRegistry&);
  /// Not implemented
  AliMUONPainterRegistry& operator=(const AliMUONPainterRegistry&);
  
private:
  static AliMUONPainterRegistry* fgInstance; ///< unique instance
  TObjArray* fPainterMatrices; ///< painter matrices
  TGPopupMenu* fHistoryMenu; ///< history menu
  TGMenuBar* fMenuBar; ///< Menu bar where to put the history menu
  Int_t fHistoryCounter; ///< index to get back history menu items
  
  ClassDef(AliMUONPainterRegistry,5) // Registry for AliMUONVPainter related stuff
};

#endif
