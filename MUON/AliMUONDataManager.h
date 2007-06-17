#ifndef ALIMUONDATAMANAGER_H
#define ALIMUONDATAMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONDataManager
/// \brief Utility class to ease access to data stores.
/// 
// Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliLoader;
class AliMUONVStore;
class TTree;

class AliMUONDataManager : public TObject
{
public:

  AliMUONDataManager(const char* file="galice.root");
  virtual ~AliMUONDataManager();

  AliMUONVStore* ReadConnectable(Int_t event, const char* tree, const char* what);

  /// Whether we were properly initialized or not
  Bool_t IsValid() const { return fIsValid; }
  
  Int_t NumberOfEvents() const;
  
  Int_t Load(Int_t event);
  void Load(const char* tree);
  TTree* Tree(const char* tree);
  void Unload(const char* tree);
  
private:
  /// Not implemented
  AliMUONDataManager(const AliMUONDataManager&);
  /// Not implemented
  AliMUONDataManager& operator=(const AliMUONDataManager&);
  
  AliLoader* fLoader; //!< Our loader to access trees
  Int_t fCurrentEvent; //!< Current loaded event
  Bool_t fIsValid; //!< Whether we were properly initialized or not
  static Int_t fgCount; //!< instance counter to be able to build a unique folder name
  
  ClassDef(AliMUONDataManager,0) // Utility class to ease data store access
};

#endif
