#ifndef ALIMUONTREEMANAGER_H
#define ALIMUONTREEMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUONTreeManager
/// \brief Helper class to ease TTree (MUON) branches manipulations
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVStore;
class TTree;
class TBranch;

class AliMUONTreeManager : public TObject
{
public:

  AliMUONTreeManager();
  virtual ~AliMUONTreeManager();

  void AddClassName(TTree& tree, const char* pattern, 
                    const char* className) const;
    
  Bool_t MakeBranch(TTree& tree, const char* storeClassName,
                    const char* branchClassName, const char* branchName, 
                    void* address,
                    Int_t bufferSize = 4000, Int_t splitLevel = 99) const;
  
  Bool_t SetAddress(TTree& tree, const char* branchName, void* address) const;
  
  TObject* CreateObject(const TTree& tree, const char* detail) const;
  
  void UpdateBranchStatuses(TTree& tree, const char* pattern) const;
  
  /** Debug method to get an event, but checking beforehand that all selected
    branches do have a non-zero address set (otherwise we leak memory). 
    */
  void GetEvent(TTree& tree, Int_t event) const;
  
  /// Debug method to show the tree branch statuses and addresses.
  void ShowStatus(TTree& tree) const;

private:

  TString GetClassName(const TTree& tree, const char* pattern,
                       Bool_t makeDefault) const;
    
  TString DefaultClassName(const char* treename, const char* pattern) const;
  
  ClassDef(AliMUONTreeManager,0) // Helper class to handle MUON TTrees
};

#endif
